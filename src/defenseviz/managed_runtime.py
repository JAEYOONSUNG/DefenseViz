"""A user-owned Docker runtime and reference cache for first installations.

The public source checkout supplies the build recipe. Provider files, built
references and run receipts stay outside that checkout. No shell interpolation
or host R installation is needed; existing explicit native engines take priority.
"""
from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import uuid

PLATFORM = 'linux/amd64'
ENTRYPOINT = '/opt/defenseviz/bin/defenseviz'
COPY_FILES = ('Dockerfile.dnmbsuite', 'Dockerfile.dnmbsuite.dockerignore',
              'pyproject.toml', 'README.md', 'LICENSE',
              'scripts/install_native_runtime.py', 'scripts/prepare_hmm_resources.py')
COPY_TREES = ('src', 'docker/dnmbsuite')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n', encoding='utf-8')
    os.replace(temporary, path)


def cache_directory(explicit=None):
    value = explicit or os.environ.get('DEFENSEVIZ_CACHE')
    if value:
        return Path(value).expanduser().resolve()
    return (Path(os.environ.get('XDG_CACHE_HOME', Path.home() / '.cache')) / 'defenseviz').expanduser().resolve()


def _supported_host():
    if os.name != 'posix':
        raise RuntimeError('Managed Docker setup supports Linux/macOS and Windows through WSL. Run this command inside WSL on Windows.')


def _docker():
    executable = shutil.which('docker')
    if executable is None:
        raise RuntimeError('Docker is required for managed setup. Install and start Docker Desktop (macOS) or Docker Engine (Linux), then run defenseviz setup again.')
    return executable


def _command(argv, record_dir, name, *, cwd=None):
    """Retain exact command, binary logs, elapsed time and exit status."""
    record_dir = Path(record_dir)
    started = time.monotonic()
    receipt = dict(argv=[str(v) for v in argv], cwd=str(Path(cwd or Path.cwd()).resolve()),
                   started_utc=datetime.now(timezone.utc).isoformat(), completed=False)
    out = record_dir / (name + '.stdout')
    err = record_dir / (name + '.stderr')
    try:
        with out.open('wb') as stdout, err.open('wb') as stderr:
            result = subprocess.run(receipt['argv'], cwd=cwd, stdout=stdout, stderr=stderr)
        receipt.update(completed=True, exit_code=result.returncode,
                       stdout_sha256=digest(out), stderr_sha256=digest(err))
        if result.returncode:
            detail = err.read_bytes()[-3000:].decode('utf-8', errors='replace')
            raise RuntimeError(f'{name} failed (exit {result.returncode}); logs: {record_dir}\n{detail}')
        return out.read_text(encoding='utf-8', errors='replace')
    except OSError as exc:
        receipt['error'] = str(exc)
        raise
    finally:
        receipt.update(elapsed_seconds=time.monotonic()-started,
                       ended_utc=datetime.now(timezone.utc).isoformat())
        write_json(record_dir / (name + '.json'), receipt)


@contextmanager
def _setup_lock(cache):
    # Kernel locks release on interruption; a killed setup cannot strand a lock.
    import fcntl
    with (cache / '.setup.lock').open('a+b') as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError('Another DefenseViz setup is using this cache; wait for it to finish.') from exc
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def source_state(source):
    source = Path(source).expanduser().resolve()
    selected = [source / n for n in COPY_FILES]
    for tree in COPY_TREES:
        selected.extend(p for p in (source / tree).rglob('*')
                        if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'
                        and not any(part.endswith('.egg-info') for part in p.parts))
    if any(not p.is_file() for p in selected[:len(COPY_FILES)]):
        raise ValueError('Source must be a complete DefenseViz checkout containing Dockerfile.dnmbsuite. Run defenseviz setup --source /path/to/DefenseViz.')
    if any(p.is_symlink() for p in selected):
        raise ValueError('Managed runtime build inputs must be regular files, not symbolic links')
    hashes = {p.relative_to(source).as_posix(): digest(p) for p in sorted(selected)}
    # The caller and the image must run the same package implementation.
    installed = Path(__file__).resolve().parent
    for relative, expected in hashes.items():
        if not relative.startswith('src/defenseviz/'):
            continue
        counterpart = installed / relative.removeprefix('src/defenseviz/')
        if not counterpart.is_file() or digest(counterpart) != expected:
            raise ValueError('Installed DefenseViz differs from --source; reinstall this checkout with python -m pip install . before setup')
    content = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
    commit = None
    status = None
    if shutil.which('git') and (source / '.git').exists():
        result = subprocess.run(['git', '-C', str(source), 'rev-parse', 'HEAD'], capture_output=True, text=True)
        if result.returncode == 0:
            commit = result.stdout.strip()
            state = subprocess.run(['git', '-C', str(source), 'status', '--porcelain=v1', '--untracked-files=all'], capture_output=True, text=True)
            if state.returncode:
                raise RuntimeError('Cannot record source working-tree state')
            status = state.stdout
    return dict(path=str(source), commit=commit, working_tree=status, files=hashes,
                content_sha256=content, revision=commit or ('source-' + content))


def _inventory(root):
    root = Path(root)
    return {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob('*')) if p.is_file()}


def _installed_files():
    root = Path(__file__).resolve().parent
    return {p.relative_to(root).as_posix(): digest(p) for p in sorted(root.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts and p.suffix != '.pyc'}


def _verify_installed_package(profile):
    expected = {p.removeprefix('src/defenseviz/'): value
                for p, value in profile.get('source', {}).get('files', {}).items()
                if p.startswith('src/defenseviz/')}
    if not expected or expected != _installed_files():
        raise ValueError('Installed DefenseViz differs from its managed image. Run defenseviz setup --source /path/to/the/current/DefenseViz checkout before analysis.')


def _verify_profile(profile, cache):
    if profile.get('schema_version') != 1 or profile.get('ready') is not True or profile.get('kind') != 'docker':
        raise ValueError('Managed runtime profile is incomplete; run defenseviz setup')
    import re
    if not re.fullmatch(r'sha256:[0-9a-f]{64}', profile.get('image_id', '')):
        raise ValueError('Managed runtime has an invalid image identity')
    resources = Path(profile['resources']).resolve()
    if not resources.is_relative_to(cache) or not resources.is_dir():
        raise ValueError('Managed reference directory is missing or outside its cache; run defenseviz setup')
    hashes = profile.get('resource_sha256', {})
    if not hashes or _inventory(resources) != hashes:
        raise ValueError('Managed database files changed or are incomplete; inspect the cache and run defenseviz setup --refresh-databases')
    engine = json.loads((resources / 'runtime-engine.json').read_text())
    for key, expected in {'reference_cache': '/resources/rebase',
                          'hmm_database': '/resources/rm-hmm/RM_HMMs.hmm',
                          'profile_metadata': '/resources/rm-hmm/profile_roles.tsv'}.items():
        if engine.get(key) != expected:
            raise ValueError('Managed database configuration does not match its mount layout')
    return profile


def load_profile(cache=None):
    cache = cache_directory(cache)
    path = cache / 'runtime.json'
    if not path.exists():
        return None
    data = path.read_bytes()
    profile = _verify_profile(json.loads(data), cache)
    profile['_profile_sha256'] = hashlib.sha256(data).hexdigest()
    return profile


def _mount(source, target, readonly=True):
    source = str(Path(source).resolve())
    # Docker's --mount grammar treats commas as separators, not path characters.
    if ',' in source or ',' in str(target) or '\n' in source:
        raise ValueError('Docker mount paths cannot contain commas or newlines')
    return ['--mount', f'type=bind,src={source},dst={target}' + (',readonly' if readonly else '')]


def _container(docker, image):
    return [docker, 'run', '--rm', '--platform', PLATFORM, '--network', 'none',
            '--read-only', '--tmpfs', '/tmp:rw,exec,size=1073741824',
            '--user', f'{os.getuid()}:{os.getgid()}',
            '--env', 'DEFENSEVIZ_MANAGED_CONTAINER=1', '--entrypoint', ENTRYPOINT]


def setup_managed(*, source=None, cache=None, refresh=False, cpu=4):
    _supported_host()
    if type(cpu) is not int or not 1 <= cpu <= 256:
        raise ValueError('cpu must be an integer from 1 to 256')
    cache = cache_directory(cache)
    cache.mkdir(parents=True, exist_ok=True)
    docker = _docker()
    source = source_state(source or Path.cwd())
    with _setup_lock(cache):
        attempt = cache / 'setups' / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + uuid.uuid4().hex[:8])
        attempt.mkdir(parents=True)
        receipt = dict(schema_version=1, ready=False, source=source, cache=str(cache), cpu=cpu,
                       started_utc=datetime.now(timezone.utc).isoformat())
        write_json(attempt / 'setup.json', receipt)
        try:
            info = json.loads(_command([docker, 'info', '--format', '{{json .}}'], attempt, 'docker_info'))
            if info.get('OSType') != 'linux':
                raise RuntimeError('DefenseViz requires the Linux Docker engine')
            previous = None if refresh else load_profile(cache)
            if previous and previous['source']['content_sha256'] == source['content_sha256'] and previous.get('cpu') == cpu:
                inspected = json.loads(_command([docker, 'image', 'inspect', previous['image_id']], attempt, 'image_inspect'))
                if inspected[0]['Id'] != previous['image_id']:
                    raise ValueError('Cached runtime image identity changed')
                receipt.update(ready=True, reused=True, provider_downloads=0, profile_sha256=digest(cache / 'runtime.json'))
                return dict(ready=True, reused=True, provider_downloads=0, cache=str(cache), setup_receipt=str(attempt / 'setup.json'))
            from .resource_download import fetch_provider_sources
            print('DefenseViz: downloading/verifying public reference sources...', file=sys.stderr, flush=True)
            provider_root = cache / 'provider-sources'
            provider = fetch_provider_sources(provider_root, refresh=refresh)
            receipt['providers'] = provider
            for key in ('rebase_source', 'hmm', 'hmm_metadata', 'off_target', 'receipt'):
                if not Path(provider[key]).resolve().is_relative_to(provider_root):
                    raise ValueError('Provider source is outside its managed cache')
            tag = 'defenseviz-managed:' + source['content_sha256'][:24]
            print('DefenseViz: building the public native runtime; logs: ' + str(attempt), file=sys.stderr, flush=True)
            _command([docker, 'build', '--platform', PLATFORM, '--build-arg', 'DEFENSEVIZ_REVISION=' + source['revision'],
                      '-f', 'Dockerfile.dnmbsuite', '-t', tag, '.'], attempt, 'image_build', cwd=source['path'])
            if source_state(source['path'])['content_sha256'] != source['content_sha256']:
                raise ValueError('Source changed during the runtime build; rerun setup from stable source')
            inspected = json.loads(_command([docker, 'image', 'inspect', tag], attempt, 'image_inspect'))[0]
            image = inspected['Id']
            receipt['image_id'] = image
            argv = _container(docker, image) + _mount(provider_root, '/providers') + _mount(attempt, '/setup', False)
            argv += [image, 'prepare-resources', '--engine-config', '/opt/defenseviz/container/engine.json', '--out', '/setup/databases']
            for option, key in (('--rebase-source', 'rebase_source'), ('--hmm', 'hmm'), ('--hmm-metadata', 'hmm_metadata'), ('--off-target', 'off_target')):
                argv += [option, '/providers/' + Path(provider[key]).resolve().relative_to(provider_root).as_posix()]
            argv += ['--source-label', provider['source_label']]
            print('DefenseViz: preparing and validating REBASE and R-M HMM indexes...', file=sys.stderr, flush=True)
            _command(argv, attempt, 'prepare_resources')
            databases = attempt / 'databases'
            engine = json.loads((databases / 'engine.json').read_text())
            engine.update(reference_cache='/resources/rebase', hmm_database='/resources/rm-hmm/RM_HMMs.hmm',
                          profile_metadata='/resources/rm-hmm/profile_roles.tsv', cpu=cpu)
            write_json(databases / 'runtime-engine.json', engine)
            argv = _container(docker, image) + _mount(databases, '/resources') + _mount(attempt, '/setup', False)
            argv += [image, 'check-resources', '--engine-config', '/resources/runtime-engine.json', '--out', '/setup/resource-check']
            _command(argv, attempt, 'check_resources')
            profile = dict(schema_version=1, ready=True, kind='docker', platform=PLATFORM, image_id=image,
                           source=source, resources=str(databases), resource_sha256=_inventory(databases),
                           providers=provider, cpu=cpu, setup_receipt=str(attempt / 'setup.json'))
            _verify_profile(profile, cache)
            write_json(cache / 'runtime.json', profile)
            receipt.update(ready=True, reused=False, profile_sha256=digest(cache / 'runtime.json'))
            return dict(ready=True, reused=False, cache=str(cache), image_id=image,
                        setup_receipt=str(attempt / 'setup.json'), next_command='defenseviz from-suite /path/to/suite-parent --out results')
        except Exception as exc:
            receipt['error'] = str(exc)
            raise
        finally:
            receipt['ended_utc'] = datetime.now(timezone.utc).isoformat()
            write_json(attempt / 'setup.json', receipt)


def managed_analysis(args, source_root, settings, project_refs):
    """Return None for native/stored workflows, otherwise run the saved image."""
    if os.environ.get('DEFENSEVIZ_MANAGED_CONTAINER') == '1' or settings.get('evaluation') == 'stored':
        return None
    if any(settings.get(name) is not None for name in ('rm_evidence', 'reference_cache', 'rscript')):
        # Explicit saved-evidence/native overrides retain their existing semantics.
        return None
    from .native_pipeline import engine_config_path
    if engine_config_path(settings.get('engine_config')) is not None:
        return None
    profile = load_profile()
    if profile is None:
        if settings.get('evaluation') == 'full':
            raise RuntimeError('No full engine is configured. Run defenseviz setup --source /path/to/DefenseViz first, or supply --engine-config.')
        return None
    _verify_installed_package(profile)
    _supported_host()
    if settings.get('external_type_projects'):
        raise ValueError('A managed run uses its cached runtime and references. Custom native executables, reference caches or external projects require an explicit native --engine-config.')
    docker = _docker()
    target = Path(args.out).expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise FileExistsError('Output must be a new directory: ' + str(target))
    target = target.resolve()
    bindings = {Path(source_root).resolve(), *(Path(r['path']).resolve() for r in project_refs)}
    root_arg = Path(args.root).resolve()
    bindings.add(root_arg)
    from .suite_project import PATH_OPTIONS
    for name in PATH_OPTIONS - {'engine_config', 'reference_cache'}:
        if settings.get(name) is not None:
            bindings.add(Path(settings[name]).resolve())
    for path in bindings:
        if not path.exists():
            raise FileNotFoundError('Missing managed input: ' + str(path))
        if target == path or target.is_relative_to(path) or path.is_relative_to(target):
            raise ValueError('Managed results must be outside the input directory')
    # Work in a separate directory so the analysis itself still sees a new out.
    target.parent.mkdir(parents=True, exist_ok=True)
    run = target.parent / ('.' + target.name + '.defenseviz-' + uuid.uuid4().hex[:8])
    run.mkdir()
    write_json(run / 'runtime_profile.json', profile)
    argv = _container(docker, profile['image_id'])
    for path in sorted(bindings):
        # Same absolute POSIX path preserves relative paths inside project files.
        if any(path != parent and parent.is_dir() and path.is_relative_to(parent) for parent in bindings):
            continue
        if str(path) in {'/', '/resources', '/output', '/tmp', '/opt'}:
            raise ValueError('Choose a dedicated input directory rather than a container system path')
        argv += _mount(path, str(path))
    argv += _mount(profile['resources'], '/resources') + _mount(run, '/output', False)
    argv += [profile['image_id'], 'from-suite', str(root_arg), '--out', '/output/results',
             '--engine-config', '/resources/runtime-engine.json', '--evaluation', 'full']
    for name in ('manifest', 'supplement', 'rm_evidence', 'system_review_rules', 'system_review_evidence', 'layout'):
        value = getattr(args, name, None)
        if value is not None:
            argv += ['--' + name.replace('_', '-'), str(Path(value).resolve())]
    if getattr(args, 'title', None) is not None:
        argv += ['--title', args.title]
    if args.strict_inputs is not None:
        argv.append('--strict-inputs' if args.strict_inputs else '--no-strict-inputs')
    print('DefenseViz: using cached databases; running fresh full analysis. Logs: ' + str(run), file=sys.stderr, flush=True)
    output = _command(argv, run, 'analysis')
    result = run / 'results'
    if not (result / 'receipt.json').is_file() or not (result / 'figures/publication.pdf').is_file():
        raise RuntimeError('Analysis did not produce its final receipt and figures; inspect ' + str(run))
    _verify_profile(profile, cache_directory())
    if target.exists() or target.is_symlink():
        raise FileExistsError('Output appeared during analysis; completed files remain at ' + str(result))
    result.rename(target)
    write_json(run / 'managed_run.json', dict(completed=True, output=str(target.resolve()),
               image_id=profile['image_id'], runtime_profile_sha256=profile['_profile_sha256'],
               result_receipt_sha256=digest(target / 'receipt.json'), provider_downloads=0,
               root=str(root_arg), database_inputs_unchanged=True))
    return dict(completed=True, output=str(target.resolve()), managed_runtime=True,
                receipt=str(target / 'receipt.json'), execution_receipt=str(run / 'analysis.json'),
                figures=str(target / 'figures/publication.pdf'), provider_downloads=0)
