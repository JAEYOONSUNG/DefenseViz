"""Build a wheel, install into a clean venv, exercise CLI outside source tree.

Requires network access to obtain declared dependencies into the new environment.
No historical analysis directory or PYTHONPATH is used.
"""
import json
import os
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def main():
    env = os.environ.copy()
    for name in ('PYTHONPATH', 'PYTHONHOME'):
        env.pop(name, None)
    env['PYTHONNOUSERSITE'] = '1'
    def run(args, cwd):
        p = subprocess.run([str(x) for x in args], cwd=cwd, env=env, text=True, capture_output=True)
        if p.returncode:
            raise RuntimeError(f'{args}\n{p.stdout}\n{p.stderr}')
        return p.stdout
    with tempfile.TemporaryDirectory(prefix='defenseviz-install-') as temp:
        work = Path(temp)
        run([sys.executable, '-m', 'build', '--wheel', '--outdir', work/'wheels', ROOT], work)
        wheel = next((work/'wheels').glob('*.whl'))
        venv.EnvBuilder(with_pip=True).create(work/'venv')
        bindir = work/'venv'/('Scripts' if os.name=='nt' else 'bin')
        python = bindir/('python.exe' if os.name=='nt' else 'python')
        run([python, '-m', 'pip', 'install', '--disable-pip-version-check', wheel], work)
        assert 'DefenseViz 2.0.0' in run([python, '-m', 'defenseviz', 'version'], work)
        run([python, '-m', 'defenseviz', 'demo', '--out', work/'demo'], work)
        run([python, '-m', 'defenseviz', 'validate', work/'demo/inputs/project.json'], work)
        run([python, '-m', 'defenseviz', 'run', work/'demo/inputs/project.json', '--out', work/'rerun'], work)
        a=json.loads((work/'demo/results/analysis.json').read_text())
        b=json.loads((work/'rerun/analysis.json').read_text())
        assert a==b, 'Demo and standalone run must be analytically identical'
        assert (work/'rerun/figure2_ab.pdf').stat().st_size>1000
        assert list((work/'rerun').glob('*.svg'))
        receipt=json.loads((work/'rerun/receipt.json').read_text())
        assert receipt['detector_runs']==0 and receipt['exit_code']==0
        assert receipt['software']['git']['status']=='not_a_source_checkout'
        print('CLEAN_INSTALL_OK')

if __name__=='__main__':main()
