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
import zipfile
from email.parser import Parser
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
        # Synthetic packaging checks must not pick up the user's real databases.
        env['DEFENSEVIZ_CACHE'] = str(work/'unconfigured-managed-cache')
        run([sys.executable, '-m', 'build', '--wheel', '--outdir', work/'wheels', ROOT], work)
        wheel = next((work/'wheels').glob('*.whl'))
        with zipfile.ZipFile(wheel) as archive:
            metadata_path = next(p for p in archive.namelist() if p.endswith('.dist-info/METADATA'))
            version = Parser().parsestr(archive.read(metadata_path).decode('utf-8'))['Version']
        venv.EnvBuilder(with_pip=True).create(work/'venv')
        bindir = work/'venv'/('Scripts' if os.name=='nt' else 'bin')
        python = bindir/('python.exe' if os.name=='nt' else 'python')
        run([python, '-m', 'pip', 'install', '--disable-pip-version-check', str(wheel) + '[suite]'], work)
        assert f'DefenseViz {version}' == run([python, '-m', 'defenseviz', 'version'], work).strip()
        managed_help = run([python, '-m', 'defenseviz', 'setup', '--help'], work)
        assert '--refresh-databases' in managed_help and '--source' in managed_help
        run([python, '-c',
             'import sys; from defenseviz.managed_runtime import source_state; '
             'from defenseviz.resource_download import fetch_provider_sources; '
             'assert source_state(sys.argv[1])["files"]', ROOT], work)
        setup_help = run([python, '-m', 'defenseviz', 'prepare-resources', '--help'], work)
        assert '--rebase-source' in setup_help and '--hmm-metadata' in setup_help
        run([python, '-c', 'from defenseviz.resource_setup import prepare_resources; from defenseviz.hmm_preparation import prepare'], work)
        catalog = json.loads(run([python, '-m', 'defenseviz', 'catalog', 'Wadjet', '--out', work/'catalog'], work))
        assert catalog['row_counts']['wiki_systems'] > 0
        assert catalog['row_counts']['model_definitions'] > 0
        assert (work/'catalog/model_components.tsv').is_file()
        assert (work/'catalog/catalog.json').is_file()
        run([python, '-c',
             "from importlib.resources import files; "
             "p=files('defenseviz').joinpath('data','notices'); "
             "assert all(p.joinpath(n).is_file() for n in "
             "('REFERENCE_CATALOGUE.txt','defensefinder-wiki-LICENSE.txt',"
             "'defense-finder-models-LICENSE.txt','CasFinder-LICENSE.txt'))"], work)
        run([python, '-c',
             "from importlib.resources import files; "
             "from defenseviz.native_pipeline import _verify_bundled_sources; "
             "native=files('defenseviz').joinpath('data','native'); "
             "assert native.joinpath('native_metadata.R').is_file(); "
             "assert native.joinpath('full_rm.R').is_file(); "
             "assert _verify_bundled_sources(native); "
             "assert files('defenseviz').joinpath('data','notices','DNMB_native_core.txt').is_file()"], work)
        run([python, '-m', 'defenseviz', 'demo', '--out', work/'demo'], work)
        run([python, '-m', 'defenseviz', 'validate', work/'demo/inputs/project.json'], work)
        run([python, '-m', 'defenseviz', 'run', work/'demo/inputs/project.json', '--out', work/'rerun'], work)
        a=json.loads((work/'demo/results/analysis.json').read_text(encoding='utf-8'))
        b=json.loads((work/'rerun/analysis.json').read_text(encoding='utf-8'))
        assert a==b, 'Demo and standalone run must be analytically identical'
        assert (work/'rerun/system_information.json').is_file()
        assert (work/'rerun/tables/system_information.tsv').is_file()
        assert (work/'rerun/figure2_ab.pdf').stat().st_size>1000
        assert list((work/'rerun').glob('*.svg'))
        receipt=json.loads((work/'rerun/receipt.json').read_text(encoding='utf-8'))
        assert receipt['detector_runs']==0 and receipt['exit_code']==0
        assert receipt['software']['git']['status']=='not_a_source_checkout'
        run([python, '-m', 'defenseviz', 'demo', '--suite', '--out', work/'suite-demo'], work)
        native=work/'suite-demo/inputs'
        run([python, '-m', 'defenseviz', 'inspect-suite', native/'runs', '--manifest', native/'cohort.json',
             '--supplement', native/'review.json'], work)
        run([python, '-m', 'defenseviz', 'from-suite', native/'runs', '--manifest', native/'cohort.json',
             '--supplement', native/'review.json', '--layout', native/'layout.json', '--evaluation', 'stored',
             '--title', 'Synthetic defense annotation comparison', '--out', work/'suite-rerun'], work)
        first=json.loads((work/'suite-demo/results/analysis.json').read_text(encoding='utf-8'))
        second=json.loads((work/'suite-rerun/analysis.json').read_text(encoding='utf-8'))
        assert first==second, 'Native demo and standalone suite run must have identical analytical results'
        assert (work/'suite-rerun/system_information.json').is_file()
        assert (work/'suite-rerun/tables/system_information.tsv').is_file()
        assert (work/'suite-rerun/figures/publication.pdf').stat().st_size>1000
        assert list((work/'suite-rerun/figures').glob('*.svg'))
        native_receipt=json.loads((work/'suite-rerun/receipt.json').read_text(encoding='utf-8'))
        from openpyxl import load_workbook
        workbook_path = work/'suite-rerun/tables/defense_annotations.xlsx'
        assert workbook_path.is_file()
        workbook = load_workbook(workbook_path, read_only=True, data_only=False)
        assert 'Classified CDSs' in workbook.sheetnames and 'R-M retained' in workbook.sheetnames
        assert sum(1 for _ in workbook['Classified CDSs'].iter_rows()) - 1 == second['global_summary']['integrated_CDS']
        assert sum(1 for _ in workbook['R-M retained'].iter_rows()) - 1 == second['global_summary']['retained_CDS']
        assert sum(1 for _ in workbook['Priority CDSs'].iter_rows()) - 1 == second['global_summary']['priority_CDS']
        workbook.close()
        for name in ('defense_annotations.xlsx', 'classified_CDSs.tsv', 'retained_RM_annotations.tsv', 'prioritized_annotations.tsv'):
            assert 'tables/' + name in native_receipt['outputs']
        assert native_receipt['software']['git']['status']=='not_a_source_checkout'
        assert native_receipt['detector_runs']==0 and native_receipt['exit_code']==0
        run([python, ROOT/'examples/suite_rm_evidence/make_example.py', work/'legacy-inputs'], work)
        legacy = work/'legacy-inputs'
        run([python, '-m', 'defenseviz', 'from-suite', legacy/'runs', '--manifest', legacy/'cohort.json',
             '--rm-evidence', legacy/'rm-evidence.json', '--out', work/'reevaluated'], work)
        evaluated = json.loads((work/'reevaluated/rm_evaluation.json').read_text(encoding='utf-8'))
        assert evaluated['scoring_runs'] == 1 and evaluated['policy_evaluations'] == 5
        assert sum(r['state'] == 'retained' for r in evaluated['records']) == 5
        eval_receipt = json.loads((work/'reevaluated/receipt.json').read_text(encoding='utf-8'))
        assert eval_receipt['detector_runs'] == 0 and eval_receipt['scoring_runs'] == 1
        assert eval_receipt['software']['git']['status'] == 'not_a_source_checkout'
        novel = work/'novel-system-inputs'
        run([python, ROOT/'examples/system_review/make_example.py', novel, '--verify'], work)
        run([python, '-m', 'defenseviz', 'from-suite', novel/'runs', '--manifest', novel/'cohort.json',
             '--system-review-rules', novel/'review-rules.json', '--out', work/'novel-reviewed'], work)
        review = json.loads((work/'novel-reviewed/system_review.json').read_text(encoding='utf-8'))
        assert review['accounting']['canonical_detection_count'] == 3
        assert review['accounting']['unaccounted_detection_count'] == 0
        states = {r['sample_id']: r['axes']['component_role']['status'] for r in review['gene_reviews']}
        assert states == {'example-supported': 'supported', 'example-unknown': 'unassessed',
                          'example-unversioned': 'unassessed'}
        assert (work/'novel-reviewed/tables/system_review_summary.tsv').is_file()
        print('CLEAN_INSTALL_OK')

if __name__=='__main__':main()
