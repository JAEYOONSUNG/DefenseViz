"""Command-line interface; all processing is local and writes fresh outputs."""
from __future__ import annotations
import argparse
import json
import sys
import tempfile
from pathlib import Path
from . import __version__


def parser():
    p = argparse.ArgumentParser(prog='defenseviz', description='Import saved DNMB results and draw editable defense annotation figures.')
    p.add_argument('--version', action='version', version=f'DefenseViz {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    v = sub.add_parser('validate', help='Validate inputs and scientific evidence contracts without output files')
    v.add_argument('config', type=Path)
    r = sub.add_parser('run', help='Evaluate saved inputs and export PDF/SVG plus auditable tables')
    r.add_argument('config', type=Path); r.add_argument('--out', type=Path, required=True)
    d = sub.add_parser('demo', help='Create a synthetic portable project and execute it')
    d.add_argument('--out', type=Path, required=True)
    d.add_argument('--suite', action='store_true', help='Demonstrate native DNMBsuite input and publication figures')
    a = sub.add_parser('assess', help='Assess supplied source-bound domain/reference evidence; no new search')
    a.add_argument('--registry', type=Path, required=True); a.add_argument('--candidates', type=Path, required=True)
    a.add_argument('--source-root', type=Path); a.add_argument('--out', type=Path, required=True)
    k = sub.add_parser('catalog', help='Browse or export bundled official defense-system reference information')
    k.add_argument('query', nargs='?', help='Optional text search for browsing; not a detection identity rule')
    k.add_argument('--out', type=Path, help='Export all matching reference records as JSON and TSV to a fresh directory')
    resources = sub.add_parser('check-resources', help='Verify DNMB databases and runtime before a full analysis')
    resources.add_argument('--engine-config', type=Path, required=True)
    resources.add_argument('--out', type=Path, help='Save the resource verification report in a new directory')
    setup = sub.add_parser('prepare-resources', help='Prepare and verify DNMB databases in a new directory')
    from .resource_setup import add_arguments
    add_arguments(setup)
    for name, help_text in [('from-suite', 'Import DNMBsuite result folders and export publication figures'),
                            ('inspect-suite', 'Validate and summarize saved DNMBsuite outputs without writing results')]:
        s = sub.add_parser(name, help=help_text)
        s.add_argument('root', type=Path, help='Completed suite folder, parent folder, or portable defenseviz_project.json')
        s.add_argument('--manifest', type=Path, help='Explicit sample paths, identities and display metadata')
        s.add_argument('--supplement', type=Path, help='Optional source-bound saved reviews')
        s.add_argument('--reference-cache', type=Path, help='Optional current DNMB REBASE cache for fresh, source-pinned homology evaluation of saved candidates')
        s.add_argument('--evaluation', choices=['auto', 'stored', 'full'], default=None,
                       help='Use configured complete DNMB execution or saved-evidence evaluation (auto); require the complete engine (full); replay only (stored)')
        s.add_argument('--engine-config', type=Path, help='DNMB runtime, full reference and native R-M HMM resource configuration')
        s.add_argument('--rm-evidence', type=Path, help='Source-bound complete R-M evidence for fresh DNMB policy evaluation')
        s.add_argument('--system-review-rules', type=Path, help='Optional exact, versioned rules for reviewing newly observed systems')
        s.add_argument('--system-review-evidence', type=Path, help='Registry and candidate paths for source-bound domain/site review')
        s.add_argument('--rscript', help='Rscript executable; explicit CLI value overrides primary and external project settings')
        s.add_argument('--strict-inputs', action=argparse.BooleanOptionalAction, default=None, help='Reject conflicting saved artifact/run metadata; explicit CLI choice overrides primary and external project settings')
        if name == 'from-suite':
            s.add_argument('--out', type=Path, required=True)
            s.add_argument('--title', default=None)
            s.add_argument('--layout', type=Path, help='Publication layout options JSON')
    sub.add_parser('version', help='Print the installed package version')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'version':
            print(f'DefenseViz {__version__}'); return 0
        if args.command == 'check-resources':
            from .resources import verify_resources
            print(json.dumps(verify_resources(args.engine_config, args.out), ensure_ascii=True, allow_nan=False))
            return 0
        if args.command == 'prepare-resources':
            from .resource_setup import from_args
            print(json.dumps(from_args(args, argv=['defenseviz', *(argv or sys.argv[1:])]),
                             ensure_ascii=True, allow_nan=False))
            return 0
        if args.command == 'catalog':
            if args.out is not None:
                from .knowledge_output import export_catalog
                response = export_catalog(args.out, query=args.query,
                                          argv=['defenseviz', *(argv or sys.argv[1:])])
            else:
                from .knowledge import search_catalog
                response = search_catalog(args.query)
            print(json.dumps(response, ensure_ascii=True, allow_nan=False))
            return 0
        if args.command in {'from-suite', 'inspect-suite'}:
            from .suite_pipeline import prepare_suite, run_suite, suite_warnings
            from .suite_project import resolve_suite_project, ANALYSIS_OPTIONS, bind_project_inputs
            kwargs = dict(manifest=args.manifest, supplement=args.supplement,
                          rm_evidence=args.rm_evidence, evaluation=args.evaluation, reference_cache=args.reference_cache, engine_config=args.engine_config,
                          system_review_rules=args.system_review_rules,
                          system_review_evidence=args.system_review_evidence,
                          rscript=args.rscript, strict_inputs=args.strict_inputs)
            overrides = {**kwargs, 'title': getattr(args, 'title', None), 'layout': getattr(args, 'layout', None)}
            source_root, settings, project_refs = resolve_suite_project(args.root, overrides)
            kwargs = {k: v for k, v in settings.items() if k in ANALYSIS_OPTIONS}
            external_overrides = {k: v for k, v in {'rscript': args.rscript, 'strict_inputs': args.strict_inputs}.items() if v is not None}
            if args.command == 'from-suite':
                response = run_suite(source_root, args.out, **kwargs,
                                     title=settings.get('title', 'DNMB integrated defense comparison'),
                                     layout=settings.get('layout'), project_refs=project_refs,
                                     external_type_projects=settings.get('external_type_projects'),
                                     external_execution_overrides=external_overrides,
                                     argv=['defenseviz', *(argv or sys.argv[1:])])
            else:
                from .analysis.saved import analyze_saved
                from .suite_pipeline import analyze_suite
                with tempfile.TemporaryDirectory(prefix='defenseviz-inspect-') as work:
                    data = prepare_suite(source_root, **kwargs, evidence_workspace=Path(work))
                    bind_project_inputs(data, project_refs)
                    data, result = analyze_suite(data)
                    external_execution = []
                    if settings.get('external_type_projects'):
                        from .suite_external import prepare_external_type_cases
                        _, external_execution, _, _ = prepare_external_type_cases(
                            settings['external_type_projects'], data,
                            Path(work) / 'external_type_cases', execution_overrides=external_overrides)
                from .knowledge import describe_systems
                information = describe_systems(data, result)
                response = {'valid': True, 'samples': len(data['samples']),
                            'external_type_case_execution': external_execution,
                            'summary': result.get('global_summary', {}),
                            'availability': result.get('availability', {}),
                            'system_review': result.get('system_review', {}).get('accounting', {}),
                            'system_information': information['accounting'],
                            'rm_evaluation': ({k: data['rm_evaluation'][k] for k in
                                ('operation', 'policy_evaluations', 'policy_sha256', 'scoring_runs', 'thresholds_changed')}
                                if data.get('rm_evaluation') else None),
                            'warnings': suite_warnings(data, result),
                            'conflicts': data.get('suite', {}).get('conflicts', [])}
            print(json.dumps(response, ensure_ascii=True, allow_nan=False))
        elif args.command == 'validate':
            from .io import load_project
            from .pipeline import evaluate_project
            data = load_project(args.config)
            result = evaluate_project(data)
            print(json.dumps({'valid': True, 'samples': len(data['samples']), 'call_rows': len(data['calls']), 'evidence_records': len(data['evidence_records']), 'warnings': result.get('warnings', [])}, ensure_ascii=True))
        elif args.command == 'assess':
            from .assessment import assess
            print(json.dumps(assess(args.registry, args.candidates, args.out, source_root=args.source_root, argv=['defenseviz', *(argv or sys.argv[1:])]), ensure_ascii=True))
        elif args.command == 'demo':
            if args.suite:
                from .suite_demo import demo_suite
                response = demo_suite(args.out)
            else:
                from .demo import demo
                response = demo(args.out)
            print(json.dumps(response, ensure_ascii=True))
        else:
            from .pipeline import run
            print(json.dumps(run(args.config, args.out, argv=['defenseviz', *(argv or sys.argv[1:])]), ensure_ascii=True))
        return 0
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as exc:
        print(f'defenseviz: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
