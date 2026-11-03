"""Command-line interface; all processing is local and writes fresh outputs."""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from . import __version__


def parser():
    p = argparse.ArgumentParser(prog='defenseviz', description='Integrate saved defense annotations, reevaluate DNMB evidence and draw Figure 2 a/b.')
    p.add_argument('--version', action='version', version=f'DefenseViz {__version__}')
    sub = p.add_subparsers(dest='command', required=True)
    v = sub.add_parser('validate', help='Validate inputs and scientific evidence contracts without output files')
    v.add_argument('config', type=Path)
    r = sub.add_parser('run', help='Evaluate saved inputs and export PDF/SVG plus auditable tables')
    r.add_argument('config', type=Path); r.add_argument('--out', type=Path, required=True)
    d = sub.add_parser('demo', help='Create a synthetic portable project and execute it')
    d.add_argument('--out', type=Path, required=True)
    a = sub.add_parser('assess', help='Assess supplied source-bound domain/reference evidence; no new search')
    a.add_argument('--registry', type=Path, required=True); a.add_argument('--candidates', type=Path, required=True)
    a.add_argument('--source-root', type=Path); a.add_argument('--out', type=Path, required=True)
    sub.add_parser('version', help='Print the installed package version')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        if args.command == 'version':
            print(f'DefenseViz {__version__}'); return 0
        if args.command == 'validate':
            from .io import load_project
            from .pipeline import evaluate_project
            data = load_project(args.config)
            result = evaluate_project(data)
            print(json.dumps({'valid': True, 'samples': len(data['samples']), 'call_rows': len(data['calls']), 'evidence_records': len(data['evidence_records']), 'warnings': result.get('warnings', [])}, ensure_ascii=False))
        elif args.command == 'assess':
            from .assessment import assess
            print(json.dumps(assess(args.registry, args.candidates, args.out, source_root=args.source_root, argv=['defenseviz', *(argv or sys.argv[1:])]), ensure_ascii=False))
        elif args.command == 'demo':
            from .demo import demo
            print(json.dumps(demo(args.out), ensure_ascii=False))
        else:
            from .pipeline import run
            print(json.dumps(run(args.config, args.out, argv=['defenseviz', *(argv or sys.argv[1:])]), ensure_ascii=False))
        return 0
    except (ValueError, OSError, KeyError, TypeError, RuntimeError) as exc:
        print(f'defenseviz: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
