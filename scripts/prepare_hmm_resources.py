#!/usr/bin/env python3
"""Compatibility launcher for the packaged exact HMM preparation implementation."""
from pathlib import Path
import sys

# A source checkout remains usable without installing its Python package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from defenseviz.hmm_preparation import *  # noqa: F403,E402

if __name__ == '__main__':
    raise SystemExit(main())
