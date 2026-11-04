"""Generate invented legacy exports and explicit Evidence; never use real data."""
from dataclasses import asdict
from pathlib import Path
import csv
import hashlib
import json
import sys

from defenseviz.analysis.policy import Policy
from defenseviz.io.suite import load_suite
from defenseviz.suite_demo import write_suite_project


def make_example(output):
    output = Path(output).absolute()
    inputs = write_suite_project(output)
    hypotheses = {}
    artifacts = {}
    for folder in sorted(inputs['root'].iterdir()):
        path = folder / 'dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv'
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        for row in rows:
            hypotheses[(folder.name, row['query'])] = (row['final_family'], row['final_role'] or 'unknown')
        # This creates an old-format synthetic fixture, never edits user inputs.
        with path.open('w') as handle:
            writer = csv.DictWriter(handle, ['query', 'blast_match', 'blast_identity'], delimiter='\t', lineterminator='\n')
            writer.writeheader()
            writer.writerows({k: row[k] for k in writer.fieldnames} for row in rows)
        artifacts[folder.name] = path
    data = load_suite(inputs['root'], manifest=inputs['manifest'])
    proteins = {(g['sample_id'], g['gene_id']): g['protein_sha256'] for g in data['genes']}
    rows = []
    for (sid, gid), (family, role) in hypotheses.items():
        rows.append(dict(sample_id=sid, gene_id=gid, protein_sha256=proteins[(sid, gid)], evidence=dict(
            family=family, role=role, valid_sequence=True, identity=0.9,
            query_coverage=0.95, reference_coverage=0.95, evalue=1e-30,
            role_compatible=True, family_compatible=True, hmm_role_ambiguous=False,
            profiles=[], raw_nuclease_classes=[], legacy_nuclease=False, legacy_motor=False,
            legacy_mtase_motif=False, old_strong=True, old_moderate=False, old_profile=False,
            annotation=True, old_motif=False, old_context=True, strict_structure=False,
            alternative_annotation=False, conflict=False, partial=False, divergent_context=False,
            curated_identity=None, curated_query_coverage=None, curated_reference_coverage=None,
            curated_evalue=None, curated_role_compatible=False, curated_family_compatible=False,
            reference_nuclease_anchor_support=False, reference_nuclease_anchor_profile_indices=[])))
    source = output / 'synthetic-evidence-source.json'
    source.write_text(json.dumps({'description': 'Invented numerical and boolean assessments for demonstration only. No biological activity is asserted.', 'records': rows}, indent=2) + '\n')
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    refs = [{'path': source.name, 'sha256': sha(source)}]
    evidence = dict(schema_version=1, selected_policy=asdict(Policy()),
        input_sha256={p.relative_to(output).as_posix(): sha(p) for p in artifacts.values()},
        records=[dict(row, source_refs=refs) for row in rows],
        candidate_scope=[dict(sample_id=sid, status='complete',
            gene_ids=[g for s,g in hypotheses if s == sid], source_refs=refs) for sid in artifacts])
    path = output / 'rm-evidence.json'
    path.write_text(json.dumps(evidence, indent=2) + '\n')
    return dict(**inputs, rm_evidence=path)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python make_example.py NEW_INPUT_DIRECTORY')
    print(json.dumps({k: str(v) for k,v in make_example(sys.argv[1]).items()}))
