"""External legacy exports preserve taxonomy and coherent reference subjects."""
import csv
import hashlib

import pytest

from defenseviz.io.suite import load_suite


def table(path, headers, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(rows)


def fixture(tmp_path):
    root = tmp_path / 'legacy-isolate'
    table(root / 'dnmb_module_rebasefinder/rebasefinder_input.tsv',
          ['locus_tag', 'translation'], [['geneA', 'MAAA']])
    return root


def references(root, values):
    path = root / 'dnmb_module_rebasefinder/DNMB_REBASEfinder_augmented_hits.tsv'
    table(path, list(values), [list(values.values())])
    return path


def ref_status(data):
    return next(r['status'] for r in data['coverage']
                if r['caller'] == 'REBASE' and r['endpoint'] == 'reference_matches')


@pytest.mark.parametrize('label', ['typeI_context:R:sequence_motif', 'typeIII_context:R:annotation'])
def test_context_marker_uses_its_explicit_supplemental_subject_pair(tmp_path, label):
    root = fixture(tmp_path)
    source = references(root, dict(query='geneA', source='rebasefinder',
        evidence_mode='operon_context_rebase_blast', hit_label=label, blast_identity='0.41',
        family_id='Type I', type='Type I', enzyme_role='R',
        supplemental_blast_match='Actual.Reference', supplemental_blast_identity='0.73'))
    original = source.read_bytes()
    data = load_suite(root)
    ref, = data['suite']['reference_matches']
    assert (ref['reference_id'], ref['identity'], ref['match_origin']) == ('Actual.Reference', .73, 'supplemental_blast')
    assert ref_status(data) == 'complete'
    call, = [c for c in data['calls'] if c['caller'] == 'REBASE']
    assert call['reference_id'] == 'Actual.Reference' and call['match_identity'] == .73
    assert call['family'] == 'Unresolved' and call['rm_type'] is None and call['role'] is None
    assert call['reference_family_source'] == 'unassessed'
    assert call['native_fields']['hit_label'] == label
    assert call['native_fields']['family_id'] == 'Type I'
    assert len(data['saved_decisions']) == 1 and data['saved_decisions'][0]['state'] == 'unassessed'
    assert source.read_bytes() == original
    assert ref['source_ref']['sha256'] == hashlib.sha256(original).hexdigest()


@pytest.mark.parametrize('identity,state', [('NA', 'complete'), ('0.41', 'partial')])
def test_context_marker_without_subject_does_not_invent_reference(tmp_path, identity, state):
    root = fixture(tmp_path)
    label = 'typeI_context:R:sequence_motif'
    references(root, dict(query='geneA', source='rebasefinder', evidence_mode='operon_context_rebase_blast',
        hit_label=label, blast_identity=identity, support='supplemental_REBASE_BLAST=ProseOnly; identity=70%'))
    data = load_suite(root)
    assert not data['suite']['reference_matches']
    assert not [c for c in data['calls'] if c['caller'] == 'REBASE']
    assert data['suite']['features'][0]['native_fields']['hit_label'] == label
    assert len(data['saved_decisions']) == 1
    assert ref_status(data) == state


@pytest.mark.parametrize('subject,identity', [('Actual.Reference', 'NA'), ('NA', '0.73')])
def test_incomplete_supplemental_pair_is_partial_and_does_not_borrow(tmp_path, subject, identity):
    root = fixture(tmp_path)
    references(root, dict(query='geneA', source='rebasefinder', evidence_mode='operon_context_rebase_blast',
        hit_label='typeIII_context:R:annotation', blast_identity='0.99',
        supplemental_blast_match=subject, supplemental_blast_identity=identity))
    data = load_suite(root)
    assert not data['suite']['reference_matches'] and ref_status(data) == 'partial'
    assert len(data['saved_decisions']) == 1


def test_real_primary_reference_still_precedes_supplemental_pair(tmp_path):
    root = fixture(tmp_path)
    references(root, dict(query='geneA', source='rebasefinder', evidence_mode='operon_context_rebase_blast',
        hit_label='Primary.Reference', blast_identity='.6',
        supplemental_blast_match='Supplemental.Reference', supplemental_blast_identity='.8'))
    data = load_suite(root)
    ref, = data['suite']['reference_matches']
    assert (ref['reference_id'], ref['identity'], ref['match_origin']) == ('Primary.Reference', .6, 'native_blast')


@pytest.mark.parametrize('fqn,parent,family,category,rm_type', [
    ('defense-finder-models/DefenseFinder/Wadjet/Wadjet_II', 'Wadjet', 'Wadjet_II', 'Unresolved', None),
    ('defense-finder-models/RM/RM/RM_Type_IV_1', 'RM', 'RM_Type_IV_1', 'R-M', 'Type IV'),
    ('CasFinder/CAS_Class1-Subtype-I-A', 'CAS', 'CAS_Class1-Subtype-I-A', 'CRISPR-Cas', None),
    ('defense-finder-models/DefenseFinder/NewFamily/NewFamily_I', 'NewFamily', 'NewFamily_I', 'Unresolved', None),
])
def test_legacy_df_projects_missing_taxonomy_from_exact_fqn_without_editing_native_fields(tmp_path, fqn, parent, family, category, rm_type):
    root = fixture(tmp_path)
    source = root / 'dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv'
    fields = dict(hit_id='geneA', sys_id='systemA', model_fqn=fqn, gene_name='Unknown_component')
    table(source, list(fields), [list(fields.values())])
    original = source.read_bytes()
    data = load_suite(root)
    call, = [c for c in data['calls'] if c['caller'] == 'DefenseFinder']
    assert call['native_fields'] == fields and call['native_type'] is None
    assert call['native_family'] == family
    assert (call['category'], call['rm_type']) == (category, rm_type)
    assert call['native_taxonomy_projection']['type'] == parent
    assert call['native_taxonomy_projection']['subtype'] == family
    assert call['role'] is None and not data['effector_evidence'] and not data['saved_decisions']
    assert call['source_ref']['sha256'] == hashlib.sha256(original).hexdigest()
    assert source.read_bytes() == original


@pytest.mark.parametrize('fqn,subtype', [
    ('custom/Wadjet/Wadjet_II', ''),
    ('defense-finder-models/DefenseFinder/Wadjet/', ''),
    ('defense-finder-models/DefenseFinder/Wadjet/Wadjet_II', 'Wadjet_I'),
])
def test_missing_df_type_does_not_guess_unknown_or_conflicting_taxonomy(tmp_path, fqn, subtype):
    root = fixture(tmp_path)
    source = root / 'dnmb_module_defensefinder/defensefinder_best_solution_genes.tsv'
    table(source, ['hit_id', 'sys_id', 'model_fqn', 'subtype'], [['geneA', 'systemA', fqn, subtype]])
    with pytest.raises(ValueError, match='model_fqn'):
        load_suite(root)


def test_suite_gene_workbook_keeps_literal_hyperlink_transform_provenance(tmp_path):
    openpyxl = pytest.importorskip('openpyxl')
    root = tmp_path / 'legacy-workbook'
    root.mkdir()
    (root / 'dnmb_module_rebasefinder').mkdir()
    source = root / 'isolate_total.xlsx'
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = '1.GenBank_table'
    sheet.append(['locus_tag', 'translation', 'external_link'])
    sheet.append(['geneA', 'MAAA', '=HYPERLINK("https://example.org/entry","Entry")'])
    workbook.save(source)
    original = source.read_bytes()
    data = load_suite(root)
    assert data['suite']['annotations'][0]['native_fields']['external_link'] == 'Entry'
    transform, = data['suite']['table_transforms']
    assert transform['cell'] == 'C2' and transform['label'] == 'Entry'
    assert transform['sha256'] == hashlib.sha256(original).hexdigest()
    assert source.read_bytes() == original
