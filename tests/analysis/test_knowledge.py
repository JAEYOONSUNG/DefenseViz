"""Independent synthetic counterexamples for metadata-only catalogue linkage."""
from copy import deepcopy
import hashlib
import json

import pytest

from defenseviz.knowledge import (build_catalog, catalog_tables, describe_systems,
                                  information_tables, load_catalog, search_catalog, _json)


def ref(path):
    return dict(path=path, sha256=hashlib.sha256(path.encode()).hexdigest())


def wiki_row(name, *, slug=None, alias=None):
    slug = slug or name
    return dict(id='defensefinder-wiki:' + slug, title=name, slug=slug,
                source=ref('content/3.defense-systems/' + slug + '.md'),
                aliases=[] if alias is None else [dict(value=alias, source_field='refseqSystemName', source_line=3)],
                mechanism_summary='Synthetic reference description; not evidence for the query.',
                facts=dict(effector=dict(state='unknown', value=None)))


def definition(name, *, namespace='DefenseFinder', family=None, source='defense-finder-models'):
    path = '/'.join(v for v in (namespace, family, name) if v) if source != 'CasFinder' else name
    return dict(definition_id=source + '/' + path, source_id=source, source_definition_key=path,
                package_version='synthetic-1', source_namespace=namespace, source_family=family,
                model_name=name, defense_scope='anti_defense' if namespace == 'ADF' else 'host_defense',
                source_ref=ref('definitions/' + path + '.xml'), exact_fqn_aliases=[],
                components=[dict(xpath='/model[1]/gene[1]', native_component_id='not_a_nuclease_assertion',
                                 declared_presence='mandatory', relationship='top_level')])


@pytest.fixture
def payloads():
    wiki = dict(schema_version=1, source=dict(repository='https://example.invalid/wiki', commit='a' * 40),
                systems=[wiki_row('Future'), wiki_row('WikiOnly'), wiki_row('Cas', alias='CAS'),
                         wiki_row('CaseSensitive'), wiki_row('CoCoNut_I-A')])
    sources = [dict(source_id=name, repository='https://example.invalid/' + name, commit='b' * 40,
                    package_version='synthetic-1') for name in ('defense-finder-models', 'CasFinder')]
    models = dict(schema_version=1, sources=sources, definitions=[
        definition('Future_I', family='Future'), definition('Future_II', family='Future'),
        definition('ModelOnly'), definition('SameName', family='First'),
        definition('SameName', family='Second'), definition('AntiThing', namespace='ADF', family='Anti_Test'),
        definition('CaseSensitive'), definition('CoCoNut_I-A'), definition('CAS_Test', namespace='CasFinder', source='CasFinder')])
    return wiki, models


@pytest.fixture
def catalog(payloads):
    return build_catalog(*payloads)


def call(name='Future_I', *, fqn=True, provider='DefenseFinder', **extra):
    native = dict(subtype=name)
    if fqn:
        native['model_fqn'] = ('defense-finder-models/DefenseFinder/Future/' + name) if fqn is True else fqn
    return dict(sample_id='s1', gene_id='g1', caller=provider, family='display label', native_family=name,
                native_fields=native, native_model='component_HMM', system_id='system-1', run_id='run-1',
                replicon_id='chr', evidence_kind='final_member', selected=True, **extra)


def info(catalog, row, result=None):
    return describe_systems(dict(calls=[row]), {} if result is None else result, catalog)


def resolution(catalog, row):
    return info(catalog, row)['detections'][0]


def review_for(calls):
    # Construct exact canonical IDs independently of the knowledge module.
    def identity(value):
        raw = json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()
        return 'detection:' + hashlib.sha256(raw).hexdigest()
    units = []
    for number, row in enumerate(calls):
        key = [row['sample_id'], row['gene_id']]
        did = identity([key, row.get('caller'), row.get('run_id'), row.get('replicon_id'), row.get('system_id'),
                        row.get('evidence_kind', 'final_member'), row.get('native_family', row.get('family')),
                        row.get('family'), row.get('native_model'), row.get('reference_id'), row.get('selected', True)])
        units.append(dict(unit_id='unit' + str(number), unit_kind='native_system', sample_id=key[0],
                          provider=row['caller'], run_id=row.get('run_id'), replicon_id=row.get('replicon_id'),
                          native_system_id=row.get('system_id'), member_keys=[key], selected_member_keys=[key],
                          native_families=[row['native_family']], native_models=[row['native_model']], detection_ids=[did],
                          axes=dict(activity=dict(status='unassessed'))))
    return dict(system_review=dict(units=units))


def test_all_source_records_are_distinct_and_model_only_preserved(catalog):
    assert catalog['coverage']['wiki_pages'] == 5
    assert catalog['coverage']['model_definitions'] == 9
    assert any(r['model_name'] == 'ModelOnly' for r in catalog['model_definitions'])
    assert not any(r['wiki_id'].endswith('WikiOnly') for r in catalog['links'])
    future = [r for r in catalog['links'] if r['wiki_id'].endswith(':Future')]
    assert len({r['definition_id'] for r in future}) == 2
    assert {r['relationship'] for r in future} == {'parent_family_context'}


@pytest.mark.parametrize('version,expected', [(None, 'version_unavailable'), ('synthetic-1', 'same_declared_package_version'), ('older', 'different_version')])
def test_exact_fqn_keeps_version_relation_separate(catalog, version, expected):
    r = resolution(catalog, call(model_version=version))
    assert r['identity_status'] == 'exact_fqn'
    assert r['version_relation'] == expected
    assert r['wiki_ids'] == ['defensefinder-wiki:Future']
    assert r['query_function_inferred'] is False


@pytest.mark.parametrize('fields', [dict(software_version='synthetic-1'),
                                  dict(model_version='synthetic-1', model_version_source='native_row.database_version'),
                                  dict(database_version='synthetic-1'), dict(model_sha256='f' * 64)])
def test_incompatible_version_or_hash_field_never_binds_definition(catalog, fields):
    assert resolution(catalog, call(**fields))['version_relation'] == 'version_unavailable'


def test_contradictory_definition_hash_and_version_is_explicit(catalog):
    row = call(model_version='old')
    model = next(r for r in catalog['model_definitions'] if r['model_name'] == 'Future_I')
    row['native_fields']['model_definition_sha256'] = model['source_ref']['sha256']
    assert resolution(catalog, row)['version_relation'] == 'provenance_conflict'


def test_exact_hash_and_commit_have_distinct_strongest_binding(catalog):
    row = call()
    row['native_fields']['model_package_commit'] = 'b' * 40
    assert resolution(catalog, row)['version_relation'] == 'same_source_commit'
    row['native_fields']['model_definition_sha256'] = catalog['model_definitions'][0]['source_ref']['sha256']
    assert resolution(catalog, row)['version_relation'] == 'same_definition_hash'


def test_conflicting_reported_model_versions_are_not_silently_prioritized(catalog):
    row = call(model_version='old-1', source_version=dict(models_version='old-2'))
    r = resolution(catalog, row)
    assert r['version_relation'] == 'provenance_conflict'
    assert {v['value'] for v in r['native_model_versions_reported']} == {'old-1', 'old-2'}


def test_explicit_anti_namespace_cannot_get_host_model_from_name(catalog):
    row = call(fqn=False); row['native_fields']['namespace'] = 'ADF'
    r = resolution(catalog, row)
    assert r['identity_status'] == 'conflicting_identity'
    assert r['defense_scope'] == 'anti_defense' and not r['definition_ids']


def test_explicit_unresolved_fqn_does_not_fallback_to_native_family(catalog):
    r = resolution(catalog, call(fqn='defense-finder-models/DefenseFinder/unknown/Future_I'))
    assert r['identity_status'] == 'unmatched'
    assert r['definition_ids'] == r['wiki_ids'] == []


def test_native_subtype_contradiction_is_not_forgiven(catalog):
    row = call(); row['native_fields']['subtype'] = 'Future_II'
    r = resolution(catalog, row)
    assert r['identity_status'] == 'conflicting_identity'
    assert not r['definition_ids'] and len(r['candidate_definition_ids']) == 1


def test_native_family_contradiction_is_not_forgiven(catalog):
    row = call(); row['native_family'] = 'Future_II'
    assert resolution(catalog, row)['identity_status'] == 'conflicting_identity'


def test_exact_unique_label_has_context_not_fqn_identity(catalog):
    r = resolution(catalog, call(fqn=False))
    assert r['identity_status'] == 'exact_label_context'
    assert len(r['definition_ids']) == 1


def test_same_leaf_in_multiple_namespaces_is_ambiguous(catalog):
    r = resolution(catalog, call('SameName', fqn=False))
    assert r['identity_status'] == 'ambiguous'
    assert len(r['candidate_definition_ids']) == 2 and not r['definition_ids']


@pytest.mark.parametrize('name', ['casesensitive', 'CoCoNut_I_A', 'Future', 'component_HMM'])
def test_no_case_punctuation_parent_or_component_model_guess(catalog, name):
    r = resolution(catalog, call(name, fqn=False))
    assert r['definition_ids'] == []
    assert r['identity_status'] == ('wiki_context_only' if name == 'Future' else 'unmatched')


def test_display_alias_does_not_supply_native_identity(catalog):
    row = call('Unknown', fqn=False); row['family'] = 'Future_I'
    assert not resolution(catalog, row)['definition_ids']


@pytest.mark.parametrize('provider', ['PADLOC', 'DefensePredictor', 'REBASE'])
def test_other_providers_do_not_borrow_defensefinder_model(catalog, provider):
    r = resolution(catalog, call(provider=provider))
    assert r['identity_status'] == 'provider_out_of_scope' and not r['definition_ids']
    r = resolution(catalog, call('Future', provider=provider, fqn=False))
    assert r['identity_status'] == 'wiki_context_only' and r['wiki_ids'] and not r['definition_ids']


def test_anti_namespace_scope_is_context_and_unknown_anti_name_is_unresolved(catalog):
    r = resolution(catalog, call('AntiThing', fqn='defense-finder-models/ADF/Anti_Test/AntiThing'))
    assert r['defense_scope'] == 'anti_defense' and r['identity_status'] == 'exact_fqn'
    assert resolution(catalog, call('AntiUnknown', fqn=False))['defense_scope'] == 'unresolved'


def test_casfinder_is_separate_package(catalog):
    r = resolution(catalog, call('CAS_Test', fqn='CasFinder/CAS_Test'))
    assert r['definition_ids'] == ['CasFinder/CAS_Test']
    r = resolution(catalog, call('Future_I', fqn='different-package/DefenseFinder/Future/Future_I'))
    assert r['identity_status'] == 'provider_out_of_scope'


def test_wiki_only_and_model_only_are_preserved(catalog):
    r = resolution(catalog, call('WikiOnly', fqn=False)); assert r['wiki_ids'] and not r['definition_ids']
    r = resolution(catalog, call('ModelOnly', fqn='defense-finder-models/DefenseFinder/ModelOnly'))
    assert r['definition_ids'] and not r['wiki_ids']


def test_source_declared_alias_is_wiki_context_only(catalog):
    r = resolution(catalog, call('CAS', fqn=False))
    assert r['wiki_ids'] == ['defensefinder-wiki:Cas'] and not r['definition_ids']


def test_documented_fqn_alias_needs_exact_proof(payloads):
    wiki, models = deepcopy(payloads)
    models['definitions'][0]['exact_fqn_aliases'] = [dict(value='old-package/exact-Future', source_ref=ref('release-notes.md'))]
    r = resolution(build_catalog(wiki, models), call(fqn='old-package/exact-Future'))
    assert r['identity_status'] == 'documented_alias'


def test_same_gene_two_runs_keeps_independent_unit_context(catalog):
    first, second = call(), call('Future_II')
    second['run_id'] = 'run-2'
    calls = [first, second]; result = review_for(calls)
    out = describe_systems(dict(calls=calls), result, catalog)
    assert len(out['units']) == len(out['detections']) == 2
    assert out['units'][0]['definition_ids'] != out['units'][1]['definition_ids']
    assert all(u['coding_sequences'] == 1 for u in out['units'])


def test_many_components_many_links_do_not_multiply_cds(catalog):
    first, second = call(), call(); second['gene_id'] = 'g2'
    calls = [first, second]; result = review_for(calls)
    a, b = result['system_review']['units']
    a['member_keys'] += b['member_keys']; a['selected_member_keys'] += b['selected_member_keys']; a['detection_ids'] += b['detection_ids']
    result['system_review']['units'] = [a]
    out = describe_systems(dict(calls=calls), result, catalog)
    assert len(out['units']) == 1 and out['units'][0]['coding_sequences'] == 2
    assert out['accounting']['detection_count'] == 2


def test_review_only_unit_has_explicit_no_native_identity(catalog):
    result = review_for([call()]); u = deepcopy(result['system_review']['units'][0])
    u.update(unit_id='review-only', unit_kind='review_only_gene', detection_ids=[], member_keys=[['s1', 'orphan']])
    result['system_review']['units'].append(u)
    out = describe_systems(dict(calls=[call()]), result, catalog)
    assert out['units'][1]['identity_statuses'] == ['no_native_model_identity']
    assert out['units'][1]['coding_sequences'] == 1


def test_normalized_calls_without_review_do_not_invent_native_systems(catalog):
    out = info(catalog, call())
    assert out['units'][0]['unit_kind'] == 'ungrouped_source_annotation'
    assert out['units'][0]['native_system_id'] is None
    assert out['units'][0]['reported_native_system_id'] == 'system-1'


@pytest.mark.parametrize('failure', ['missing', 'duplicate', 'wrong_gene'])
def test_review_detection_accounting_fails_closed(catalog, failure):
    calls = [call()]; result = review_for(calls)
    if failure == 'missing':calls = []
    elif failure == 'duplicate':calls.append(deepcopy(calls[0]))
    else:calls[0]['gene_id'] = 'other'
    with pytest.raises(ValueError, match='reconcile|unaccounted'):
        describe_systems(dict(calls=calls), result, catalog)


def test_metadata_cannot_promote_scientific_fields_or_mutate_inputs(catalog):
    data = dict(calls=[call()], genes=[dict(sample_id='s1', gene_id='g1')])
    result = review_for(data['calls']); result.update(gene_ledger=[dict(effector=False, associated=False, priority=False)], cells=[dict(N=1, effector_CDS=0)])
    before = deepcopy((data, result, catalog))
    out = describe_systems(data, result, catalog)
    assert (data, result, catalog) == before
    assert out['scientific_membership_changed'] is False
    assert out['query_function_inferred'] is False
    assert 'axes' not in out['units'][0]


def test_search_is_forgiving_but_does_not_change_exact_runtime_binding(catalog):
    search = search_catalog('casesensitive', catalog)
    assert search['search']['matched_model_definitions'] == 1
    assert resolution(catalog, call('casesensitive', fqn=False))['identity_status'] == 'unmatched'


def test_model_search_retains_parent_context_without_all_sibling_models(catalog):
    model_id = 'defense-finder-models/DefenseFinder/Future/Future_I'
    out = search_catalog(model_id, catalog)
    assert out['search']['mode'] == 'exact_record_id'
    assert [r['model_name'] for r in out['model_definitions']] == ['Future_I']
    assert [r['title'] for r in out['wiki_systems']] == ['Future']
    assert out['search']['direct_wiki_ids'] == []
    assert out['search']['linked_context_wiki_ids'] == ['defensefinder-wiki:Future']
    assert {r['definition_id'] for r in out['links']} == {model_id}


def test_explicit_related_system_names_are_context_never_aliases(payloads):
    wiki, models = deepcopy(payloads)
    wiki['systems'][0]['related_names'] = [dict(value='FutureRelated', relationship='group_member', source_lines=[8],
                                              scope='descriptive_name_relationship_not_runtime_alias')]
    cat = build_catalog(wiki, models)
    r = resolution(cat, call('FutureRelated', fqn=False))
    assert r['identity_status'] == 'unmatched' and not r['definition_ids'] and not r['wiki_ids']
    assert r['related_wiki_context'][0]['wiki_id'] == 'defensefinder-wiki:Future'
    assert r['related_wiki_context'][0]['relationship'] == 'group_member'
    assert resolution(cat, call('futurerelated', fqn=False))['related_wiki_context'] == []
    out = info(cat, call('FutureRelated', fqn=False))
    assert out['accounting']['related_wiki_context_detections'] == 1
    assert out['accounting']['wiki_linked_detections'] == 0
    assert out['units'][0]['related_reference_mechanisms']


def test_related_component_name_is_not_system_context(payloads):
    wiki, models = deepcopy(payloads)
    wiki['systems'][0]['related_names'] = [dict(value='sensor', relationship='component_name', source_lines=[8],
                                              scope='descriptive_name_relationship_not_runtime_alias')]
    assert resolution(build_catalog(wiki, models), call('sensor', fqn=False))['related_wiki_context'] == []


def test_component_metadata_disagreement_is_visible_without_changing_identity(catalog):
    r = resolution(catalog, call(model_version='synthetic-1'))
    assert r['identity_status'] == 'exact_fqn'
    assert r['current_definition_component_context']['state'] == 'not_listed_in_current_definition'
    assert 'native_component_not_listed_in_current_definition' in r['reason_codes']
    row = call(); row['native_fields']['hit_gene_ref'] = 'not_a_nuclease_assertion'
    r = resolution(catalog, row)
    assert r['current_definition_component_context']['state'] == 'partial_declared_component_names'
    assert r['current_definition_component_context']['not_declared_names'] == ['component_HMM']
    assert r['query_function_inferred'] is False


def test_provider_wide_version_cannot_bind_a_different_model_package(catalog):
    row = call('CAS_Test', fqn='CasFinder/CAS_Test', model_version='synthetic-1',
               model_version_source='saved_module_metadata.models_version',
               source_version=dict(models_version='synthetic-1'))
    r = resolution(catalog, row)
    assert r['version_relation'] == 'version_unavailable'
    assert 'version_scope_unbound' in r['reason_codes']
    row['source_version']['models_package'] = 'CasFinder'
    assert resolution(catalog, row)['version_relation'] == 'same_declared_package_version'


def test_duplicate_raw_identity_conflict_is_not_hidden(catalog):
    row = call(); second = deepcopy(row['native_fields'])
    second['model_fqn'] = 'defense-finder-models/ADF/Future/Future_I'
    row['native_observations'] = [row['native_fields'], second]
    r = resolution(catalog, row)
    assert r['identity_status'] == 'conflicting_identity' and not r['definition_ids']


def test_package_resource_load_uses_exact_bytes_and_no_repository(payloads, tmp_path, monkeypatch):
    from defenseviz import knowledge
    data = tmp_path / 'data'; data.mkdir()
    for name, payload in zip(knowledge.RESOURCE_NAMES, payloads):
        (data / name).write_bytes(json.dumps(payload).encode())
    monkeypatch.setattr(knowledge, 'files', lambda name: tmp_path)
    out = load_catalog()
    assert out['coverage']['model_definitions'] == 9
    for resource in out['provenance']['resources']:
        assert resource['sha256'] == hashlib.sha256((data / resource['resource'].split('/')[-1]).read_bytes()).hexdigest()


def test_tables_retain_components_and_exact_source_records(catalog):
    tables = catalog_tables(catalog)
    assert len(tables['wiki_systems']) == 5 and len(tables['model_definitions']) == 9
    assert len(tables['model_components']) == 9
    assert tables['wiki_systems'][0]['effector_state'] == 'unknown'
    assert tables['wiki_systems'][0]['source_commit'] == 'a' * 40
    assert 'source_mechanism_statement_state' in tables['wiki_systems'][0]
    out = info(catalog, call()); flat = information_tables(out)
    assert flat['units'][0]['reference_titles'] == 'Future'
    assert isinstance(flat['units'][0]['reference_mechanism_summaries'], str)
    assert [{k:v for k,v in row.items() if k != 'native_model_version_first_reported'} for row in flat['detections']] == out['detections']
    flat['units'].clear(); assert out['units']


@pytest.mark.parametrize('mutation', ['duplicate_model','bad_id','unknown_source','scope_conflict','version_conflict','bad_hash','unproven_alias','duplicate_xpath','duplicate_wiki','unproven_wiki_alias'])
def test_invalid_catalogue_is_rejected(payloads, mutation):
    wiki, models = deepcopy(payloads); row = models['definitions'][0]
    if mutation == 'duplicate_model':models['definitions'].append(deepcopy(row))
    elif mutation == 'bad_id':row['definition_id'] = 'collapsed/Future_I'
    elif mutation == 'unknown_source':row['source_id'] = 'unknown'
    elif mutation == 'scope_conflict':models['definitions'][5]['defense_scope'] = 'host_defense'
    elif mutation == 'version_conflict':row['package_version'] = 'wrong'
    elif mutation == 'bad_hash':row['source_ref']['sha256'] = 'invalid'
    elif mutation == 'unproven_alias':row['exact_fqn_aliases'] = [dict(value='alias')]
    elif mutation == 'duplicate_xpath':row['components'].append(deepcopy(row['components'][0]))
    elif mutation == 'duplicate_wiki':wiki['systems'].append(deepcopy(wiki['systems'][0]))
    else:wiki['systems'][0]['aliases'] = [dict(value='guessed', source_field='fuzzy', source_line=1)]
    with pytest.raises(ValueError):build_catalog(wiki, models)


@pytest.mark.parametrize('raw', ['{"a":1,"a":2}', '{"a":NaN}', '{"a":Infinity}'])
def test_ambiguous_or_nonfinite_json_is_rejected(raw):
    with pytest.raises(ValueError):_json(raw)


def test_forged_combined_derived_links_are_not_trusted(catalog):
    catalog['links'].append(dict(definition_id='missing', wiki_id='missing'))
    with pytest.raises(ValueError, match='links or coverage'):
        info(catalog, call())


@pytest.mark.parametrize('target', ['wiki', 'model', 'namespace', 'source_path', 'model_name'])
def test_source_inventory_or_identity_contradiction_is_rejected(payloads, target):
    wiki, models = deepcopy(payloads)
    if target == 'wiki':wiki['source']['coverage'] = dict(system_pages=999)
    elif target == 'model':models['coverage'] = dict(definition_count=999)
    elif target == 'namespace':models['definitions'][0]['source_namespace'] = 'ADF'
    elif target == 'source_path':models['definitions'][0]['source_ref']['path'] = 'somewhere-else.xml'
    else:models['definitions'][0]['model_name'] = 'different'
    with pytest.raises(ValueError):build_catalog(wiki, models)


def test_boolean_schema_version_is_not_accepted(catalog):
    catalog['schema_version'] = True
    with pytest.raises(ValueError):info(catalog, call())


def test_bundled_official_sources_and_exact_parent_context_work_outside_repository(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cat = load_catalog()
    assert cat['coverage']['wiki_pages'] == cat['wiki_source']['coverage']['system_pages'] > 0
    assert cat['coverage']['model_definitions'] == cat['source_inventory_coverage']['definition_count'] > 0
    assert set(cat['coverage']['model_definitions_by_source']) == {'defense-finder-models', 'CasFinder'}
    exact = search_catalog('defense-finder-models/DefenseFinder/Wadjet/Wadjet_II', cat)
    assert exact['search']['mode'] == 'exact_record_id'
    assert [r['model_name'] for r in exact['model_definitions']] == ['Wadjet_II']
    avs = search_catalog('defense-finder-models/DefenseFinder/Avs/Avs10', cat)
    assert any(r['title'] == 'Avs' for r in avs['wiki_systems'])
    assert avs['search']['linked_context_wiki_ids']
    assert len(cat['provenance']['resources']) == 2
    assert all(r['resource'].startswith('package:defenseviz/data/') for r in cat['provenance']['resources'])


@pytest.mark.parametrize('value', ['', ' ', 'NA', 'NaN', 'nan', 'null', 'None', '-'])
def test_producer_missing_fqn_sentinels_are_nonfatal_label_context(catalog, value):
    row = call(); row['native_fields']['model_fqn'] = value
    r = resolution(catalog, row)
    assert r['identity_status'] == 'exact_label_context'
    assert r['native_model_fqn'] == value and r['native_model_fqn_normalized'] is None
    assert 'native_fqn_field_present_but_missing' in r['reason_codes']


@pytest.mark.parametrize('value', [3.1, 1, True, ['model'], {'name': 'model'}])
def test_malformed_optional_fqn_is_nonfatal_and_cannot_fallback(catalog, value):
    row = call(); row['native_fields']['model_fqn'] = value
    r = resolution(catalog, row)
    assert r['identity_status'] == 'unmatched' and not r['definition_ids']
    assert r['native_model_fqn'] == value
    assert 'optional_native_identity_metadata_malformed' in r['reason_codes']


@pytest.mark.parametrize('field,value', [('model_version', 3.1), ('model_version', {}),
                                       ('model_version_source', []), ('model_package', {}),
                                       ('source_version', ['invalid']), ('model_version', '')])
def test_malformed_or_blank_optional_versions_do_not_abort_existing_analysis(catalog, field, value):
    row = call(); row[field] = value
    r = resolution(catalog, row)
    assert r['identity_status'] == 'exact_fqn'
    assert r['version_relation'] == 'version_unavailable'
    assert r['optional_version_metadata_diagnostics']


@pytest.mark.parametrize('field,value', [('model_definition_sha256', ['bad']),
                                       ('model_definition_sha256', 'invalid'),
                                       ('model_package_commit', 3), ('model_package_commit', 'not-a-commit')])
def test_malformed_optional_binding_digest_is_diagnostic_not_fatal(catalog, field, value):
    row = call(model_version='synthetic-1'); row['native_fields'][field] = value
    r = resolution(catalog, row)
    assert r['identity_status'] == 'exact_fqn' and r['version_relation'] == 'version_unavailable'
    assert 'optional_model_version_metadata_malformed' in r['reason_codes']


def test_provider_specific_accounting_is_exhaustive_with_overlapping_links(catalog):
    calls = [call(), call('Future', provider='PADLOC', fqn=False), call('Unknown', provider='PADLOC', fqn=False)]
    out = describe_systems(dict(calls=calls), {}, catalog)
    by = out['accounting']['by_provider']
    assert by['DefenseFinder']['detection_count'] == 1 and by['DefenseFinder']['model_linked_detections'] == 1
    assert by['DefenseFinder']['wiki_linked_detections'] == 1
    assert by['PADLOC']['detection_count'] == 2 and by['PADLOC']['wiki_linked_detections'] == 1
    assert by['PADLOC']['detections_without_metadata'] == 1
    assert sum(sum(r['identity_status_counts'].values()) for r in by.values()) == len(calls)
    assert sum(sum(r['version_relation_counts'].values()) for r in by.values()) == len(calls)


def test_native_and_catalogue_scope_remain_separate_for_label_only_context(catalog):
    r = resolution(catalog, call(fqn=False))
    assert r['identity_status'] == 'exact_label_context'
    assert r['native_declared_defense_scope'] == r['defense_scope'] == 'unresolved'
    assert r['catalogue_definition_defense_scope'] == 'host_defense'
    assert r['catalogue_definition_defense_scope_basis'].startswith('exact_label_context:')


def test_cross_provider_wiki_context_is_explicit_label_coincidence(catalog):
    r = resolution(catalog, call('Future', provider='PADLOC', fqn=False))
    assert r['wiki_identity_context_basis'] == 'cross_provider_exact_label_coincidence_not_equivalence'
    assert r['identity_status'] == 'wiki_context_only' and not r['definition_ids']


def test_unbound_component_name_comparison_is_not_architecture_evidence(catalog):
    r = resolution(catalog, call())
    assert 'component_name_comparison_unbound_to_run_version' in r['reason_codes']
    assert 'native_component_not_listed_in_current_definition' not in r['reason_codes']
    assert 'not_evidence_of_architecture_change' in r['current_definition_component_context']['scope']


def test_package_agreement_exports_source_basis_without_definition_claim(catalog):
    r = resolution(catalog, call(model_version='synthetic-1'))
    assert r['version_relation'] == 'same_declared_package_version'
    assert r['package_version_agreement_scope'] == 'package_declaration_only_not_exact_definition_binding'
    assert 'catalogue_package_version_basis' in r
    assert information_tables(info(catalog, call(model_version='synthetic-1')))['detections'][0]['native_model_version_first_reported'] == 'synthetic-1'


def test_real_review_reconciles_gene_contig_explicit_replicon_and_missing_gene(catalog):
    from defenseviz.analysis.system_review import review_all_systems
    a, b, c = call(), call(), call()
    a.update(gene_id='from-gene', replicon_id=None)
    b.update(gene_id='explicit', replicon_id='explicit-replicon')
    c.update(gene_id='missing-gene', replicon_id=None)
    data = dict(samples=[dict(sample_id='s1')], calls=[a, b, c], genes=[
        dict(sample_id='s1', gene_id='from-gene', contig='gene-contig'),
        dict(sample_id='s1', gene_id='explicit', contig='ignored-contig')])
    result = dict(system_review=review_all_systems(data, {}))
    out = describe_systems(data, result, catalog)
    assert {r['gene_id']:r['replicon_id'] for r in out['detections']} == {
        'from-gene':'gene-contig', 'explicit':'explicit-replicon', 'missing-gene':None}
    assert {r['detection_id'] for r in out['detections']} == {r['detection_id'] for r in result['system_review']['detection_accounting']}


@pytest.mark.parametrize('fqn', ['defense-finder-models/DefenseFinder/Missing/Future', 3.1])
def test_independent_wiki_label_context_never_resolves_an_invalid_model_identity(catalog, fqn):
    row = call('Future', fqn=False); row['native_fields']['model_fqn'] = fqn
    r = resolution(catalog, row)
    assert r['identity_status'] == 'unmatched' and not r['definition_ids']
    assert r['wiki_ids'] == ['defensefinder-wiki:Future']
    assert r['wiki_identity_context_basis'] == 'same_provider_exact_native_label_context'
    assert 'independent_wiki_label_context_does_not_resolve_model_identity' in r['reason_codes']
