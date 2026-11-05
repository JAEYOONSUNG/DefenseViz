"""Pinned official reference context; never a functional evidence classifier.

Wiki pages, source model definitions and saved detections have separate identities.
All functions are pure with respect to caller inputs. Runtime lookup is exact;
forgiving search only returns browsing candidates.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import hashlib
from importlib.resources import files
import json
import re
from .analysis.system_review import declared_defense_scope, detection_identity


RESOURCE_NAMES = ('defensefinder_wiki.json', 'defensefinder_models.json')
_CONTEXT = 'current_source_metadata_not_query_function_validation'
_PROVIDER_PACKAGE_VERSION_SOURCE = 'defense-finder-models'
_PRODUCER_MISSING = {'', 'NA', 'NaN', 'nan', 'null', 'None', '-'}
_RELATED_SYSTEM_RELATIONS = {'group_member', 'related_system_name', 'former_name',
                             'source_alternate_name', 'subtype_related_name', 'subtype_alternate_name'}


def _text(value, where):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{where} must be nonempty text')
    return value


def _rows(value, where):
    if not isinstance(value, list) or any(not isinstance(r, dict) for r in value):
        raise ValueError(f'{where} must be a list of objects')
    return value


def _sha(value, where, size=64):
    if not isinstance(value, str) or not re.fullmatch('[0-9a-f]{' + str(size) + '}', value):
        raise ValueError(f'{where} must be a lowercase {size}-digit digest')
    return value


def _ref(value, where):
    if not isinstance(value, dict):
        raise ValueError(f'{where} must be source metadata')
    _text(value.get('path'), where + '.path')
    _sha(value.get('sha256'), where + '.sha256')
    return value


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)


def _id(prefix, value):
    return prefix + ':' + hashlib.sha256(_canonical(value).encode()).hexdigest()


def _optional_text(value, field, diagnostics):
    """Optional producer metadata may be absent or malformed without aborting analysis."""
    if value is None:
        return None, 'missing'
    if isinstance(value, str) and value.strip() in _PRODUCER_MISSING:
        diagnostics.append(dict(field=field, reason='producer_missing_sentinel', raw_value=value))
        return None, 'missing'
    if not isinstance(value, str):
        diagnostics.append(dict(field=field, reason='malformed_optional_text', raw_value=deepcopy(value)))
        return None, 'malformed'
    return value, 'provided'


def _optional_mapping(value, field, diagnostics):
    if value is None:
        return {}
    if not isinstance(value, dict):
        diagnostics.append(dict(field=field, reason='malformed_optional_mapping', raw_value=deepcopy(value)))
        return {}
    return value


def _unique(rows, key, where):
    result = {}
    for row in rows:
        identity = _text(row.get(key), where + '.' + key)
        if identity in result:
            raise ValueError(f'Duplicate {where} identity: {identity}')
        result[identity] = row
    return result


def _wiki_labels(row):
    return [(row['title'], 'title'), (row['slug'], 'slug')] + [
        (a['value'], 'alias:' + a['source_field']) for a in row.get('aliases', [])]


def build_catalog(wiki, models, *, resources=None):
    """Validate two source snapshots and build exact descriptive-context links."""
    for value in (wiki, models):
        if not isinstance(value, dict) or type(value.get('schema_version')) is not int or value['schema_version'] != 1:
            raise ValueError('Unsupported official catalogue schema_version')
    source = wiki.get('source')
    if not isinstance(source, dict):
        raise ValueError('Wiki source metadata missing')
    _text(source.get('repository'), 'wiki.repository')
    _sha(source.get('commit'), 'wiki.commit', 40)
    wrows = _rows(wiki.get('systems'), 'wiki.systems')
    if source.get('coverage', {}).get('system_pages', len(wrows)) != len(wrows):
        raise ValueError('Wiki source inventory count does not reconcile')
    windex = _unique(wrows, 'id', 'wiki')
    slugs = set()
    for row in wrows:
        slug = _text(row.get('slug'), 'wiki.slug')
        if slug in slugs or row['id'] != 'defensefinder-wiki:' + slug:
            raise ValueError('Duplicate or inconsistent wiki slug identity')
        slugs.add(slug)
        _text(row.get('title'), 'wiki.title'); _ref(row.get('source'), 'wiki.source')
        if not isinstance(row.get('facts', {}), dict):
            raise ValueError('Wiki facts must be source-specific objects')
        for fact in row.get('facts', {}).values():
            if not isinstance(fact, dict) or fact.get('state') not in {'reported', 'unknown', 'not_provided'}:
                raise ValueError('Invalid wiki fact state')
        for alias in _rows(row.get('aliases', []), 'wiki.aliases'):
            _text(alias.get('value'), 'wiki.alias.value')
            if alias.get('source_field') != 'refseqSystemName' or type(alias.get('source_line')) is not int or alias['source_line'] < 1:
                raise ValueError('Wiki alias needs exact source-field provenance')
        for related in _rows(row.get('related_names', []), 'wiki.related_names'):
            _text(related.get('value'), 'related name.value')
            if not isinstance(related.get('source_lines'), list) or not related['source_lines'] or any(type(n) is not int or n < 1 for n in related['source_lines']):
                raise ValueError('Related name needs exact source-line provenance')
    sources = _rows(models.get('sources'), 'models.sources')
    sindex = _unique(sources, 'source_id', 'model source')
    for row in sources:
        _text(row.get('repository'), 'model source.repository')
        _sha(row.get('commit'), 'model source.commit', 40)
        _text(row.get('package_version'), 'model source.package_version')
    defs = _rows(models.get('definitions'), 'model definitions')
    coverage = models.get('coverage', {})
    if coverage.get('definition_count', len(defs)) != len(defs):
        raise ValueError('Model source definition count does not reconcile')
    by_source = dict(Counter(row.get('source_id') for row in defs))
    if coverage.get('source_definition_counts', by_source) != by_source:
        raise ValueError('Model per-source definition counts do not reconcile')
    configurations = _rows(models.get('package_configurations', []), 'package_configurations')
    if (coverage.get('package_configuration_count', len(configurations)) != len(configurations)
            or coverage.get('xml_file_count', len(defs) + len(configurations)) != len(defs) + len(configurations)):
        raise ValueError('Model configuration/XML inventory does not reconcile')
    mindex = _unique(defs, 'definition_id', 'model')
    aliases = {}
    for row in defs:
        sid = row.get('source_id')
        if sid not in sindex:
            raise ValueError('Model definition has unknown source_id')
        key = _text(row.get('source_definition_key'), 'model.source_definition_key')
        if row['definition_id'] != sid + '/' + key or key.startswith('/') or '..' in key.split('/'):
            raise ValueError('Model identity is not its exact source-qualified path')
        if row.get('package_version', sindex[sid]['package_version']) != sindex[sid]['package_version']:
            raise ValueError('Model package version contradicts source metadata')
        _text(row.get('model_name'), 'model.model_name')
        _text(row.get('source_namespace'), 'model.source_namespace')
        if row['model_name'] != key.split('/')[-1]:
            raise ValueError('Model name contradicts exact definition path')
        namespace = 'CasFinder' if sid == 'CasFinder' else key.split('/')[0]
        if row['source_namespace'] != namespace:
            raise ValueError('Model namespace contradicts exact definition path')
        if row.get('defense_scope') not in {'host_defense', 'anti_defense', 'unresolved'}:
            raise ValueError('Invalid source defense scope')
        if row['source_namespace'] == 'ADF' and row['defense_scope'] != 'anti_defense':
            raise ValueError('ADF source scope conflict')
        if namespace in {'DefenseFinder', 'RM', 'CasFinder'} and row['defense_scope'] != 'host_defense':
            raise ValueError('Host-defense source namespace scope conflict')
        _ref(row.get('source_ref'), 'model.source_ref')
        if row['source_ref']['path'] != 'definitions/' + key + '.xml':
            raise ValueError('Definition source path contradicts model identity')
        paths = set()
        for component in _rows(row.get('components'), 'model.components'):
            path = _text(component.get('xpath'), 'component.xpath')
            _text(component.get('native_component_id'), 'component.native_component_id')
            if path in paths:
                raise ValueError('Duplicate component XPath within model')
            paths.add(path)
        for alias in _rows(row.get('exact_fqn_aliases', []), 'model.exact_fqn_aliases'):
            name = _text(alias.get('value'), 'model alias.value')
            _ref(alias.get('source_ref'), 'model alias.source_ref')
            if name in mindex or name in aliases:
                raise ValueError('Ambiguous or shadowing exact FQN alias')
            aliases[name] = row['definition_id']
    labels = defaultdict(set)
    for row in wrows:
        for label, field in _wiki_labels(row):
            labels[label].add((row['id'], field))
    links = []
    for row in defs:
        for field in ('model_name', 'source_family'):
            value = row.get(field)
            for wid, wiki_field in sorted(labels.get(value, ())):
                link = dict(definition_id=row['definition_id'], wiki_id=wid,
                            relationship='parent_family_context' if field == 'source_family' else 'same_declared_name_context',
                            model_field=field, wiki_field=wiki_field, exact_value=value,
                            model_source_ref=deepcopy(row['source_ref']), wiki_source_ref=deepcopy(windex[wid]['source']))
                link['link_id'] = _id('catalog_link', link)
                links.append(link)
    result = dict(schema_version=1, catalogue_kind='official_reference_context',
                  wiki_source=deepcopy(source), model_sources=deepcopy(sources),
                  wiki_systems=deepcopy(wrows), model_definitions=deepcopy(defs), links=links,
                  package_configurations=deepcopy(models.get('package_configurations', [])),
                  source_inventory_coverage=deepcopy(coverage),
                  provenance=dict(resources=deepcopy(resources or []), context_scope=_CONTEXT,
                                  query_function_inferred=False),
                  policy=dict(runtime_identity='exact_source_qualified', search='browsing_candidates_only',
                              row_counts='source_records_not_unique_biological_systems',
                              scientific_membership_changed=False))
    result['coverage'] = dict(wiki_pages=len(wrows), model_definitions=len(defs),
                              package_configurations=len(result['package_configurations']),
                              model_definitions_by_source=dict(Counter(r['source_id'] for r in defs)),
                              wiki_pages_linked_to_models=len({r['wiki_id'] for r in links}),
                              model_definitions_linked_to_wiki=len({r['definition_id'] for r in links}))
    return result


def _json(raw):
    def pairs(values):
        obj = {}
        for k, v in values:
            if k in obj:
                raise ValueError('Duplicate catalogue JSON key: ' + k)
            obj[k] = v
        return obj
    def bad(value):
        raise ValueError('Nonfinite catalogue JSON value: ' + value)
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=bad)


def load_catalog():
    """Load and validate bundled snapshots without network or repository access."""
    payloads, resources = [], []
    for name in RESOURCE_NAMES:
        raw = files('defenseviz').joinpath('data', name).read_bytes()
        payloads.append(_json(raw))
        resources.append(dict(resource='package:defenseviz/data/' + name,
                              sha256=hashlib.sha256(raw).hexdigest()))
    return build_catalog(*payloads, resources=resources)


def _validated(catalog):
    if catalog is None:
        return load_catalog()
    if not isinstance(catalog, dict) or type(catalog.get('schema_version')) is not int or catalog['schema_version'] != 1:
        raise ValueError('Invalid combined catalogue')
    # Rebuild derived context, then reject caller-supplied contradictions.
    rebuilt = build_catalog(dict(schema_version=1, source=catalog.get('wiki_source'), systems=catalog.get('wiki_systems')),
                         dict(schema_version=1, sources=catalog.get('model_sources'), definitions=catalog.get('model_definitions'),
                              package_configurations=catalog.get('package_configurations', []),
                              coverage=catalog.get('source_inventory_coverage', {})),
                         resources=catalog.get('provenance', {}).get('resources', []))
    if catalog.get('links') != rebuilt['links'] or catalog.get('coverage') != rebuilt['coverage']:
        raise ValueError('Combined catalogue links or coverage contradict source records')
    return rebuilt


def search_catalog(query=None, catalog=None):
    """Browse source records; case-insensitive search never assigns identities."""
    catalog = _validated(catalog)
    if query is None:
        return catalog
    query = _text(query, 'catalogue query')
    needle = query.casefold()
    matched_wiki = [r for r in catalog['wiki_systems'] if query == r['id']]
    matched_models = [r for r in catalog['model_definitions'] if query == r['definition_id']]
    exact = bool(matched_wiki or matched_models)
    if not exact:
        matched_wiki = [r for r in catalog['wiki_systems'] if needle in _canonical(r).casefold()]
        matched_models = [r for r in catalog['model_definitions'] if needle in _canonical(r).casefold()]
    wids, mids = {r['id'] for r in matched_wiki}, {r['definition_id'] for r in matched_models}
    # Expand exactly one relationship from direct hits. A context-only parent
    # page must not pull in all sibling subtypes of a specifically queried model.
    context_wids = {r['wiki_id'] for r in catalog['links'] if r['definition_id'] in mids} - wids
    context_mids = {r['definition_id'] for r in catalog['links'] if r['wiki_id'] in wids} - mids
    result = deepcopy(catalog)
    result['wiki_systems'] = [r for r in catalog['wiki_systems'] if r['id'] in wids | context_wids]
    result['model_definitions'] = [r for r in catalog['model_definitions'] if r['definition_id'] in mids | context_mids]
    result['links'] = [r for r in catalog['links'] if r['wiki_id'] in wids | context_wids and r['definition_id'] in mids | context_mids]
    result['search'] = dict(query=query, mode='exact_record_id' if exact else 'case_insensitive_substring_candidates_only',
                            matched_wiki_pages=len(wids), matched_model_definitions=len(mids),
                            direct_wiki_ids=sorted(wids), direct_definition_ids=sorted(mids),
                            linked_context_wiki_ids=sorted(context_wids), linked_context_definition_ids=sorted(context_mids),
                            total_wiki_pages=len(catalog['wiki_systems']), total_model_definitions=len(catalog['model_definitions']))
    return result


def _versions(call, model, source):
    diagnostics = []
    native = _optional_mapping(call.get('native_fields'), 'native_fields', diagnostics)
    sv = _optional_mapping(call.get('source_version'), 'source_version', diagnostics)
    origin, _ = _optional_text(call.get('model_version_source', 'canonical_model_version'), 'model_version_source', diagnostics)
    origin = origin or 'unavailable'
    fallback = origin.endswith('database_version') or origin.endswith('db_version')
    package_claims = set()
    for field, value in [('model_package', call.get('model_package')), ('native_fields.model_package', native.get('model_package')), ('source_version.models_package', sv.get('models_package'))]:
        value, _ = _optional_text(value, field, diagnostics)
        if value is not None:package_claims.add(value)
    package_conflict = bool(package_claims - {model['source_id']})
    scoped = package_claims == {model['source_id']}
    # This saved module field names the DF package, not every package loaded by
    # DefenseFinder. Agreement is package-level, never a definition hash binding.
    provider_wide = model['source_id'] != _PROVIDER_PACKAGE_VERSION_SOURCE and not scoped
    version_reports, ignored = [], []
    for field, value in [('canonical.model_version', None if fallback else call.get('model_version')),
                         ('source_version.models_version', sv.get('models_version')),
                         ('native_fields.model_version', native.get('model_version'))]:
        value, _ = _optional_text(value, field, diagnostics)
        if value is not None:
            if provider_wide and (field == 'source_version.models_version' or field == 'canonical.model_version' and origin == 'saved_module_metadata.models_version'):
                ignored.append(field)
                continue
            version_reports.append(dict(source_field=field, value=value))
    reported = version_reports[0]['value'] if version_reports else None
    # A component HMM model_sha256 is not the system XML definition hash.
    definition_hash, _ = _optional_text(native.get('model_definition_sha256'), 'native_fields.model_definition_sha256', diagnostics)
    package_commit, _ = _optional_text(native.get('model_package_commit'), 'native_fields.model_package_commit', diagnostics)
    for field, value, size in [('model_definition_sha256', definition_hash, 64), ('model_package_commit', package_commit, 40)]:
        if value is not None and not re.fullmatch('[0-9a-f]{' + str(size) + '}', value):
            diagnostics.append(dict(field='native_fields.' + field, reason='malformed_optional_digest', raw_value=value))
            if field == 'model_definition_sha256':definition_hash = None
            else:package_commit = None
    checks = []
    for item in version_reports:
        checks.append(('same_declared_package_version', item['value'] == source['package_version']))
    if definition_hash is not None:
        checks.append(('same_definition_hash', definition_hash == model['source_ref']['sha256']))
    if package_commit is not None:
        checks.append(('same_source_commit', package_commit == source['commit']))
    if any(d['reason'].startswith('malformed_') for d in diagnostics):
        relation = 'version_unavailable'
    elif package_conflict or len({v['value'] for v in version_reports}) > 1 or (any(v for _, v in checks) and not all(v for _, v in checks)):
        relation = 'provenance_conflict'
    elif checks and not any(v for _, v in checks):
        relation = 'different_version'
    elif checks:
        relation = next(label for label in ('same_definition_hash', 'same_source_commit', 'same_declared_package_version') if (label, True) in checks)
    else:
        relation = 'version_unavailable'
    return dict(version_relation=relation, native_model_version=reported,
                optional_version_metadata_diagnostics=diagnostics,
                native_model_versions_reported=version_reports,
                ignored_unscoped_version_fields=ignored, model_package_claims=sorted(package_claims),
                model_version_scope='explicit_package' if scoped else 'provider_package' if not provider_wide else 'per_record_declaration_only',
                native_model_version_source=origin, native_database_version=call.get('database_version', sv.get('db_version')),
                native_software_version=call.get('software_version', sv.get('software_version')),
                native_model_definition_sha256=definition_hash, native_model_package_commit=package_commit,
                catalogue_package_version=source['package_version'], catalogue_source_commit=source['commit'],
                catalogue_package_version_basis=deepcopy(source.get('version_basis')),
                package_version_agreement_scope='package_declaration_only_not_exact_definition_binding',
                catalogue_definition_sha256=model['source_ref']['sha256'])


def _indexes(catalog):
    models = {r['definition_id']: r for r in catalog['model_definitions']}
    names, aliases, wiki, related = defaultdict(list), {}, defaultdict(set), defaultdict(list)
    for r in models.values():
        names[r['model_name']].append(r['definition_id'])
        for a in r.get('exact_fqn_aliases', []):aliases[a['value']] = r['definition_id']
    for r in catalog['wiki_systems']:
        for label, _ in _wiki_labels(r):wiki[label].add(r['id'])
        for relationship in r.get('related_names', []):
            if relationship.get('relationship') in _RELATED_SYSTEM_RELATIONS and relationship.get('scope') == 'descriptive_name_relationship_not_runtime_alias':
                related[relationship['value']].append(dict(wiki_id=r['id'], **deepcopy(relationship), source_ref=deepcopy(r['source'])))
    links = defaultdict(set)
    for r in catalog['links']:links[r['definition_id']].add(r['wiki_id'])
    return models, names, aliases, wiki, links, {r['source_id']: r for r in catalog['model_sources']}, related


def _resolve(call, indexes):
    models, names, aliases, wiki, links, sources, related = indexes
    diagnostics = []
    native = _optional_mapping(call.get('native_fields'), 'native_fields', diagnostics)
    provider = call.get('caller')
    raw_family = call.get('native_family', call.get('family'))
    family, _ = _optional_text(raw_family, 'native_family', diagnostics)
    raw_fqn = native.get('model_fqn', call.get('model_fqn'))
    fqn, fqn_state = _optional_text(raw_fqn, 'native_fields.model_fqn', diagnostics)
    raw_namespace = call.get('native_namespace', native.get('namespace'))
    namespace, _ = _optional_text(raw_namespace, 'native_namespace', diagnostics)
    scope_call = dict(call, native_fields=dict(native, model_fqn=fqn, namespace=namespace), native_namespace=namespace)
    scope, scope_reason = declared_defense_scope(scope_call)
    result = dict(native_family=raw_family, native_component_model=call.get('native_model'), native_model_fqn=raw_fqn,
                  native_model_fqn_normalized=fqn, native_namespace=raw_namespace,
                  native_version_claims=dict(model_version=call.get('model_version'), model_version_source=call.get('model_version_source'),
                                             database_version=call.get('database_version'), software_version=call.get('software_version'),
                                             source_version=deepcopy(call.get('source_version') or {})),
                  identity_status='unmatched', version_relation='not_applicable',
                  definition_ids=[], candidate_definition_ids=[], wiki_ids=sorted(wiki.get(family, ())),
                  wiki_identity_context_basis='cross_provider_exact_label_coincidence_not_equivalence' if provider != 'DefenseFinder' else 'same_provider_exact_native_label_context',
                  related_wiki_context=deepcopy(related.get(family, [])),
                  defense_scope=scope, native_declared_defense_scope=scope, native_declared_defense_scope_basis=scope_reason,
                  catalogue_definition_defense_scope=None, catalogue_definition_defense_scope_basis=None,
                  defense_scope_basis='native_declared_namespace_separate_from_catalogue_context',
                  optional_identity_metadata_diagnostics=diagnostics,
                  reason_codes=[], context_scope=_CONTEXT, query_function_inferred=False)
    observations = call.get('native_observations') or []
    if not isinstance(observations, list):
        diagnostics.append(dict(field='native_observations', reason='malformed_optional_observations', raw_value=deepcopy(observations)))
        observations = []
    observed = defaultdict(set)
    for index, row in enumerate([native, *observations]):
        row = _optional_mapping(row, 'native_observations[' + str(index) + ']', diagnostics)
        for field in ('model_fqn', 'subtype', 'namespace'):
            value, _ = _optional_text(row.get(field), 'native_observations[' + str(index) + '].' + field, diagnostics)
            if value is not None:observed[field].add(value)
    conflicts = [field for field, values in observed.items() if len(values) > 1]
    canonical_fqn, _ = _optional_text(call.get('model_fqn'), 'canonical.model_fqn', diagnostics)
    if canonical_fqn is not None and fqn is not None and canonical_fqn != fqn:
        conflicts.append('canonical.model_fqn')
    if fqn_state == 'malformed' or any(d['reason'].startswith('malformed_') for d in diagnostics):
        result['reason_codes'].append('optional_native_identity_metadata_malformed')
        if result['wiki_ids']:
            result['reason_codes'].append('independent_wiki_label_context_does_not_resolve_model_identity')
        return result
    if any(d['reason'] == 'producer_missing_sentinel' and 'model_fqn' in d['field'] for d in diagnostics):
        result['reason_codes'].append('native_fqn_field_present_but_missing')
    if conflicts:
        result.update(identity_status='conflicting_identity', conflicting_native_fields=conflicts)
        result['reason_codes'].append('inconsistent_native_identity_observations')
        return result
    if provider != 'DefenseFinder':
        result['identity_status'] = 'provider_out_of_scope'
        result['reason_codes'].append('no_model_catalogue_for_provider')
        result['wiki_ids'] = sorted(wiki.get(family, ()))
        if result['wiki_ids']:result['identity_status'] = 'wiki_context_only'
        return result
    chosen = None
    if fqn is not None:
        chosen = fqn if fqn in models else aliases.get(fqn)
        if chosen is None:
            result['reason_codes'].append('explicit_model_identity_unresolved')
            if result['wiki_ids']:
                result['reason_codes'].append('independent_wiki_label_context_does_not_resolve_model_identity')
            if fqn.split('/')[0] not in sources:
                result['identity_status'] = 'provider_out_of_scope'
            return result
        result['identity_status'] = 'exact_fqn' if chosen == fqn else 'documented_alias'
        model = models[chosen]
        subtype, _ = _optional_text(native.get('subtype'), 'native_fields.subtype', diagnostics)
        allowed_family = {model['model_name'], model.get('source_family')}
        if (subtype is not None and subtype != model['model_name']) or (family is not None and family not in allowed_family):
            result.update(identity_status='conflicting_identity', candidate_definition_ids=[chosen])
            result['reason_codes'].append('native_label_contradicts_explicit_model_identity')
            return result
    else:
        candidates = names.get(family, [])
        result['candidate_definition_ids'] = sorted(candidates)
        if len(candidates) > 1:
            result['identity_status'] = 'ambiguous'
            result['reason_codes'].append('native_label_has_multiple_model_definitions')
        elif len(candidates) == 1:
            chosen = candidates[0]; result['identity_status'] = 'exact_label_context'
            result['reason_codes'].append('native_fqn_unavailable_label_context_only')
        else:
            result['reason_codes'].append('not_in_catalogue_snapshot')
        result['wiki_ids'] = sorted(wiki.get(family, ()))
        if chosen is None and result['wiki_ids'] and not candidates:
            result['identity_status'] = 'wiki_context_only'
    if chosen is not None:
        model = models[chosen]
        if namespace is not None and namespace != model['source_namespace']:
            result.update(identity_status='conflicting_identity', candidate_definition_ids=[chosen], definition_ids=[])
            result['reason_codes'].append('native_namespace_contradicts_model_definition')
            return result
        result['definition_ids'] = [chosen]
        label_wiki = bool(result['wiki_ids'])
        result['wiki_ids'] = sorted(set(result['wiki_ids']) | links[chosen])
        if links[chosen]:
            result['wiki_identity_context_basis'] = 'exact_label_and_current_model_source_context' if label_wiki else 'current_model_source_name_context'
        result['catalogue_definition_defense_scope'] = model['defense_scope']
        result['catalogue_definition_defense_scope_basis'] = result['identity_status'] + ':exact_source_definition_namespace'
        result.update(_versions(call, model, sources[model['source_id']]))
        declared = {r['native_component_id'] for r in model['components']}
        observed_fields = {}
        for field, value in [('native_model', call.get('native_model')), ('native_fields.hit_gene_ref', native.get('hit_gene_ref'))]:
            value, _ = _optional_text(value, field, diagnostics)
            if value is not None:observed_fields[field] = value
        observed = set(observed_fields.values())
        listed = sorted(observed & declared)
        unlisted = sorted(observed - declared)
        result['current_definition_component_context'] = dict(
            state='partial_declared_component_names' if listed and unlisted else 'exact_declared_component_name' if listed else 'not_listed_in_current_definition' if observed else 'component_identity_unavailable',
            observed_fields=observed_fields, observed_names=sorted(observed), exact_declared_names=listed, not_declared_names=unlisted,
            saved_run_version_relation=result['version_relation'],
            scope='current_snapshot_names_only_not_evidence_of_architecture_change_or_role_validation')
        if unlisted:
            result['reason_codes'].append('component_name_comparison_unbound_to_run_version' if result['version_relation'] in {'different_version', 'version_unavailable', 'provenance_conflict'} else 'native_component_not_listed_in_current_definition')
        if any(d['reason'].startswith('malformed_') for d in result['optional_version_metadata_diagnostics']):
            result['reason_codes'].append('optional_model_version_metadata_malformed')
        if result['ignored_unscoped_version_fields']:
            result['reason_codes'].append('version_scope_unbound')
        if result['version_relation'] in {'different_version', 'version_unavailable', 'provenance_conflict'}:
            result['reason_codes'].append('current_definition_not_bound_to_saved_run_version')
    return result


def describe_systems(data, result, catalog=None):
    """Describe every saved detection and review unit without changing analysis."""
    catalog = _validated(catalog); indexes = _indexes(catalog)
    wiki_by_id = {r['id']: r for r in catalog['wiki_systems']}
    calls = _rows(data.get('calls', []), 'calls')
    review = result.get('system_review')
    units, expected = {}, {}
    if review is not None:
        for unit in _rows(review.get('units'), 'system_review.units'):
            uid = unit['unit_id']
            if uid in units:raise ValueError('Duplicate system review unit')
            units[uid] = {k: deepcopy(unit.get(k)) for k in ('unit_id', 'unit_kind', 'sample_id', 'provider', 'run_id', 'replicon_id', 'native_system_id', 'member_keys', 'selected_member_keys', 'native_families', 'native_models')}
            units[uid]['detection_ids'] = list(unit['detection_ids'])
            for did in unit['detection_ids']:
                if did in expected:raise ValueError('Detection belongs to multiple primary review units')
                expected[did] = uid
    genes = {(g['sample_id'], g['gene_id']): g for g in data.get('genes', [])}
    detections, seen = [], set()
    for number, call in enumerate(calls):
        key = (call['sample_id'], call['gene_id'])
        replicon = call.get('replicon_id')
        if replicon is None:replicon = genes.get(key, {}).get('contig')
        identity = detection_identity(call, replicon)
        if review is not None:
            did = _id('detection', identity)
            if did not in expected or did in seen:
                raise ValueError('Catalogue detection identity does not reconcile with system review')
            uid = expected[did]
        else:
            did = _id('source_annotation', [number, identity]); uid = did
            units[uid] = dict(unit_id=uid, unit_kind='ungrouped_source_annotation', sample_id=key[0], provider=call.get('caller'),
                              run_id=call.get('run_id'), replicon_id=replicon, native_system_id=None,
                              reported_native_system_id=call.get('system_id'), member_keys=[list(key)],
                              selected_member_keys=[list(key)] if call.get('selected', True) else [],
                              native_families=[call.get('native_family', call.get('family'))], native_models=[call.get('native_model')], detection_ids=[did])
        seen.add(did)
        row = dict(detection_id=did, unit_id=uid, sample_id=key[0], gene_id=key[1], provider=call.get('caller'),
                   run_id=call.get('run_id'), replicon_id=replicon, native_system_id=call.get('system_id'),
                   evidence_kind=call.get('evidence_kind', 'final_member'), selected=call.get('selected', True),
                   source_refs=deepcopy(call.get('source_refs') or ([call['source_ref']] if call.get('source_ref') else [])))
        row.update(_resolve(call, indexes)); detections.append(row)
    if review is not None and seen != set(expected):
        raise ValueError('System review contains unaccounted catalogue detections')
    by_unit = defaultdict(list)
    for row in detections:by_unit[row['unit_id']].append(row)
    for uid, unit in units.items():
        rows = by_unit[uid]
        unit.update(coding_sequences=len({tuple(k) for k in unit.get('member_keys') or []}),
                    selected_coding_sequences=len({tuple(k) for k in unit.get('selected_member_keys') or []}),
                    definition_ids=sorted({x for r in rows for x in r['definition_ids']}),
                    wiki_ids=sorted({x for r in rows for x in r['wiki_ids']}),
                    identity_statuses=sorted({r['identity_status'] for r in rows}) or ['no_native_model_identity'],
                    version_relations=sorted({r['version_relation'] for r in rows}) or ['not_applicable'],
                    native_declared_defense_scopes=sorted({r['native_declared_defense_scope'] for r in rows}),
                    catalogue_definition_defense_scopes=sorted({r['catalogue_definition_defense_scope'] for r in rows if r['catalogue_definition_defense_scope'] is not None}),
                    reason_codes=sorted({x for r in rows for x in r['reason_codes']}) or ([] if rows else ['no_native_detection']),
                    context_scope=_CONTEXT, query_function_inferred=False)
        unit['reference_titles'] = [wiki_by_id[wid]['title'] for wid in unit['wiki_ids']]
        unit['reference_mechanisms'] = [dict(wiki_id=wid,
                                            summary=wiki_by_id[wid].get('mechanism_summary'),
                                            context=deepcopy(wiki_by_id[wid].get('mechanism_context')))
                                        for wid in unit['wiki_ids']]
        related_context = {_canonical(r): r for row in rows for r in row['related_wiki_context']}
        unit['related_wiki_context'] = [deepcopy(related_context[k]) for k in sorted(related_context)]
        unit['related_reference_mechanisms'] = [dict(wiki_id=r['wiki_id'], relationship=r['relationship'],
                                                    summary=wiki_by_id[r['wiki_id']].get('mechanism_summary'),
                                                    context=deepcopy(wiki_by_id[r['wiki_id']].get('mechanism_context')))
                                                for r in unit['related_wiki_context']]
    mids = {x for r in detections for x in r['definition_ids']}
    wids = {x for r in detections for x in r['wiki_ids']} | {x['wiki_id'] for r in detections for x in r['related_wiki_context']}
    accounting = dict(input_call_count=len(calls), detection_count=len(detections), unit_count=len(units),
                      detections_without_metadata=sum(not r['definition_ids'] and not r['wiki_ids'] and not r['related_wiki_context'] for r in detections),
                      model_linked_detections=sum(bool(r['definition_ids']) for r in detections),
                      wiki_linked_detections=sum(bool(r['wiki_ids']) for r in detections),
                      related_wiki_context_detections=sum(bool(r['related_wiki_context']) for r in detections),
                      identity_status_counts=dict(Counter(r['identity_status'] for r in detections)),
                      version_relation_counts=dict(Counter(r['version_relation'] for r in detections)),
                      provider_counts=dict(Counter(r['provider'] for r in detections)), unaccounted_detections=0)
    accounting['link_count_semantics'] = 'model_wiki_and_related_context_counts_overlap_not_recovery_metrics'
    accounting['by_provider'] = {}
    for provider in sorted({r['provider'] for r in detections}, key=str):
        rows = [r for r in detections if r['provider'] == provider]
        accounting['by_provider'][provider] = dict(
            detection_count=len(rows), identity_status_counts=dict(Counter(r['identity_status'] for r in rows)),
            version_relation_counts=dict(Counter(r['version_relation'] for r in rows)),
            model_linked_detections=sum(bool(r['definition_ids']) for r in rows),
            wiki_linked_detections=sum(bool(r['wiki_ids']) for r in rows),
            related_wiki_context_detections=sum(bool(r['related_wiki_context']) for r in rows),
            detections_without_metadata=sum(not r['definition_ids'] and not r['wiki_ids'] and not r['related_wiki_context'] for r in rows))
    return dict(schema_version=1, operation='descriptive_context_for_saved_annotations',
                accounting=accounting, provenance=deepcopy(catalog['provenance']),
                units=list(units.values()), detections=detections,
                wiki_systems=[deepcopy(r) for r in catalog['wiki_systems'] if r['id'] in wids],
                model_definitions=[deepcopy(r) for r in catalog['model_definitions'] if r['definition_id'] in mids],
                links=[deepcopy(r) for r in catalog['links'] if r['wiki_id'] in wids and r['definition_id'] in mids],
                sources=[dict(deepcopy(catalog['wiki_source']), source_id='defensefinder-wiki'), *deepcopy(catalog['model_sources'])],
                scientific_membership_changed=False, query_function_inferred=False)


def catalog_tables(catalog):
    """Flat record projections retaining source-specific structured facts."""
    rows = dict(wiki_systems=[], model_definitions=[], model_components=[],
                links=deepcopy(catalog['links']), sources=[dict(deepcopy(catalog['wiki_source']), source_id='defensefinder-wiki'), *deepcopy(catalog['model_sources'])],
                package_configurations=deepcopy(catalog.get('package_configurations', [])))
    for system in catalog['wiki_systems']:
        row = deepcopy(system)
        for field in ('sensor', 'activator', 'effector', 'pfam'):
            fact = system.get('facts', {}).get(field, {})
            row[field] = fact.get('value')
            row[field + '_state'] = fact.get('state', 'not_provided')
        row.update(source_mechanism_statement_state=system.get('mechanism_context', {}).get('state', 'not_provided'),
                   citation_dois='; '.join(r['doi'] for r in system.get('references', []) if r.get('doi')),
                   wiki_url=system['source'].get('website_url'), source_url=system['source'].get('url'),
                   source_path=system['source']['path'], source_sha256=system['source']['sha256'],
                   source_commit=catalog['wiki_source']['commit'],
                   reported_example_count=len(system.get('protein_examples', [])))
        rows['wiki_systems'].append(row)
    sources = {r['source_id']: r for r in catalog['model_sources']}
    for model in catalog['model_definitions']:
        row = {k: deepcopy(v) for k, v in model.items() if k != 'components'}
        row.update(source_url=model['source_ref'].get('url'), source_path=model['source_ref']['path'],
                   source_sha256=model['source_ref']['sha256'], source_commit=sources[model['source_id']]['commit'],
                   component_definition_count=len(model['components']))
        rows['model_definitions'].append(row)
        for component in model['components']:
            rows['model_components'].append(dict(definition_id=model['definition_id'], source_id=model['source_id'],
                                                 model_name=model['model_name'], source_path=model['source_ref']['path'],
                                                 source_sha256=model['source_ref']['sha256'], source_commit=sources[model['source_id']]['commit'],
                                                 **deepcopy(component)))
    return rows


def information_tables(info):
    result = {name: deepcopy(info[name]) for name in ('units', 'detections', 'sources', 'wiki_systems', 'model_definitions', 'links')}
    wiki = {r['id']: r for r in info['wiki_systems']}
    models = {r['definition_id']: r for r in info['model_definitions']}
    for row in result['units']:
        primary = [wiki[wid] for wid in row['wiki_ids']]
        related = [wiki[wid] for wid in sorted({r['wiki_id'] for r in row['related_wiki_context']})]
        row['reference_titles'] = '; '.join(r['title'] for r in primary)
        row['reference_mechanism_summaries'] = '; '.join(r['title'] + ': ' + r['mechanism_summary'] for r in primary if r.get('mechanism_summary'))
        row['reference_source_mechanism_statement_states'] = '; '.join(r['title'] + ': ' + r.get('mechanism_context', {}).get('state', 'not_provided') for r in primary)
        row['reference_wiki_urls'] = '; '.join(r['source']['website_url'] for r in primary if r['source'].get('website_url'))
        row['reference_source_urls'] = '; '.join(sorted({r['source']['url'] for r in primary if r['source'].get('url')} |
                                                       {models[mid]['source_ref']['url'] for mid in row['definition_ids'] if models[mid]['source_ref'].get('url')}))
        row['related_reference_titles'] = '; '.join(r['title'] for r in related)
        row['related_mechanism_summaries'] = '; '.join(r['title'] + ': ' + r['mechanism_summary'] for r in related if r.get('mechanism_summary'))
        row['related_reference_wiki_urls'] = '; '.join(r['source']['website_url'] for r in related if r['source'].get('website_url'))
    # Reuse the human-readable source projections without reclassifying records.
    wiki_source = next(r for r in info['sources'] if r['source_id'] == 'defensefinder-wiki')
    projected = catalog_tables(dict(wiki_systems=info['wiki_systems'], wiki_source=wiki_source,
                                    model_definitions=info['model_definitions'],
                                    model_sources=[r for r in info['sources'] if r['source_id'] != 'defensefinder-wiki'], links=info['links']))
    result.update({k: projected[k] for k in ('wiki_systems', 'model_definitions', 'model_components')})
    for row in result['detections']:
        row['native_model_version_first_reported'] = row.get('native_model_version')
    return result
