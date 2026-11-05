"""Frozen cognate Type-II divergent-neighborhood evidence, with source bounds.

This is the selected 800-bp branch of DNMB's existing context rule. It does not
infer an operon, recognition specificity, catalytic function or activity.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import math
import hashlib
import re

from .analysis.policy import Policy

RULE_VERSION = "selected_cognate_divergent_800_v1"
RULE_PROVENANCE = {'policy/source/divergent_context.py': '1be7ba45272e71e32f19ea3d944314e17cb8e34fe9e8cc459fd853a8c3420481', 'policy/build_fixed_context.py': 'c6ee7bf8c27f882e7e590e1124f73c809a23fca19ab05d651de0320811e0bcfc', 'DNMB/db_module_rebasefinder.R': 'b054397e6951d05c9f1fb38a2c4630538fc481ba19edce62de6839fe54355997'}
MAX_INTERVENING_GENES = 1
RM_TYPES = {"Type I", "Type II", "Type III", "Type IV"}
MISSING = {"", "NA", "NaN", "nan", "None", "null", "-"}


def _value(value):
    return None if value is None or isinstance(value, str) and value.strip() in MISSING else value


def _number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _bool(value):
    if type(value) is bool:
        return value
    if str(value).lower() in {"true", "1"}:
        return True
    if str(value).lower() in {"false", "0"}:
        return False
    return None


def _and(*values):
    return False if any(v is False for v in values) else None if any(v is None for v in values) else True


def _refs(rows):
    out = []
    for row in rows:
        for ref in [row.get("source_ref"), *row.get("source_refs", [])]:
            if isinstance(ref, dict) and ref.get("path") and re.fullmatch(r"[a-f0-9]{64}", str(ref.get("sha256", ""))):
                if ref not in out:
                    out.append(dict(ref))
    return out


def cognate_subjects(r_subject, m_subject):
    """Exact frozen relationship between database records, not an allowlist."""
    if not isinstance(r_subject, str) or not isinstance(m_subject, str) or not _value(r_subject) or not _value(m_subject):
        return False
    r_stem = r_subject[2:] if r_subject.startswith("R.") else r_subject
    return m_subject.startswith("M.") and m_subject[2:] == r_stem


def divergent_geometry(left, right, intervening, gap_max=None):
    """Inclusive CDS endpoints; exactly head-to-head and at most one intervening gene."""
    gap_max = Policy().max_divergent_gap_bp if gap_max is None else gap_max
    if gap_max != Policy().max_divergent_gap_bp:
        raise ValueError("Automatic neighborhood evidence requires the frozen 800-bp gap")
    if not 0 <= intervening <= MAX_INTERVENING_GENES:
        return False
    if left.get("contig") != right.get("contig"):
        return False
    if left.get("direction", left.get("strand")) != "-" or right.get("direction", right.get("strand")) != "+":
        return False
    bounds = [_number(r.get(k)) for r in (left, right) for k in ("start", "end")]
    if any(v is None for v in bounds):
        return False
    a, b, c, d = bounds
    return a <= b < c <= d and c - b - 1 <= gap_max


def _inventory(data, sid):
    """Prove the complete imported gene inventory against a pinned annotation table."""
    genes = [g for g in data.get("genes", []) if g.get("sample_id") == sid]
    ids = [g.get("gene_id") for g in genes]
    valid = bool(genes) and all(ids) and len(ids) == len(set(ids))
    for gene in genes:
        a, b = _number(gene.get("start")), _number(gene.get("end"))
        valid &= bool(_value(gene.get("contig")) and a is not None and b is not None
                      and a.is_integer() and b.is_integer() and 1 <= a <= b
                      and gene.get("strand", gene.get("direction")) in {"+", "-"})
    refs = []
    matching_source = None
    for source in (data.get("suite") or {}).get("sources", []):
        if source.get("sample_id") != sid or source.get("role") != "gene_annotations" or source.get("selected") is False:
            continue
        if source.get("row_count") != len(genes) or not re.fullmatch(r"[a-f0-9]{64}", str(source.get("sha256", ""))):
            continue
        row_numbers = []
        for gene in genes:
            matches = [r for r in gene.get("source_refs", []) if r.get("path") == source.get("path") and r.get("sha256") == source.get("sha256") and type(r.get("row")) is int]
            if len(matches) != 1:
                break
            row_numbers.append(matches[0]["row"])
        if len(row_numbers) == len(genes) and set(row_numbers) == set(range(2, len(genes) + 2)):
            matching_source = source
            refs = [{k: source[k] for k in ("path", "sha256")}]
            if source.get("sheet"):
                refs[0]["sheet"] = source["sheet"]
            break
    complete = bool(valid and matching_source)
    positions, contigs = {}, defaultdict(list)
    if complete:
        for gene in genes:
            contigs[gene["contig"]].append(gene)
        for local in contigs.values():
            local.sort(key=lambda g: (float(g["start"]), float(g["end"]), g["gene_id"]))
            for i, gene in enumerate(local):
                positions[gene["gene_id"]] = (i, gene)
    return dict(complete=complete, positions=positions, source_refs=refs,
                gene_count=len(genes), reason="complete_pinned_gene_annotation_inventory" if complete else "gene_inventory_or_coordinate_binding_incomplete")


def _alignment(data, record, features, alignments, gene):
    """Measure the selected curated reference, never replace it with a raw hit."""
    sha = record.get("protein_sha256")
    if record.get("protein_binding_status") != "verified" or sha != gene.get("protein_sha256", gene.get("sequence_sha256")):
        return dict(status=None, subject=None, sources=[], reason="protein_binding_unassessed")
    real_features = [r for r in features if _refs([r])]
    choices = []
    for row in real_features:
        n = row.get("native_fields", {})
        subject = _value(n.get("curated_blast_resolved_subject_id", n.get("curated_blast_match")))
        if isinstance(subject,str) and subject:
            choices.append((row,n,subject))
    if len(choices) > 1:
        return dict(status=None, subject=None, sources=_refs(real_features), reason="ambiguous_curated_reference_selection")
    candidates = set()
    for row in real_features:
        f = row.get("native_fields", {})
        for name in ("hit_label", "blast_match", "raw_hit_label"):
            value = _value(f.get(name))
            if isinstance(value, str) and value and not value.startswith(("typeI_context:", "typeII_context:", "typeIII_context:", "typeIV_context:")):
                candidates.add(value)
    if choices:
        # The frozen branch joins cognate names and measurements from this
        # selected alignment. A fresh raw hit is a different hypothesis.
        candidates = {choices[0][2]}
    current = [r for r in alignments if r.get("operation") in {"stored_exact_database_alignment", "fresh_reference_alignment"}
               and (r.get("operation") == "stored_exact_database_alignment" or r.get("command_exit_code") == 0)
               and r.get("protein_sha256") == sha
               and re.fullmatch(r"[a-f0-9]{64}", str(r.get("reference_sequence_sha256", "")))
               and _refs([r]) and _refs([{"source_refs": r.get("reference_metadata_source_refs", [])}])
               and r.get("native_fields", {}).get("subject_id") in candidates]
    if len(current) > 1:
        return dict(status=None, subject=None, sources=_refs(current), reason="ambiguous_bound_reference_alignments")
    if current:
        row = current[0]
        n = row["native_fields"]
        refs = _refs([row, {"source_refs": row["reference_metadata_source_refs"]}])
        subject = n.get("resolved_subject_id", n["subject_id"])
        if not isinstance(subject,str) or not _value(subject):
            return dict(status=None,subject=None,sources=refs,reason="reference_subject_identifier_invalid")
        if row.get("alignment_status") == "no_alignment":
            return dict(status=None, subject=None, measured_subject=subject, sources=refs,
                        reason="subject_limited_reference_no_hit_not_global_absence")
        qlen = _number(row.get("query_sequence_length", n.get("query_length", gene.get("protein_length"))))
        rlen = _number(n.get("reference_length"))
        spans = []
        for coverage,a,b,length in (("query_coverage", "qstart", "qend", qlen), ("reference_coverage", "sstart", "send", rlen)):
            explicit = _number(n.get(coverage))
            if explicit is not None:
                spans.append(explicit)
                continue
            # Enclosing coordinates of multiple HSPs include unaligned gaps.
            # Only a single-HSP receipt can use an interval as its coverage.
            if _number(n.get("dnmb_hsp_count")) not in {None, 1}:
                spans.append(None)
                continue
            first, last = _number(n.get(a)), _number(n.get(b))
            span = abs(last-first)+1 if first is not None and last is not None else _number(n.get("length"))
            spans.append(min(1., span/length) if length and span is not None else None)
        nums = [_number(n.get("identity")), *spans, _number(n.get("evalue"))]
        reference_family, reference_role = _value(n.get("reference_family")), _value(n.get("reference_role"))
        role_match = None if reference_role is None else reference_role == record.get("role")
        family_match = None if reference_family is None else reference_family == "Type II"
        source_kind = row["operation"]
    else:
        if len(choices) != 1:
            return dict(status=None, subject=None, sources=_refs(real_features), reason="coherent_curated_reference_unavailable")
        row,n,subject = choices[0]
        refs = _refs([row])
        nums = [_number(n.get("curated_blast_" + s)) for s in ("identity", "query_coverage", "reference_coverage", "evalue")]
        reference_family, reference_role = _value(n.get("curated_blast_family")), _value(n.get("curated_blast_role"))
        role_match = _and(None if reference_role is None else reference_role == record.get("role"), _bool(n.get("curated_blast_role_compatible")))
        family_match = _and(None if reference_family is None else reference_family == "Type II", _bool(n.get("curated_blast_family_compatible")))
        source_kind = "stored_curated_alignment"
    # Exact original prerequisites: no inferred Type-II family, fused role or
    # generic restriction annotation can substitute for this same alignment.
    thresholds = [.5, .8, .8, 1e-5]
    measured = [None if v is None else (0 <= v <= limit if i==3 else limit <= v <= 1) for i,(v,limit) in enumerate(zip(nums,thresholds))]
    status = _and(role_match, family_match, *measured)
    measured_subject = subject
    if source_kind == "fresh_reference_alignment" and status is not True:
        # An unsuccessful nominated-subject measurement does not reproduce
        # the native best-reference search. It cannot exclude other homologs,
        # or use an unselected name to exclude a potential cognate pair.
        status, subject = None, None
    return dict(status=status, subject=subject, sources=refs, reason="exact_reference_alignment_predicate",
                measured_subject=measured_subject,
                reference_selection="frozen_curated" if choices else "nominated_subject",
                source_kind=source_kind, reference_family=reference_family, reference_role=reference_role,
                measurements=dict(zip(("identity","query_coverage","reference_coverage","evalue"), nums)))


def assess_sequence_partial(record, gene, features, annotations):
    """Evaluate the original minimum-length/annotation rule on exact stored protein."""
    sha = record.get("protein_sha256")
    if record.get("evidence", {}).get("valid_sequence") is not True or record.get("protein_binding_status") != "verified":
        return dict(status=None, reason="exact_intact_protein_unavailable", source_refs=[])
    bound, sequences, text = [], [], []
    annotation_columns = False
    for row in [*annotations, *features]:
        f = row.get("native_fields", {})
        sequence = _value(f.get("translation"))
        if not sequence or not _refs([row]):
            continue
        sequence = "".join(str(sequence).split()).upper()
        if hashlib.sha256(sequence.encode()).hexdigest() != sha:
            continue
        sequences.append(sequence)
        bound.append(row)
        for name in ("product", "note", "support", "REBASEfinder_support"):
            if name in f:
                annotation_columns = True
                if _value(f[name]):
                    text.append(str(f[name]))
    if not sequences or len(set(sequences)) != 1 or sha != gene.get("protein_sha256", gene.get("sequence_sha256")):
        return dict(status=None, reason="exact_stored_translation_unavailable", source_refs=[])
    # Native add_sequence_partial_status copies hit.support onto the gene row.
    for row in features:
        if _refs([row]) and _value(row.get("native_fields", {}).get("support")):
            text.append(str(row["native_fields"]["support"]))
            bound.append(row)
    family,role = record.get("family"),record.get("role")
    minimum = {("Type I","R"):700,("Type I","M"):350,("Type I","S"):300,
               ("Type III","R"):650,("Type III","M"):300,("Type IV","R"):100}.get((family,role),
               {"RM":250,"M":200,"R":150,"S":250}.get(role,100))
    # Bind the unchanged source translation above, then apply the frozen native
    # amino-acid normalization for feature measurements (e.g. terminal '*').
    normalized_sequence = re.sub(r"[^A-Z]", "", sequences[0])
    length = len(normalized_sequence)
    annotated = bool(re.search(r"partial|truncat|fragment|incomplete|missing[ -]?(n|c)[ -]?terminus|frameshift|pseudogene", " ".join(text), re.I))
    status = "partial_or_short" if length < minimum or annotated else "full_length" if annotation_columns else None
    return dict(status=status, length_aa=length, expected_min_aa=minimum, annotated_partial=annotated,
                normalized_sequence_sha256=hashlib.sha256(normalized_sequence.encode()).hexdigest(),
                sequence_normalization="native_uppercase_A_to_Z_only_after_raw_hash_binding",
                reason="frozen_native_sequence_partial_table" if status is not None else "partial_annotation_scope_unavailable",
                source_refs=_refs(bound), operation="fresh_stored_sequence_feature_evaluation")


def _enrich_partial(data, records, genes, features, annotations):
    for key,record in records.items():
        if key not in genes:
            continue
        source_states = {_value(r.get("native_fields", {}).get("partial_status")) for r in features[key]}
        source_states.discard(None)
        native_hypotheses = {(_value(r.get("native_fields", {}).get("family_id", r.get("native_fields", {}).get("rm_type"))),
                              _value(r.get("native_fields", {}).get("enzyme_role", r.get("native_fields", {}).get("subunit")))) for r in features[key]}
        hypothesis_matches = native_hypotheses == {(record.get("family"),record.get("role"))}
        if hypothesis_matches and len(source_states) == 1 and next(iter(source_states)) in {"full_length","partial_or_short"}:
            status = next(iter(source_states))
            explicit = any(re.search(r"annotation suggests|pseudogene|frameshift|internal stop|contig edge", str(r.get("native_fields", {}).get("partial_reason", "")), re.I) for r in features[key])
            assessment = dict(status=status,annotated_partial=explicit,source_refs=_refs(features[key]),reason="preserved_native_partial_status")
        elif hypothesis_matches and len(source_states) > 1:
            assessment = dict(status=None,source_refs=_refs(features[key]),reason="contradictory_native_partial_status")
        else:
            assessment = assess_sequence_partial(record,genes[key],features[key],annotations[key])
        assessment["native_hypothesis_compatible"] = hypothesis_matches
        record["sequence_partial_assessment"] = assessment
        e = record["evidence"]
        if e.get("partial") is not None or assessment["status"] is None:
            continue
        rawpartial = assessment["status"] == "partial_or_short"
        qcov,rcov = e.get("query_coverage"),e.get("reference_coverage")
        rescue = _and(e.get("old_strong"),None if qcov is None else qcov>=.9,None if rcov is None else rcov>=.9)
        explicit = assessment.get("annotated_partial") is True
        penalty = False if not rawpartial else True if explicit or rescue is False else None if rescue is None else False
        e["partial"] = penalty
        unknown = set(record.get("unknown_axes",[]))
        if penalty is None:unknown.add("partial")
        else:unknown.discard("partial")
        record["unknown_axes"] = sorted(unknown)
        refs = assessment["source_refs"] + e_source_refs(record,("old_strong","query_coverage","reference_coverage"))
        refs = _refs([{"source_refs":refs}])
        record.setdefault("axis_sources",{})["partial"] = dict(rule="frozen_native_partial_status_and_reference_rescue",
            status="unknown" if penalty is None else "evaluated",source_refs=refs,literal_status=assessment["status"],
            annotated_partial=explicit,reference_rescue=rescue)
        record["source_refs"] = _refs([{"source_refs":[*record.get("source_refs",[]),*refs]}])


def e_source_refs(record, names):
    return [ref for name in names for ref in record.get("axis_sources",{}).get(name,{}).get("source_refs",[])]


def _eligibility(data, record, features, alignments, gene):
    e = record.get("evidence", {})
    family, role = record.get("family", e.get("family")), record.get("role", e.get("role"))
    if family in RM_TYPES and family != "Type II" or role in {"S", "RM"}:
        return dict(status=False, subject=None, sources=[], reason="outside_type_ii_separate_R_M_branch")
    if family != "Type II" or role not in {"R", "M"}:
        return dict(status=None, subject=None, sources=[], reason="candidate_family_or_role_unresolved")
    alignment = _alignment(data, record, features, alignments, gene)
    partial_states = {_value(r.get("native_fields", {}).get("partial_status")) for r in features}
    partial_states.discard(None)
    literal = record.get("sequence_partial_assessment", {}).get("status")
    full_length = (literal == "full_length" if literal is not None else None) if "sequence_partial_assessment" in record else next(iter(partial_states)) == "full_length" if len(partial_states) == 1 else None
    # Unlike the partial penalty, this branch requires literal full_length;
    # reference-rescued partial annotations do not meet the frozen criterion.
    conflict = e.get("conflict")
    alternative = e.get("alternative_annotation")
    conditions = [alignment["status"], full_length,
                  None if conflict is None else not conflict,
                  None if alternative is None else not alternative,
                  e.get("valid_sequence")]
    alignment.update(status=_and(*conditions), full_length_status=full_length,
                     conflict=conflict, alternative_annotation=alternative,
                     sources=alignment["sources"] + _refs(features))
    return alignment


def enrich_neighborhood_evidence(data: dict, compiled: dict) -> dict:
    """Return copied compiled evidence plus exact selected-branch pair provenance."""
    result = deepcopy(compiled)
    genes = {(g["sample_id"],g["gene_id"]): g for g in data.get("genes", [])}
    tables = {}
    for name in ("features", "reference_alignments", "annotations"):
        tables[name] = defaultdict(list)
        for row in (data.get("suite") or {}).get(name, []):
            tables[name][(row.get("sample_id"),row.get("gene_id"))].append(row)
    samples = {r["sample_id"] for r in result.get("records", [])}
    inventory = {sid:_inventory(data,sid) for sid in samples}
    scopes = {s["sample_id"]:s for s in compiled.get("candidate_scope", [])}
    records = {(r["sample_id"],r["gene_id"]):r for r in result.get("records", [])}
    _enrich_partial(data,records,genes,tables["features"],tables["annotations"])
    eligibility = {}
    for key, record in records.items():
        eligibility[key] = (_eligibility(data, record, tables["features"][key], tables["reference_alignments"][key], genes[key])
                            if key in genes else dict(status=None, subject=None, sources=[], reason="canonical_gene_missing"))
    pairs, states, pair_keys = [], {}, set()
    for key, record in records.items():
        e = record["evidence"]
        if type(e.get("divergent_context")) is bool:
            states[key] = dict(value=e["divergent_context"], reason="preserved_explicit_native_context", partners=[])
            continue
        own = eligibility[key]
        inv = inventory[key[0]]
        if own["status"] is False:
            states[key] = dict(value=False, reason="measured_frozen_prerequisite_not_met", partners=[])
            continue
        if not inv["complete"] or record.get("family") != "Type II" or record.get("role") not in {"R","M"}:
            states[key] = dict(value=None, reason="unassessed_inventory_or_assignment", partners=[])
            continue
        unknown_pair = False
        partners = []
        for otherkey, other in records.items():
            if otherkey == key or otherkey[0] != key[0]:
                continue
            other_role = other.get("role")
            if other_role in {"R","M","S","RM"} and other_role == record["role"] or other_role in {"S","RM"}:
                continue
            mate = eligibility[otherkey]
            if mate["status"] is False:
                continue
            a, b = inv["positions"].get(key[1]), inv["positions"].get(otherkey[1])
            if a is None or b is None:
                unknown_pair = True
                continue
            left,right = (a,b) if a[0] < b[0] else (b,a)
            intervening = abs(a[0]-b[0])-1
            if not divergent_geometry(left[1], right[1], intervening):
                continue
            rkey,mkey = (key,otherkey) if record["role"] == "R" else (otherkey,key)
            rs,ms = eligibility[rkey].get("subject"), eligibility[mkey].get("subject")
            if rs and ms and not cognate_subjects(rs,ms):
                if any(eligibility[k].get("reference_selection") == "nominated_subject" for k in (rkey,mkey)):
                    # A positive hit to one nominee is not proof that no
                    # different reference could form a cognate pair.
                    unknown_pair = True
                continue
            if own["status"] is not True or mate["status"] is not True or not rs or not ms:
                unknown_pair = True
                continue
            partners.append(otherkey[1])
            if (rkey,mkey) not in pair_keys:
                pair_keys.add((rkey,mkey))
                refs = inv["source_refs"] + eligibility[rkey]["sources"] + eligibility[mkey]["sources"] + _refs([genes[rkey],genes[mkey]])
                pairs.append(dict(sample_id=key[0],r_gene_id=rkey[1],m_gene_id=mkey[1],
                                  r_protein_sha256=records[rkey]["protein_sha256"],m_protein_sha256=records[mkey]["protein_sha256"],
                                  r_subject=rs,m_subject=ms,gap_bp=int(right[1]["start"])-int(left[1]["end"])-1,
                                  intervening_genes=intervening,source_refs=_refs([{"source_refs":refs}]),
                                  branch="divergent_cognate_neighborhood_not_same_operon",activity_inferred=False))
        complete_scope = scopes.get(key[0],{}).get("status") == "complete"
        value = True if partners else None if unknown_pair or not complete_scope else False
        states[key] = dict(value=value,reason="cognate_divergent_pair" if partners else "potential_pair_unassessed" if value is None else "complete_scoped_neighborhood_no_qualifying_pair",partners=sorted(partners))
    for key, state in states.items():
        record = records[key]
        if state["reason"] == "preserved_explicit_native_context":
            continue
        record["evidence"]["divergent_context"] = state["value"]
        unknown = set(record.get("unknown_axes",[]))
        if state["value"] is None:
            unknown.add("divergent_context")
        else:
            unknown.discard("divergent_context")
        record["unknown_axes"] = sorted(unknown)
        refs = inventory[key[0]]["source_refs"] + eligibility[key]["sources"]
        for gid in state["partners"]:
            refs += eligibility[(key[0],gid)]["sources"] + _refs([genes[(key[0],gid)]])
        refs = _refs([{"source_refs":refs}])
        record.setdefault("axis_sources",{})["divergent_context"] = dict(rule=RULE_VERSION,
            status="unknown" if state["value"] is None else "evaluated", reason=state["reason"],
            partners=state["partners"],source_refs=refs)
        record["source_refs"] = _refs([{"source_refs":[*record.get("source_refs",[]),*refs]}])
    result["neighborhood_evaluation"] = dict(rule_version=RULE_VERSION,rule_provenance=RULE_PROVENANCE,max_gap_bp=Policy().max_divergent_gap_bp,
        max_intervening_genes=MAX_INTERVENING_GENES,identity_min=.5,query_coverage_min=.8,reference_coverage_min=.8,evalue_max=1e-5,
        thresholds_changed=False,activity_inferred=False,detector_runs=0,pairs=pairs,
        inventory=[dict(sample_id=sid,**{k:v for k,v in inv.items() if k!="positions"}) for sid,inv in sorted(inventory.items())],
        candidates=[dict(sample_id=k[0],gene_id=k[1],eligibility=v,**states[k]) for k,v in eligibility.items()])
    return result
