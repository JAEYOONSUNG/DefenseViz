"""Lossless intake of DNMB's exported R–M curation, not a new decision policy.

``curation_records`` is a canonical *observation* array. Each record preserves
the exact upstream row in ``native_fields`` and its file SHA256/row/sheet in
``source_ref``. Explicit upstream fields are projected with their native names;
missing fields stay None. No family, activity, role or retention is inferred.

``binding_status == 'bound'`` requires an exact sample/query identity, a known
current protein hash and source-bound gene metadata. ``protein_binding_origin``
distinguishes an importer-attached current protein from a matching digest or
translation actually reported by the curation source. It does not claim that a
source without a sequence measurement independently validated that sequence.
Consumers must require ``usable_for_attribution`` before using a record. This
only certifies intake consistency, not biological correctness or activity.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json
import math
import re

SCHEMA_VERSION = 1
SOURCE_KIND = "DNMB_REBASEfinder_curation"
MISSING = {"", "NA", "NaN", "nan", "null", "None", "-"}
TEXT_FIELDS = (
    "curation_tier", "curation_decision", "curation_reason", "final_family",
    "final_role", "final_component", "rm_association_class", "annotation_class",
    "evidence_axes", "exclusion_reason",
)
BOOL_FIELDS = (
    "curation_keep", "other_defense_confirmed", "other_defense_conflict",
    "evidence_conflict", "role_conflict", "unresolved_role_conflict",
    "family_conflict", "context_resolves_blast_conflict", "non_dna_mtase_conflict",
)
NUMBER_FIELDS = ("curation_score", "independent_evidence_axes")
CURATION_FIELDS = TEXT_FIELDS + BOOL_FIELDS + NUMBER_FIELDS
# A family or descriptive product alone is not an exported curation result.
CURATION_MARKERS = (
    "curation_keep", "curation_tier", "curation_decision", "curation_score",
    "rm_association_class", "other_defense_confirmed",
)


def _present(value):
    return value is not None and not (isinstance(value, str) and value.strip() in MISSING)


def _hash(value, size):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-fA-F]{%d}" % size, value) is not None


def native_curation_record(sample, gene_id, native, source_ref, gene):
    """Return a source-bound observation, or None for a non-curation row.

    Malformed typed values and conflicting identities are retained in a
    quarantined record. The caller records the issues in the suite consistency
    log and may reject them in strict mode. Native decisions remain unchanged.
    """
    if not any(_present(native.get(key)) for key in CURATION_MARKERS):
        return None
    issues = []

    def issue(code, field, observed, expected):
        issues.append(dict(code=code, field=field, observed=observed, expected=expected))

    record = dict(schema_version=SCHEMA_VERSION, source_kind=SOURCE_KIND,
                  sample_id=sample, gene_id=gene_id, source_ref=deepcopy(source_ref),
                  native_fields=deepcopy(native), decision_namespace="DNMB_native_curation",
                  downstream_policy_reevaluated=False)
    identity = [sample, gene_id, source_ref.get("sha256"), source_ref.get("sheet"), source_ref.get("row")]
    record["record_id"] = "native-curation-" + hashlib.sha256(
        json.dumps(identity, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()[:24]
    for key in TEXT_FIELDS:
        value = native.get(key)
        record[key] = value if _present(value) else None
        if record[key] is not None and not isinstance(value, str):
            issue("invalid_curation_text", key, value, "explicit text or missing")
    for key in BOOL_FIELDS:
        value = native.get(key)
        if not _present(value):
            record[key] = None
        elif isinstance(value, bool):
            record[key] = value
        elif str(value).lower() in {"true", "1", "false", "0"}:
            record[key] = str(value).lower() in {"true", "1"}
        else:
            record[key] = None
            issue("invalid_curation_boolean", key, value, "explicit true/false or missing")
    for key in NUMBER_FIELDS:
        value = native.get(key)
        record[key] = None
        if _present(value):
            try:
                number = float(value)
                if isinstance(value, bool) or not math.isfinite(number):
                    raise ValueError
                if key == "independent_evidence_axes" and (number < 0 or not number.is_integer()):
                    raise ValueError
                record[key] = int(number) if key == "independent_evidence_axes" else number
            except (TypeError, ValueError, OverflowError):
                issue("invalid_curation_number", key, value, "finite number; evidence axes a nonnegative integer")

    bound = gene is not None
    if gene is not None and (gene.get("sample_id"), gene.get("gene_id")) != (sample, gene_id):
        issue("curation_gene_identity_mismatch", "gene", [gene.get("sample_id"), gene.get("gene_id")], [sample, gene_id])
    for key in ("query", "locus_tag", "gene_id"):
        if _present(native.get(key)) and native[key] != gene_id:
            issue("curation_identity_field_conflict", key, native[key], gene_id)
    if _present(native.get("sample_id")) and native["sample_id"] != sample:
        issue("curation_sample_mismatch", "sample_id", native["sample_id"], sample)
    if not _hash(source_ref.get("sha256"), 64) or not source_ref.get("path") or type(source_ref.get("row")) is not int or source_ref["row"] < 2:
        issue("invalid_curation_source_binding", "source_ref", source_ref, "exact file SHA256, path and data row")

    gene = gene or {}
    digest = gene.get("protein_sha256")
    gene_refs = deepcopy(gene.get("source_refs", []))
    protein_origin = "importer_attached"
    measured_sequence = False
    for key, expected, size in (("protein_sha256", digest, 64), ("sequence_sha256", digest, 64),
                                ("protein_md5", gene.get("protein_md5"), 32), ("sequence_md5", gene.get("protein_md5"), 32)):
        value = native.get(key)
        if not _present(value):
            continue
        if not _hash(value, size):
            issue("invalid_curation_sequence_digest", key, value, f"{size}-character hexadecimal digest")
        elif expected is None:
            bound = False
        elif value.lower() != expected.lower():
            issue("curation_sequence_digest_mismatch", key, value, expected)
        else:
            measured_sequence = True
    if _present(native.get("translation")):
        sequence = "".join(str(native["translation"]).split()).upper()
        if not sequence or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ*" for c in sequence):
            issue("invalid_curation_translation", "translation", native["translation"], "protein sequence")
        elif digest is not None:
            actual = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            if actual != digest:
                issue("curation_sequence_digest_mismatch", "translation_sha256", actual, digest)
            else:
                measured_sequence = True
    for key, canonical in (("contig", "contig"), ("direction", "strand"), ("strand", "strand"),
                           ("start", "start"), ("end", "end"), ("protein_length", "protein_length")):
        value = native.get(key)
        if not _present(value) and not (key in {"strand", "direction"} and value == "-"):
            continue
        if key in {"start", "end", "protein_length"}:
            try:
                number = float(value)
                if isinstance(value, bool) or not math.isfinite(number) or number < 1 or not number.is_integer():
                    raise ValueError
                value = int(number)
            except (TypeError, ValueError, OverflowError):
                issue("invalid_curation_coordinates", key, value, "positive integer")
                continue
        if key in {"strand", "direction"} and value not in {"+", "-"}:
            issue("invalid_curation_coordinates", key, value, "+ or -")
        elif canonical not in gene:
            bound = False
        elif value != gene[canonical]:
            issue("curation_gene_metadata_mismatch", key, value, gene[canonical])
    if not _hash(digest, 64) or not gene_refs:
        bound = False
    if measured_sequence:
        protein_origin = "native_digest_verified"
    record.update(
        protein_sha256=digest,
        identity_binding=dict(method="exact_native_query_to_gene", protein_binding_origin=protein_origin,
                              gene_source_refs=gene_refs),
        binding_status="conflicting" if issues else "bound" if bound else "unassessed",
        usable_for_attribution=bool(bound and not issues),
        validation_issues=issues,
    )
    return record
