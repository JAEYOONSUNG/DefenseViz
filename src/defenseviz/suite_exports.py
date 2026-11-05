"""Readable, exact-key DNMBsuite annotation deliverables; no new inference."""
from __future__ import annotations

from collections import defaultdict
import csv
import math
from pathlib import Path
import re


CDS_COLUMNS = (
    "sample_id", "sample_label", "gene_id", "gene_symbol", "protein_id", "product",
    "annotation_description", "domain_annotations", "contig", "start", "end", "strand",
    "protein_length", "family_ids", "categories", "selected_sources", "native_types",
    "native_roles", "rm_type", "rm_role", "retention_state", "rm_score", "rm_score_min",
    "rm_score_max", "rm_decision_reason", "integrated_CDS", "retained_RM", "R_subunit",
    "effector_CDS", "associated_CDS", "priority_CDS", "priority_support_qualified",
    "policy_R_component_support", "effector_status", "associated_status", "motif_status",
    "review_execution", "effector_rules", "activity_assessment", "protein_sha256",
    "annotation_record_count", "excluded_protein_mismatch_records", "source_ref_count", "source_refs",
)
MISSING = {"", "NA", "N/A", "None", "null", "-"}


def _key(row):
    return row["sample_id"], row["gene_id"]


def _index(rows, label):
    result = {}
    for row in rows:
        key = _key(row)
        if key in result:
            raise ValueError(f"Duplicate {label} gene identity: {key!r}")
        result[key] = row
    return result


def _values(value):
    if isinstance(value, (list, tuple, set)):
        return [part for v in value for part in _values(v)]
    if value is None or isinstance(value, dict):
        return []
    text = str(value).strip()
    return [] if text in MISSING else [text]


def _text(values):
    return "; ".join(dict.fromkeys(part for value in values for part in _values(value)))


def _number(value):
    if value is None or isinstance(value, bool):
        return value
    try:
        n = float(value)
    except (ValueError, TypeError):
        return None
    return int(n) if math.isfinite(n) and n.is_integer() else n if math.isfinite(n) else None


def _refs(rows):
    refs = {}
    for row in rows:
        for ref in [row.get("source_ref"), *row.get("source_refs", [])]:
            if not isinstance(ref, dict) or not ref.get("path"):
                continue
            identity = tuple((k, str(ref[k])) for k in sorted(ref) if ref[k] is not None)
            refs[identity] = ref
    return list(refs.values())


def _ref_text(refs):
    return " | ".join(str(r["path"]) + " [" + "; ".join(
        f"{k}={r[k]}" for k in sorted(r) if k != "path" and r[k] is not None) + "]" for r in refs)


def _flat_fields(rows, names):
    return _text(row.get(name) for row in rows for name in names)


def _make_cds_rows(data, result):
    genes = _index(data.get("genes", []), "canonical")
    saved = _index(data.get("saved_decisions", []), "saved decision")
    ledger = _index(result.get("gene_ledger", []), "integrated ledger")
    priorities = _index(result.get("priority_ledger", []), "priority ledger")
    assessments = _index(result.get("automatic_system_review", {}).get("gene_assessments", []), "automatic assessment")
    samples = {r["sample_id"]: r for r in data.get("samples", result.get("samples", []))}
    calls, annotations, effectors = defaultdict(list), defaultdict(list), defaultdict(list)
    for row in data.get("calls", []):
        calls[_key(row)].append(row)
    for row in [*data.get("annotations", []), *data.get("suite", {}).get("annotations", [])]:
        annotations[_key(row)].append(row)
    for row in data.get("effector_evidence", []):
        effectors[_key(row)].append(row)
    categories = defaultdict(list)
    for row in result.get("inventory_rows", {}).get("categories", []):
        for key in row.get("member_keys", []):
            categories[tuple(key)].append(row.get("label", row["id"]))
    retained = {k for k, row in saved.items() if row.get("state") == "retained"}
    needed = set(ledger) | set(saved) | set(priorities)
    missing = needed - genes.keys()
    if missing:
        raise ValueError(f"Annotation export lacks canonical gene identities: {sorted(missing)[:5]!r}")
    if not retained <= ledger.keys() or not priorities.keys() <= ledger.keys():
        raise ValueError("Retained/priority identities are absent from the integrated CDS ledger")
    for key, row in ledger.items():
        if "retained" in row and bool(row["retained"]) != (key in retained):
            raise ValueError(f"Retained decision and CDS ledger disagree: {key!r}")
        if "priority" in row and bool(row["priority"]) != (key in priorities):
            raise ValueError(f"Priority list and CDS ledger disagree: {key!r}")
    made = {}
    for key in sorted(needed):
        gene, decision, entry = genes[key], saved.get(key, {}), ledger.get(key, {})
        assessment = assessments.get(key, {})
        endpoints = assessment.get("endpoints", {})
        selected = [c for c in calls[key] if c.get("selected") is not False]
        protein = gene.get("protein_sha256", gene.get("sequence_sha256"))
        good_annotations = [r for r in annotations[key] if not (protein and r.get("protein_sha256") and r["protein_sha256"] != protein)]
        fields = [gene, *[r.get("native_fields", r) for r in good_annotations]]
        ev = [*effectors[key], *entry.get("effector_evidence", [])]
        refs = _refs([gene, decision, *selected, *good_annotations, *ev])
        membership = entry.get("source_membership", {})
        row = {
            "sample_id": key[0], "sample_label": samples.get(key[0], {}).get("label", key[0]), "gene_id": key[1],
            "gene_symbol": _flat_fields(fields, ("gene", "gene_symbol", "Preferred_name")),
            "protein_id": _flat_fields(fields, ("protein_id",)),
            "product": _flat_fields(fields, ("product",)),
            "annotation_description": _flat_fields(fields, ("description", "Description", "annotation", "function")),
            "domain_annotations": _text(v for f in fields for k, v in f.items() if re.search(r"interpro|pfam|domain", k, re.I) and not isinstance(v, dict)),
            "contig": gene.get("contig"), "start": _number(gene.get("start")), "end": _number(gene.get("end")),
            "strand": gene.get("strand"), "protein_length": _number(gene.get("protein_length")),
            "family_ids": _text(entry.get("family_ids", [])), "categories": _text(categories[key]),
            "selected_sources": _text([n for n, included in membership.items() if included]),
            "native_types": _text(f"{c.get('caller', 'source')}: {c.get('rm_type')}" for c in selected if c.get("rm_type")),
            "native_roles": _text(f"{c.get('caller', 'source')}: {c.get('role')}" for c in selected if c.get("role")),
            "rm_type": decision.get("rm_type"), "rm_role": decision.get("role"), "retention_state": decision.get("state"),
            "rm_score": _number(decision.get("score_total")), "rm_score_min": _number(decision.get("score_min")),
            "rm_score_max": _number(decision.get("score_max")), "rm_decision_reason": decision.get("decision_reason"),
            "integrated_CDS": key in ledger, "retained_RM": key in retained, "R_subunit": entry.get("restriction_role"),
            "effector_CDS": entry.get("effector"), "associated_CDS": entry.get("associated"), "priority_CDS": key in priorities,
            "priority_support_qualified": entry.get("priority_support_qualified"), "policy_R_component_support": decision.get("integrated_R_component_support"),
            "effector_status": endpoints.get("effector", {}).get("status"), "associated_status": endpoints.get("associated", {}).get("status"),
            "motif_status": endpoints.get("motif", {}).get("status"), "review_execution": assessment.get("review_execution"),
            "effector_rules": _text([*[r.get("rule_id") for r in ev], *endpoints.get("effector", {}).get("rule_ids", [])]),
            "activity_assessment": entry.get("activity", assessment.get("activity")), "protein_sha256": protein,
            "annotation_record_count": len(good_annotations), "excluded_protein_mismatch_records": len(annotations[key])-len(good_annotations),
            "source_ref_count": len(refs), "source_refs": _ref_text(refs),
        }
        made[key] = row
    return made, set(ledger), retained, set(priorities)


def _inventory_rows(rows):
    for row in rows:
        result = {"family_or_category": row.get("label", row.get("id")), "genomes": row.get("G", row.get("genome_count")),
                  "CDSs": row.get("N"), "associated_CDSs": row.get("associated_CDS"), "effector_CDSs": row.get("effector_CDS"),
                  "motif_CDSs": row.get("motif_CDS"), "inventory_status": row.get("counts_status")}
        for name, count in row.get("source_counts", {}).items():
            result[name + " CDSs"] = count
        for name in ("associated", "effector", "motif"):
            result[name + "_functional_status"] = row.get("metric_status", {}).get(name)
            result[name + "_observed_count_status"] = row.get("observed_metric_status", {}).get(name)
        yield result


def _system_rows(result):
    units = result.get("system_review", {}).get("units")
    if units is None:
        units = result.get("groups_selected", [])
    for unit in units:
        members = unit.get("member_keys", [])
        row = {"sample_id": unit.get("sample_id"), "provider": unit.get("provider"),
               "system_id": unit.get("native_system_id", unit.get("group_id")), "unit_id": unit.get("unit_id"),
               "families": _text(unit.get("native_families", [unit.get("family")])),
               "unit_kind": unit.get("unit_kind", unit.get("membership_kind")), "member_CDSs": len(members),
               "selected_member_CDSs": len(unit.get("selected_member_keys", members)),
               "member_gene_ids": _text(k[1] for k in members), "defense_scope": unit.get("defense_scope"),
               "review_execution": unit.get("review_execution", {}).get("status")}
        for name, axis in unit.get("axes", {}).items():
            row[name + "_status"] = axis.get("status")
        yield row


def _type_rows(result):
    for row in result.get("type_tables", []):
        out = {k: row.get(k) for k in ("sample_id", "rm_type", "is_total", "unique_CDS", "pool_CDS", "R", "M", "S", "A",
                                           "counts_status", "context_state", "incomplete_R_context", "orphan_M", "orphan_S")}
        out.update({n + " CDSs": v for n, v in row.get("source_counts", {}).items()})
        out.update({n + "_applicability": v for n, v in row.get("role_states", {}).items()})
        yield out


def _columns(rows, default=()):
    return list(dict.fromkeys([*default, *[k for row in rows for k in row]]))


def _tsv(path, rows, columns):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def export_suite_workbook(data, result, out):
    """Write classified Excel and exact filtered annotation tables under out/tables.

    The integrated ledger defines Classified CDSs; the current saved retained
    decisions define R-M retained; the actual priority ledger defines Priority
    CDSs. Counts and predicates are copied, never inferred from a family name.
    """
    try:
        from openpyxl import Workbook
        from openpyxl.cell import WriteOnlyCell
        from openpyxl.styles import Alignment, Font, PatternFill
        from openpyxl.utils import get_column_letter
    except ImportError as exc:
        raise RuntimeError("DNMBsuite Excel export requires openpyxl. Install it with: pip install 'defenseviz[suite]'") from exc
    made, classified_keys, retained_keys, priority_keys = _make_cds_rows(data, result)
    classified = [made[k] for k in sorted(classified_keys)]
    retained = [made[k] for k in sorted(retained_keys)]
    priorities = [made[k] for k in sorted(priority_keys)]
    evaluation = [made[_key(r)] for r in data.get("saved_decisions", [])]
    overview = [
        {"metric": "Genomes", "value": len(data.get("samples", result.get("samples", []))), "scope": "Input cohort"},
        {"metric": "Classified CDSs", "value": len(classified), "scope": "Observed integrated annotation inventory"},
        {"metric": "R-M retained CDSs", "value": len(retained), "scope": "Current selected-policy retained decisions"},
        {"metric": "Priority CDSs", "value": len(priorities), "scope": "Current R-subunit or qualified-effector union"},
        {"metric": "Effector CDSs", "value": sum(r["effector_CDS"] is True for r in classified), "scope": "Qualified annotations; biological activity is not inferred"},
        {"metric": "R-M evaluated CDSs", "value": len(evaluation), "scope": "Current saved decision inventory"},
        {"metric": "Annotation evidence", "value": "Original calls and current decisions remain distinct", "scope": "Non-R-M annotations are not automatically validated by inclusion"},
        {"metric": "Detailed provenance", "value": "Flat TSV tables and analysis_bundle.json", "scope": "Long source references are preserved outside Excel"},
    ]
    for sample in data.get("samples", result.get("samples", [])):
        sid = sample["sample_id"]
        overview.append({"metric": sample.get("label", sid), "value": sum(k[0]==sid for k in classified_keys), "scope": "Classified CDSs · "+sid})
    tabs = [
        ("Overview", overview, ("metric", "value", "scope")),
        ("Classified CDSs", classified, CDS_COLUMNS[:-1]),
        ("R-M retained", retained, CDS_COLUMNS[:-1]),
        ("Priority CDSs", priorities, CDS_COLUMNS[:-1]),
        ("R-M evaluation", evaluation, CDS_COLUMNS[:-1]),
        ("Systems", list(_system_rows(result)), ("sample_id", "provider", "system_id", "member_CDSs")),
        ("Families", list(_inventory_rows(result.get("inventory_rows", {}).get("families", []))), ("family_or_category", "genomes", "CDSs")),
        ("Categories", list(_inventory_rows(result.get("inventory_rows", {}).get("categories", []))), ("family_or_category", "genomes", "CDSs")),
        ("Type summary", list(_type_rows(result)), ("sample_id", "rm_type", "unique_CDS", "R", "M", "S", "A")),
    ]
    tables = Path(out) / "tables"
    paths = {"workbook": tables / "defense_annotations.xlsx", "classified": tables / "classified_CDSs.tsv",
             "retained": tables / "retained_RM_annotations.tsv", "priority": tables / "prioritized_annotations.tsv"}
    if any(p.exists() for p in paths.values()):
        raise FileExistsError("Annotation exports already exist; use a fresh output directory")
    if any(len(rows) >= 1_048_576 for _, rows, _ in tabs):
        raise ValueError("An annotation worksheet exceeds Excel's row limit; split the cohort before exporting")
    tables.mkdir(parents=True, exist_ok=True)
    wb = Workbook(write_only=True)
    wb.properties.title = "DNMB defense annotations"
    wb.properties.description = "Source-traceable recorded annotations; no biological activity inferred"
    header_font, font = Font(name="Arial", size=10, bold=True, color="FFFFFF"), Font(name="Arial", size=10, color="24393B")
    header_fill = PatternFill("solid", fgColor="277F7B")
    band_fill = PatternFill("solid", fgColor="F1F5F5")
    for name, rows, default in tabs:
        ws = wb.create_sheet(name)
        columns = [c for c in _columns(rows, default) if c != "source_refs"]
        ws.freeze_panes = "D2" if name in {"Classified CDSs", "R-M retained", "Priority CDSs", "R-M evaluation"} else "A2"
        ws.sheet_view.showGridLines = False
        ws.auto_filter.ref = f"A1:{get_column_letter(len(columns))}{len(rows)+1}"
        ws.row_dimensions[1].height = 31
        for i, column in enumerate(columns, 1):
            ws.column_dimensions[get_column_letter(i)].width = 58 if column in {"product", "annotation_description", "domain_annotations", "scope"} else 42 if column in {"sample_label", "metric", "member_gene_ids"} else 25 if column in {"sample_id", "gene_id", "protein_id", "contig", "family_ids", "selected_sources", "native_types", "native_roles"} else 20
        header = []
        for column in columns:
            cell = WriteOnlyCell(ws, value=column.replace("_", " "))
            cell.font, cell.fill, cell.alignment = header_font, header_fill, Alignment(vertical="center", wrap_text=True)
            header.append(cell)
        ws.append(header)
        for index, row in enumerate(rows):
            cells = []
            for column in columns:
                value = row.get(column)
                if isinstance(value, (list, dict, tuple, set)):
                    raise ValueError(f"Non-flat annotation export value in {name}/{column}")
                if isinstance(value, str):
                    value = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", "", value)
                    if len(value)>32767:
                        value=value[:32680]+" … [full value in the TSV annotation tables]"
                cell = WriteOnlyCell(ws, value=value)
                if isinstance(value, str):
                    cell.data_type = "s"  # Annotation text is never an Excel formula.
                cell.font = font
                cell.alignment = Alignment(horizontal="center" if isinstance(value, bool) else "right" if isinstance(value, (int,float)) else "left", vertical="top")
                if index % 2:
                    cell.fill = band_fill
                if isinstance(value, (int,float)) and not isinstance(value, bool):
                    cell.number_format = "#,##0" if isinstance(value,int) else "0.###"
                cells.append(cell)
            ws.append(cells)
    wb.save(paths["workbook"])
    for name, rows in [("classified", classified), ("retained", retained), ("priority", priorities)]:
        _tsv(paths[name], rows, CDS_COLUMNS)
    return paths
