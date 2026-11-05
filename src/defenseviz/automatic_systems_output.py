"""Apply bounded automatic component proposals and exact-set review coverage."""
from __future__ import annotations

from copy import deepcopy

from .automatic_systems import review_systems, summarize_assessments


_AXES = ("effector", "associated", "motif")


def _update_row(row, summary):
    # A scoped status must describe exactly the numeric set already rendered.
    # This also catches accidentally joining one family's review to another.
    for axis in _AXES:
        field = axis + "_CDS"
        if field in row and row[field] != summary[field]:
            raise ValueError(f"Automatic system review disagrees with {field} membership")
    row.setdefault("metric_status", {}).update(summary["metric_status"])
    row["observed_metric_status"] = deepcopy(summary["observed_metric_status"])
    row["observed_count_scope"] = summary["observed_count_scope"]
    row["automatic_assessment"] = {
        "scope": "exact_observed_integrated_members",
        "reviewed_CDS": summary["CDS_count"],
        "unresolved_CDS": deepcopy(summary["unresolved_CDS"]),
        "not_applicable_CDS": deepcopy(summary["not_applicable_CDS"]),
        "exhaustive_functional_discovery": False,
    }


def apply_system_review(data, result):
    """Return ``(data, result)`` after automatic review and scoped aggregation.

    Canonical additions are allowed only for explicit source-bound individual
    effector rules. Native calls, saved R-M decisions and supported annotations
    are never replaced. With no such proposals, numerical memberships remain
    identical. Explicit coverage declarations retain their original precedence;
    the automatic report remains available alongside them for auditing.
    """
    report = review_systems(data, result)
    proposals = report["proposals"]["effector_evidence"]
    if proposals:
        from .analysis.saved import analyze_saved
        data.setdefault("effector_evidence", []).extend(deepcopy(proposals))
        result = analyze_saved(data)
        report = review_systems(data, result)
        if report["proposals"]["effector_evidence"]:
            raise ValueError("Automatic system proposals did not converge after canonical validation")
    report["materialization"] = {
        "new_effector_evidence_records": len(proposals),
        "canonical_reanalysis_runs": int(bool(proposals)),
        "detector_runs": 0,
    }
    options = data.get("options") or {}
    explicit = bool(options.get("assessment_coverage") or options.get("assessment_coverage_by_sample"))
    report["coverage_application"] = {
        "status": "explicit_declarations_preserved" if explicit else "applied_exact_member_scopes",
        "scope": "functional_endpoints_only_inventory_coverage_preserved",
    }
    result["automatic_system_review"] = report
    if explicit:
        return data, result

    assessments = report["gene_assessments"]
    for kind in ("categories", "families"):
        for row in result.get("inventory_rows", {}).get(kind, []):
            if "member_keys" not in row:
                raise ValueError("Inventory row lacks exact automatic-review member scope")
            _update_row(row, summarize_assessments(assessments, row["member_keys"],
                                                  inventory_status=row.get("metric_status", {}).get("inventory", row.get("counts_status", "unassessed"))))
    for row in result.get("cells", []):
        if "member_keys" not in row:
            raise ValueError("Family cell lacks exact automatic-review member scope")
        if any(key[0] != row["sample_id"] for key in row["member_keys"]):
            raise ValueError("Automatic family-cell scope crosses samples")
        _update_row(row, summarize_assessments(assessments, row["member_keys"],
                                              inventory_status=row.get("metric_status", {}).get("inventory", row.get("counts_status", "unassessed"))))
    _update_row(result["global_summary"], report["summary"])
    by_sample = {r["sample_id"]: r for r in report["counts_by_sample"]}
    for row in result.get("sample_summary", []):
        _update_row(row, by_sample[row["sample_id"]])

    availability = result.setdefault("availability", {})
    for axis in _AXES:
        availability[axis] = report["summary"]["metric_status"][axis]
    scoped = availability.setdefault("assessment_coverage_by_sample", {})
    provenance = availability.setdefault("assessment_coverage_provenance_by_sample", {})
    for sid, summary in by_sample.items():
        scoped[sid] = deepcopy(summary["metric_status"])
        provenance[sid] = {
            axis: dict(origin="automatic_exact_member_review", option_path=None,
                       source_refs=[], scope="observed_integrated_members",
                       member_count=summary["CDS_count"], unresolved_count=summary["unresolved_CDS"][axis])
            for axis in _AXES
        }
    # Keep the original declarations visible. This parallel field describes the
    # automatic classification of the observed universe, not user declarations.
    for row in result.get("system_review", {}).get("sample_coverage", []):
        if row["sample_id"] in by_sample:
            row["automatic_endpoint_coverage"] = deepcopy(by_sample[row["sample_id"]]["metric_status"])
    return data, result
