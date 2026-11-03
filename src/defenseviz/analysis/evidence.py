"""Strict public boundary for the unchanged DNMB numerical evidence predicate.

Missing numerical measurements are explicit nulls. Missing boolean assessments
are not converted to False: an incomplete record remains unassessed as a whole.
"""
from dataclasses import MISSING, asdict, fields
from math import isfinite
from .policy import Evidence, Policy, Profile, CALIBRATION_STATUS


def _number(value, name, *, fraction=False):
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        raise ValueError(f"{name} must be a finite number or null")
    if value < 0 or (fraction and value > 1):
        raise ValueError(f"{name} must be {'in [0, 1]' if fraction else 'nonnegative'}")
    return value


def parse_policy(value=None):
    """Return the historically selected preset, optionally with explicit overrides."""
    if value is None:
        p = Policy()
    elif isinstance(value, Policy):
        p = value
    elif isinstance(value, dict):
        unknown = set(value) - {f.name for f in fields(Policy)}
        if unknown:
            raise ValueError(f"Unknown policy parameters: {sorted(unknown)}")
        p = Policy(**value)
    else:
        raise ValueError("policy must be a parameter object or Policy")
    defaults = asdict(Policy())
    for name, v in asdict(p).items():
        template = defaults[name]
        if isinstance(template, bool) and type(v) is not bool:
            raise ValueError(f"policy.{name} must be boolean")
        if isinstance(template, int) and not isinstance(template, bool) and type(v) is not int:
            raise ValueError(f"policy.{name} must be an integer")
        if isinstance(template, float):
            _number(v, f"policy.{name}", fraction=name.endswith(("coverage_min", "identity_min")))
            if v is None:
                raise ValueError(f"policy.{name} cannot be null")
        if isinstance(template, str) and (not isinstance(v, str) or not v.strip()):
            raise ValueError(f"policy.{name} must be a nonempty string")
    p.validate()
    return p


def parse_profile(value, where="profile"):
    if isinstance(value, Profile):
        value = asdict(value)
    if not isinstance(value, dict) or set(value) != {f.name for f in fields(Profile)}:
        raise ValueError(f"{where} needs exactly the five Profile fields")
    for k in ("family_compatible", "role_compatible"):
        if type(value[k]) is not bool:
            raise ValueError(f"{where}.{k} must be boolean")
    for k in ("evalue", "query_coverage", "model_coverage"):
        _number(value[k], f"{where}.{k}", fraction=k != "evalue")
    return Profile(**value)


def parse_evidence(value):
    """Validate a complete Evidence object; no identifier or label is interpreted."""
    if isinstance(value, Evidence):
        value = asdict(value)
    if not isinstance(value, dict):
        raise ValueError("evidence must be an object")
    fs = {f.name: f for f in fields(Evidence)}
    unknown = set(value) - set(fs)
    missing = {k for k, f in fs.items() if f.default is MISSING and f.default_factory is MISSING} - set(value)
    if unknown or missing:
        raise ValueError(f"Evidence fields: missing={sorted(missing)}, unknown={sorted(unknown)}")
    v = dict(value)
    for k in ("family", "role"):
        if not isinstance(v[k], str) or not v[k].strip():
            raise ValueError(f"evidence.{k} must be a nonempty string")
    for k, f in fs.items():
        if k not in v:
            continue
        if f.type is bool and type(v[k]) is not bool:
            raise ValueError(f"evidence.{k} must be boolean (not missing or a string)")
    if v["family_compatible"] is not None and type(v["family_compatible"]) is not bool:
        raise ValueError("evidence.family_compatible must be boolean or null")
    for k in ("identity", "query_coverage", "reference_coverage", "evalue", "curated_identity",
              "curated_query_coverage", "curated_reference_coverage", "curated_evalue"):
        if k in v:
            _number(v[k], f"evidence.{k}", fraction=not k.endswith("evalue"))
    if not isinstance(v["profiles"], (list, tuple)):
        raise ValueError("evidence.profiles must be a list")
    v["profiles"] = tuple(parse_profile(p, f"evidence.profiles[{i}]") for i, p in enumerate(v["profiles"]))
    raw = v["raw_nuclease_classes"]
    if not isinstance(raw, (list, tuple, set, frozenset)) or any(not isinstance(x, str) for x in raw):
        raise ValueError("evidence.raw_nuclease_classes must be a list of strings")
    v["raw_nuclease_classes"] = frozenset(raw)
    indices = v.get("reference_nuclease_anchor_profile_indices", ())
    if not isinstance(indices, (list, tuple)) or any(type(i) is not int or not 0 <= i < len(v["profiles"]) for i in indices):
        raise ValueError("Anchor profile indices must refer to actual profiles")
    if len(indices) != len(set(indices)):
        raise ValueError("Anchor profile indices must not be duplicated")
    if indices and not v.get("reference_nuclease_anchor_support", False):
        raise ValueError("Anchor indices require explicit anchor support")
    if v.get("reference_nuclease_anchor_support", False) and not indices:
        raise ValueError("Anchor support requires linked profile indices")
    v["reference_nuclease_anchor_profile_indices"] = tuple(indices)
    return Evidence(**v)


def evidence_dict(value):
    out = asdict(value)
    out["profiles"] = [asdict(p) for p in value.profiles]
    out["raw_nuclease_classes"] = sorted(value.raw_nuclease_classes)
    out["reference_nuclease_anchor_profile_indices"] = list(value.reference_nuclease_anchor_profile_indices)
    return out


def policy_metadata(p):
    parameters = asdict(p)
    selected = asdict(Policy())
    overrides = {key: {"selected": selected[key], "value": value}
                 for key, value in parameters.items() if value != selected[key]}
    return {"parameters": parameters,
            "status": "custom_parameters_not_independent_validation" if overrides else CALIBRATION_STATUS,
            "preset": "custom_parameter_override" if overrides else "selected_evidence_policy",
            "is_selected_preset": not bool(overrides),
            "parameter_overrides": overrides, "original_calibration_status": CALIBRATION_STATUS,
            "biological_activity_inferred": False,
            "independent_predictive_metrics_available": False,
            "calibration": ("Custom parameter overrides of a historically calibrated base preset; the custom configuration has not undergone that selection or independent validation."
                            if overrides else "Evidence gates were selected in historical manuscript-cohort calibration; weights and retention cutoff were fixed. This is not independent validation."),
            "homology_thresholds_active": p.recompute_homology,
            "numeric_R_thresholds_active": p.numeric_r_enabled}
