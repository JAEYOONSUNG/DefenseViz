"""DNMB evidence policy over stored numerical features; no identities or labels.

The deployed predicate accepts only Evidence, Profile and Policy objects. Sample,
locus, manuscript membership and experimental labels live in the runner, never here.
"""
from dataclasses import dataclass, asdict
from math import isfinite

# Exact scientific predicate provenance; adapter validation is in evidence.py.
ORIGINAL_SOURCE = "revision027/pipeline/scoring_extension/dnmb_policy.py"
ORIGINAL_SOURCE_SHA256 = "2ade4cfac8162a7cb1db78e3b4f41b1327c3bce75a273b6adbe8d326698360fb"
# The published selected preset overrides the historical class default (0.70).
CALIBRATION_STATUS = "historically_calibrated_candidate_not_independent_validation"

MISSING = {"", "NA", "NaN", "nan", "None", "null", "-"}
RAW_NUCLEASE_CLASSES = frozenset({"PD-ExK", "HNH", "GIY-YIG", "MmeI-PDExK"})
# PLD-HKD is deliberately absent: the native rule needs two HKD occurrences plus
# domain evidence, which a raw class-presence string does not establish.


def numeric(value):
    if value is None or str(value) in MISSING or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (ValueError, TypeError):
        return None
    return value if isfinite(value) else None


def truth(value):
    if isinstance(value,bool):return value
    token=str(value).strip().lower()
    if token in {"true","t","1","1.0"}:return True
    if token in {"false","f","0","0.0","","na","nan","none","null","-"}:return False
    raise ValueError("Unsupported boolean vocabulary: "+repr(value))


@dataclass(frozen=True)
class Profile:
    family_compatible: bool
    role_compatible: bool
    evalue: float | None
    query_coverage: float | None
    model_coverage: float | None


@dataclass(frozen=True)
class Evidence:
    family: str
    role: str
    valid_sequence: bool
    identity: float | None
    query_coverage: float | None
    reference_coverage: float | None
    evalue: float | None
    role_compatible: bool
    family_compatible: bool | None
    hmm_role_ambiguous: bool
    profiles: tuple[Profile, ...]
    raw_nuclease_classes: frozenset[str]
    legacy_nuclease: bool
    legacy_motor: bool
    legacy_mtase_motif: bool
    old_strong: bool
    old_moderate: bool
    old_profile: bool
    annotation: bool
    old_motif: bool
    old_context: bool
    strict_structure: bool
    alternative_annotation: bool
    conflict: bool
    partial: bool
    divergent_context: bool = False
    curated_identity: float | None = None
    curated_query_coverage: float | None = None
    curated_reference_coverage: float | None = None
    curated_evalue: float | None = None
    curated_role_compatible: bool = False
    curated_family_compatible: bool = False
    reference_nuclease_anchor_support: bool = False
    reference_nuclease_anchor_profile_indices: tuple[int, ...] = ()


@dataclass(frozen=True)
class Policy:
    version: str = "DNMB_evidence_calibration_candidate_v1"
    recompute_homology: bool = False
    profile_wired: bool = True
    divergent_context_enabled: bool = True
    max_divergent_gap_bp: int = 800
    numeric_r_enabled: bool = True
    family_motif_credit: bool = True
    new_motif_requires_context: bool = True
    min_capped_evidence_axes: int = 1
    identity_min: float = .50
    query_coverage_min: float = .80
    reference_coverage_min: float = .80
    homology_evalue_max: float = 1e-5
    hmm_evalue_max: float = 1e-10
    hmm_query_coverage_min: float = .70
    hmm_model_coverage_min: float = .70
    r_identity_min: float = .60
    r_query_coverage_min: float = .90
    r_reference_coverage_min: float = .90
    r_evalue_max: float = 1e-5
    retain_cutoff: int = 50
    reference_strong: int = 35
    reference_moderate: int = 25
    profile: int = 25
    annotation: int = 10
    motif: int = 20
    context: int = 20
    structure: int = 25
    alternative_annotation_penalty: int = 15
    conflict_penalty: int = 15
    partial_penalty: int = 10

    def validate(self):
        for key, value in asdict(self).items():
            if key.endswith(("coverage_min", "identity_min")):
                if isinstance(value, bool) or not isfinite(value) or not 0 <= value <= 1:
                    raise ValueError("Invalid fraction: " + key)
            if key.endswith("evalue_max") and (isinstance(value, bool) or not isfinite(value) or value < 0):
                raise ValueError("Invalid evalue: " + key)
        if self.retain_cutoff != 50:
            raise ValueError("Retention cutoff is frozen at 50 for this calibration")
        if self.min_capped_evidence_axes not in (1, 2):
            raise ValueError("Only prespecified axis-count guard 1/2 is allowed")
        if self.max_divergent_gap_bp not in (0, 200, 400, 800, 1500):
            raise ValueError("Context must use a published bounded export")
        fixed = {"reference_strong":35,"reference_moderate":25,"profile":25,"annotation":10,
                 "motif":20,"context":20,"structure":25,"alternative_annotation_penalty":15,
                 "conflict_penalty":15,"partial_penalty":10}
        if any(getattr(self,k) != v for k,v in fixed.items()):
            raise ValueError("Numerical weights remain frozen; this search concerns evidence criteria")


def threshold(value, bound, lower=True):
    return value is not None and (value >= bound if lower else value <= bound)


def profile_pass(profile, policy):
    """All three numeric gates must hold for this one profile."""
    return (profile.family_compatible and profile.role_compatible
            and threshold(profile.evalue, policy.hmm_evalue_max, False)
            and threshold(profile.query_coverage, policy.hmm_query_coverage_min)
            and threshold(profile.model_coverage, policy.hmm_model_coverage_min))


def reference_pass(e, identity, query, reference, evalue, allow_missing_evalue=False,
                   allow_unknown_family=False, curated=False):
    prefix = "curated_" if curated else ""
    value = lambda name: getattr(e, prefix+name)
    return (e.valid_sequence and value("role_compatible")
            and (value("family_compatible") is True or (allow_unknown_family and value("family_compatible") is None))
            and threshold(value("identity"), identity) and threshold(value("query_coverage"), query)
            and threshold(value("reference_coverage"), reference)
            and (threshold(value("evalue"), evalue, False) or (allow_missing_evalue and value("evalue") is None)))


def evaluate(e, p):
    """Return annotation evidence and policy points; never infer biological activity."""
    if p.recompute_homology:
        strong = reference_pass(e,p.identity_min,p.query_coverage_min,p.reference_coverage_min,p.homology_evalue_max,True,True)
        moderate = not strong and reference_pass(e,.30,.70,.70,1e-3,True,True)
    else:
        strong, moderate = e.old_strong, e.old_moderate
    qualifying_profiles = tuple(i for i,h in enumerate(e.profiles) if e.valid_sequence and profile_pass(h,p))
    profile = e.old_profile or (p.profile_wired and bool(qualifying_profiles))
    context = e.old_context or (p.divergent_context_enabled and e.divergent_context)
    raw_nuclease = bool(e.raw_nuclease_classes & RAW_NUCLEASE_CLASSES)
    numeric_reference = reference_pass(e,p.r_identity_min,p.r_query_coverage_min,p.r_reference_coverage_min,p.r_evalue_max,curated=True)
    numeric_r = (p.numeric_r_enabled and e.role in {"R","RM"} and e.family in {"Type I","Type II","Type III","Type IV"}
                 and raw_nuclease and numeric_reference and bool(qualifying_profiles) and not e.hmm_role_ambiguous)
    context_r = numeric_r and context
    # A new score credit needs independent locus-context evidence. The raw motif,
    # reference and HMM admission are correlated and are never counted as three axes.
    compatible_guard = numeric_reference and bool(qualifying_profiles) and not e.hmm_role_ambiguous
    family_eligible = ((e.family == "Type II" and e.role == "R")
                       or (e.family in {"Type I","Type III"} and e.role == "R" and e.legacy_motor)
                       or (e.family == "Type II" and e.role == "RM" and e.legacy_mtase_motif))
    # Revision027 source-bound anchor route: independent mapping and exact-domain
    # checks are performed by the provenance adapter. This evidence cannot use
    # an unrelated qualified profile. Numerical policy parameters are unchanged.
    linked_profiles = e.reference_nuclease_anchor_profile_indices
    if not isinstance(e.reference_nuclease_anchor_support, bool):
        raise ValueError("Anchor support must be boolean")
    if any(isinstance(i,bool) or not isinstance(i,int) or not 0 <= i < len(e.profiles) for i in linked_profiles):
        raise ValueError("Anchor profile indices must nominate actual profiles")
    reference_anchor_r = (e.reference_nuclease_anchor_support and family_eligible and compatible_guard
                          and bool(set(linked_profiles) & set(qualifying_profiles)))
    admitted_nuclease = e.legacy_nuclease or numeric_r or reference_anchor_r
    integrated_r = e.role in {"R","RM"} and (e.legacy_nuclease or context_r or (reference_anchor_r and context))
    # Type IV is heterogeneous; retain its existing motif score unchanged instead
    # of inventing a generic Type-IV nuclease/motor criterion.
    new_motif = (p.family_motif_credit and family_eligible and admitted_nuclease and compatible_guard
                 and (context or not p.new_motif_requires_context))
    motif = e.old_motif or new_motif
    sequence_points = max(p.reference_strong*strong,p.reference_moderate*moderate,p.profile*profile,p.annotation*e.annotation)
    motif_points,context_points,structure_points = p.motif*motif,p.context*context,p.structure*e.strict_structure
    alt = p.alternative_annotation_penalty*e.alternative_annotation
    conflict = p.conflict_penalty*(e.conflict and not e.alternative_annotation)
    partial = p.partial_penalty*e.partial
    total = max(0,sequence_points+motif_points+context_points+structure_points-alt-conflict-partial)
    axes = sum((sequence_points>0,motif_points>0,context_points>0,structure_points>0))
    decision = {
        "strong":strong,"moderate":moderate,"profile_support":profile,"annotation_support":e.annotation,
        "qualifying_profile_indices":qualifying_profiles,"raw_nuclease_class_observed":raw_nuclease,
        "numeric_reference_pass":numeric_reference,"numeric_R_component_admission":numeric_r,
        "context_qualified_R_component_support":context_r,"integrated_R_component_support":integrated_r,
        "legacy_nuclease_verified":e.legacy_nuclease,"legacy_motor_verified":e.legacy_motor,
        "family_motif_eligible":family_eligible,"new_motif_credit":new_motif and not e.old_motif,
        "motif_support":motif,"context_support":context,"divergent_cognate_context":p.divergent_context_enabled and e.divergent_context,
        "strict_structure_support":e.strict_structure,"alternative_annotation":e.alternative_annotation,
        "conflict":e.conflict,"partial":e.partial,
        "score_sequence":sequence_points,"score_motif":motif_points,"score_context":context_points,
        "score_structure":structure_points,"score_alternative_penalty":alt,"score_conflict_penalty":conflict,
        "score_partial_penalty":partial,"score_total":total,"capped_evidence_axis_count":axes,
        "keep":total>=p.retain_cutoff and axes>=p.min_capped_evidence_axes,
    }
    if e.reference_nuclease_anchor_support:
        decision.update(reference_nuclease_anchor_support=reference_anchor_r,
                        reference_nuclease_anchor_profile_indices=linked_profiles,
                        reference_nuclease_anchor_motif_credit=bool(reference_anchor_r and new_motif and not e.old_motif))
    return decision


def from_feature(row, profiles=(), divergent=False):
    """Explicit allowlist adapter; identifier/product/label values are discarded."""
    f = lambda key: truth(row.get("factor_"+key))
    family = row.get("blast_family_compatible")
    family = None if family is None or family in MISSING else truth(family)
    return Evidence(
        family=row["family_id"],role=row["enzyme_role"],valid_sequence=(numeric(row.get("sequence_length_aa")) or 0)>0,
        identity=numeric(row.get("blast_identity")),query_coverage=numeric(row.get("blast_query_coverage")),
        reference_coverage=numeric(row.get("blast_reference_coverage")),evalue=numeric(row.get("blast_evalue")),
        role_compatible=truth(row.get("blast_role_compatible")),family_compatible=family,
        hmm_role_ambiguous=truth(row.get("hmm_role_ambiguous")),profiles=tuple(profiles),
        raw_nuclease_classes=frozenset(str(row.get("rease_motif_hits_raw","")).split("+")),
        legacy_nuclease=truth(row.get("rease_nuclease_verified")),legacy_motor=truth(row.get("rease_motor_verified")),
        legacy_mtase_motif=truth(row.get("mtase_motif_verified")),old_strong=f("strong"),old_moderate=f("moderate"),
        old_profile=f("profile"),annotation=f("annotation"),old_motif=f("motif"),old_context=f("context"),
        strict_structure=truth(row.get("numeric_structure_admissible_structural_evidence")),
        alternative_annotation=f("alternative_annotation"),conflict=f("conflict"),partial=f("partial"),
        divergent_context=bool(divergent),curated_identity=numeric(row.get("curated_blast_identity")),
        curated_query_coverage=numeric(row.get("curated_blast_query_coverage")),
        curated_reference_coverage=numeric(row.get("curated_blast_reference_coverage")),
        curated_evalue=numeric(row.get("curated_blast_evalue")),
        curated_role_compatible=truth(row.get("curated_blast_role_compatible")),
        curated_family_compatible=truth(row.get("curated_blast_family_compatible")))


def profile_from_row(row):
    return Profile(truth(row["native_family_compatible"]),truth(row["native_role_compatible"]),
                   numeric(row["min_sequence_evalue"]),numeric(row["query_union_coverage_fraction"]),
                   numeric(row["profile_union_coverage_fraction"]))
