"""Compile source-bound R-M policy inputs from saved DNMB suite observations.

Stored observations and exact protein motif features are evaluated locally;
no detector execution or activity inference occurs here. Missing measurements are
unknown. False denotes an applied evidence predicate, never biological absence.
Annotation expressions and motif architecture constraints preserve the original
DNMB curation rules; the numerical Policy is intentionally not modified.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import fields
from functools import lru_cache
from importlib.resources import files
import hashlib
import json
import math
import re

from .analysis.policy import Evidence, Policy

RULE_VERSION = "saved_native_rm_v2"
RULE_PROVENANCE = {
    "annotation_module": "DNMB/R/rebasefinder_curation.R",
    "annotation_source_sha256": "47c058f75373f5b0d0a504db55de6fc1dc32bfd0a577f2a675c4bd909bbc28c0",
    "motif_module": "DNMB/R/db_module_rebasefinder.R",
    "motif_architecture_source_sha256": "b054397e6951d05c9f1fb38a2c4630538fc481ba19edce62de6839fe54355997",
    "interpretation": "annotation evidence; biological activity not inferred",
}
MOTIF_PATTERNS = {'SAM': '[FYW].G.[GA]', 'Amino-IV': '[DNS]PP[YFW]', 'N5C-PC': 'PCQ', 'N5C-ENV': 'ENV', 'N5C-QxRxR': 'Q.R.R', 'MmeI-X': 'GAHY[TS]', 'MmeI-I': 'F[FYLIV]DPACG[CS]GNF', 'MmeI-PDExK': 'D.{8,14}E.{0,3}K', 'HsdR-PD': 'PD.{5,30}[DE].K', 'P-loop': '(?:[AG].{4}GK[ST]|[GAK]SGK[ST]|AAGK[ST])', 'HsdR-WB': 'DE[AHCFQ][DHQ]', 'HsdR-MIII': '[ST]AT', 'PD-ExK': 'PD.{5,30}[DE].K', 'HNH': 'H.{1,3}N.{5,25}H', 'GIY-YIG': 'GIY.{5,20}YIG', 'PLD-HKD': 'H[A-Z]K.{3,7}D', 'ResIII-WA': '(?:[AG].{4}GK[ST]|[GAK]SGK[ST]|AAGK[ST])', 'ResIII-WB': 'DE.H', 'ResIII-MIII': '[ST]AT', 'ResIII-PD': 'PD.{1,30}[DE].K'}
MOTIF_SOURCE_SHA256 = "595a77a42b9c818d3ad855c54d79acc04f353d7872ec4b06354a62562bf183c4"
RM_TYPES = {"Type I", "Type II", "Type III", "Type IV"}
ROLES = {"R", "M", "S", "RM"}
MISSING = {"", "NA", "N/A", "NaN", "nan", "null", "None", "-"}
ANNOTATION_PATTERNS = {'typeiv_mrr': 'pf04471|mrr[_ -]?cat|\\bmrr family\\b', 'typeiv_resiii_pld': 'pf04851|resiii|dexhc_re_i_iii_res', 'dna_mtase': 'dna[^|;]{0,45}(methylase|methyltransferase)|(?:methylase|methyltransferase)[^|;]{0,45}dna|adenine-specific dna-methyl|cytosine-specific dna-methyl|(?:type[ _-]?(?:i|iii|1|3)\\b|restriction[- ]modification system)[^|;]{0,45}(subunit m|mod subunit|modification methylase)|type[ _-]?ii[gls]?[^|;]{0,55}(restriction enzyme/methyltransferase|restriction-modification enzyme|methyltransferase)|restriction enzyme/methyltransferase|\\bhsdm\\b|hsdm_n|rmtype1_m|rmtype3_mod|n6_n4_mtase|methyltransf[ _-]?d12|pf02086|eco57i|c5_dna_methyl', 'rm_rease': 'restriction[- ]modification[^|;]{0,45}(endonuclease|restriction enzyme|subunit r|r subunit)|restriction (endonuclease|enzyme)|type[ _-]?(i|ii|iii|iv)[^|;]{0,40}(restriction endonuclease|restriction enzyme|res subunit|subunit r|r subunit)|type[ _-]?ii[gls]?[^|;]{0,40}restriction-modification enzyme|\\bmmei\\b|\\bhsdr\\b|hsdr_n|ecor124i|resiii|rmtype1_r|rmtype3_res|modification-dependent restriction|\\bmrr\\b|\\bmcr[abc]\\b', 'rm_specificity': 'specificity subunit|target recognition domain|(?:type[ _-]?i|restriction[- ]modification system)[^|;]{0,45}subunit s|\\bhsds\\b|methylase_s|restriction endonuclease subunit s|rmtype1_s', 'non_dna_mtase': '(?:rrna|trna|ribosomal rna|rna)[^|;]{0,55}(?:methylase|methyltransferase)|(?:methylase|methyltransferase)[^|;]{0,55}(?:rrna|trna|ribosomal rna)|protein(?:(?!dna)[^|;]){0,25}(?:methylase|methyltransferase)|(?:methylase|methyltransferase)(?:[, -]+protein(?:[- ][a-z]+)?)\\b|methylated-dna[^|;]{0,45}protein-cysteine|o-?6-alkylguanine[^|;]{0,45}(transferase|methyltransferase)|\\b(?:ogt|mgmt)\\b|serine hydroxymethyltransferase|aminomethyltransferase|glycine hydroxymethyltransferase|chemotaxis[^|;]{0,35}methyl|ubiquinone[^|;]{0,35}methyl|menaquinone[^|;]{0,35}methyl|methylthiotransferase|spou(?:_methylase)?|ftsj(?:-like)?|\\b(?:rsm|rlm|trm|prm|erm)[a-z0-9]*\\b|\\b(?:ksga|cfr|ruma|rumb|fmu|rimo|hemk)\\b|lysine methyltransferase|glutamine methyltransferase|\\bche[rb]\\b|\\bubi[eg]\\b', 'direct_repair_mtase': 'methylated-dna[^|;]{0,45}protein-cysteine|o-?6-alkylguanine[^|;]{0,45}(transferase|methyltransferase)|\\b(?:ogt|mgmt)\\b', 'other_defense': '\\bbrex\\b|brex__|\\bpglx\\b|\\bpglxi\\b|\\bbrx[a-z0-9]*\\b|\\bdisarm\\b|\\bdnd[a-z0-9]*\\b|phosphorothioat|\\bdrm[a-z0-9]*\\b|crispr|\\bcas[0-9][a-z0-9]*\\b|abortive infection|anti[- ]phage|\\bgao[_ -]?rl\\b|toxin[- ]antitoxin|\\babieii?\\b|type iv ta system|\\bmqs[ar]\\b|\\brele\\b|\\bhig[ab]\\b|\\bmaz[ef]\\b|\\bvap[bc]\\b', 'repair_nuclease': '\\bvsr\\b|very short patch|mismatch repair|endonuclease (iii|iv|v|viii)|ap endonuclease|exonuclease|\\bmut[hlms]\\b|\\buvr[abcd]\\b|\\brecbcd?\\b|\\badd[ab]\\b|dna repair nuclease|homing endonuclease|\\bribonuclease\\b', 'unrelated_helicase': 'atp-dependent rna helicase|rna helicase|\\bdbpa\\b|\\bhrpa\\b|\\bsrmb\\b|dna helicase recq|\\brecq\\b|\\buvrd\\b|\\bpcr[ab]\\b|\\brep helicase|\\bding\\b|\\bpria\\b|\\bseca\\b|\\bmfd\\b|\\brecg\\b|molecular chaperone dnak|fe-s protein assembly chaperone hsca|\\bdnak\\b|\\bhsca\\b|hsp70', 'toxin_atpase': 'pf13304|aaa_21|cog1106|k06926|putative abieii toxin|type iv ta system|pf07693|sll1717|kap p-loop', 'unrelated_metabolic': 'phosphonate[^|;]{0,35}phnh|\\bphnh\\b|serine hydroxymethyltransferase|aminomethyltransferase|glycine hydroxymethyltransferase', 'mobile_nuclease': 'transposase|insertion sequence|\\bins[a-z0-9]+\\b|integrase|recombinase|reverse gyrase|terminase|\\biscb\\b|\\btnpb\\b|homing endonuclease|anti-restriction|\\barda\\b', 'non_rm_restriction_fold': '\\byaeq\\b|pf07152|\\byran\\b|upf0102'}
ANNOTATION_COLUMNS = re.compile(
    r"^(product|gene|gene_synonym|note|Description|Preferred_name|PFAMs|KEGG_ko)$|"
    r"^EggNOG_(Description|Preferred_name|PFAMs|KEGG_ko)$|"
    r"^Signature[ ._]?(accession|description)[ ._](CDD|Pfam|NCBIfam|SUPERFAMILY|InterPro|TIGRFAM|Hamap|SMART|PIRSF|Gene3D|FunFam|PANTHER|PRINTS|ProSiteProfiles|ProSitePatterns)$|"
    r"^DefenseFinder_(system_type|system_subtype|gene_name|profiles_in_system)$|"
    r"^PADLOC_(system|protein_name|target_description|hmm_accession|hmm_name|all_domains|best_hits)$", re.I)
MOTIF_ANNOTATION_COLUMNS = re.compile(
    r"^(product|gene|note)$|description|accession|pfam|interpro|cdd|gene3d|superfamily|ncbi|smart|prosite|eggNOG|PADLOC|DefenseFinder", re.I)


def _pattern_digest(patterns):
    return hashlib.sha256(json.dumps(patterns, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


# Source-file hashes above identify the independently checked historical code;
# these hashes identify the actual installed expressions, not an external file.
MOTIF_CATALOG_SHA256 = _pattern_digest(MOTIF_PATTERNS)
ANNOTATION_PATTERNS_SHA256 = _pattern_digest(ANNOTATION_PATTERNS)
RULE_PROVENANCE.update(motif_source_sha256=MOTIF_SOURCE_SHA256,
    installed_motif_patterns_sha256=MOTIF_CATALOG_SHA256,
    installed_annotation_patterns_sha256=ANNOTATION_PATTERNS_SHA256,
    installed_text_selectors_sha256=_pattern_digest({"curation": ANNOTATION_COLUMNS.pattern, "motif": MOTIF_ANNOTATION_COLUMNS.pattern}))
_EXPECTED_PATTERN_HASHES = {
    "installed_motif_patterns_sha256": "55027a76d25fe182cb0f66f85795642618183813b0abf085aa6e66d5fca4b1f8",
    "installed_annotation_patterns_sha256": "c9ddd4c2c75462956bf1d9267aae8576d2c4e09f0cc7b0bd7550912ee537fa30",
    "installed_text_selectors_sha256": "e05a7f146cdb6ba88d36e9154fe6b73ee81fe5d1003fb822cd7446ea32d51dd4",
}
if any(RULE_PROVENANCE[key] != expected for key, expected in _EXPECTED_PATTERN_HASHES.items()):
    raise RuntimeError("Installed native R-M expressions differ from the audited pattern version")
RULE_PROVENANCE.update(
    divergent_module="policy/source/divergent_context.py",
    divergent_source_sha256="1be7ba45272e71e32f19ea3d944314e17cb8e34fe9e8cc459fd853a8c3420481",
    installed_compiler_sha256=hashlib.sha256(files("defenseviz").joinpath("automatic_rm.py").read_bytes()).hexdigest())


def _value(v):
    return None if v is None or isinstance(v, str) and v.strip() in MISSING else v


def _num(v):
    try:
        n = float(v)
        return n if math.isfinite(n) and n >= 0 else None
    except (ValueError, TypeError):
        return None


def _bool(v):
    if type(v) is bool:
        return v
    if str(v).lower() in {"true", "1", "yes"}:
        return True
    if str(v).lower() in {"false", "0", "no"}:
        return False
    return None


def _and(*values):
    return False if any(v is False for v in values) else None if any(v is None for v in values) else True


def _or(*values):
    return True if any(v is True for v in values) else None if any(v is None for v in values) else False


def _not(value):
    return None if value is None else not value


def _refs(rows):
    refs = []
    for row in rows:
        for ref in [row.get("source_ref"), *row.get("source_refs", [])]:
            if isinstance(ref, dict) and ref.get("path") and re.fullmatch(r"[a-f0-9]{64}", str(ref.get("sha256", ""))):
                if ref not in refs:
                    refs.append(dict(ref))
    return refs


def _flags(text):
    text = text.lower()
    f = {k: bool(re.search(pattern, text)) for k, pattern in ANNOTATION_PATTERNS.items()}
    gmr = bool(re.search(r"\bgmrsd\b", text)) or (bool(re.search(r"pf03235|duf262|gmrsd[_ -]?n", text)) and bool(re.search(r"pf07510|gmrsd[_ -]?c|\bgmrsd\b", text)))
    mrr = bool(re.search(r"pf04471|mrr[_ -]?cat|\bmrr family\b", text))
    respld = bool(re.search(r"pf04851|resiii|dexhc_re_i_iii_res", text)) and bool(re.search(r"pf13091|pldc_2|pld[- ]like|phospholipase d/nuclease", text))
    f["typeiv_specific"] = gmr or mrr or respld or bool(re.search(r"pvruts1i|pvurts1i|sauusi|mspji|modification-dependent restriction endonuclease", text))
    f["rm_rease"] |= f["typeiv_specific"]
    f["explicit_non_dna_mtase"] = f["non_dna_mtase"] or f["direct_repair_mtase"]
    return f


@lru_cache(maxsize=1)
def _df_models():
    catalog = json.loads(files("defenseviz").joinpath("data/defensefinder_models.json").read_text())
    return frozenset(catalog.get("indices", {}).get("exact_component_membership", {}))


def _component(provider, model):
    """Native component nomenclature, not system-name-to-effector inference."""
    model = str(model or "")
    if provider == "PADLOC":
        m = re.fullmatch(r"(REase|NEase|MTase|Specificity)_(I|II|III|IV)_\d+", model)
        if m:
            return "Type " + m[2], {"REase": "R", "NEase": "R", "MTase": "M", "Specificity": "S"}[m[1]]
    if provider == "DefenseFinder":
        # Require an exact shipped profile identity before interpreting its
        # documented component naming grammar. New opaque models stay unknown.
        if model not in _df_models():
            return None
        m = re.fullmatch(r"(?:RM|RM_Type_(?:I|II|III))__Type_(I|II|III)_(MTases|REases|REase|S)(?:_FAM_\d+|_\d+|\d+)(?:\.einsi_trimmed)?", model)
        if m:
            return "Type " + m[1], {"MTases": "M", "REases": "R", "REase": "R", "S": "S"}[m[2]]
        # A Type IV profile name alone does not distinguish catalytic from
        # binding/auxiliary components. Assignment needs the native R role or
        # independent restriction annotation.
    return None


def _union(intervals):
    total = 0
    end = 0
    for left, right in sorted(intervals):
        if left < 1 or right < left:
            return None
        total += max(0, right - max(left, end + 1) + 1)
        end = max(end, right)
    return total


def _profiles(rows, gene):
    """One profile keeps its own evalue and both coverages; no mixed maxima."""
    out, grouped = [], defaultdict(list)
    length = _num(gene.get("protein_length"))
    for row in rows:
        f = row.get("native_fields", {})
        ns = row.get("namespace", row.get("caller"))
        if ns == "DefenseFinder":
            model = f.get("gene_name", row.get("native_model"))
            hint = _component(ns, model)
            if not hint:
                continue
            reported = _num(f.get("hit_sequence_length", f.get("hit_seq_len")))
            if length and reported and length != reported:
                continue
            out.append(dict(model=model, family=hint[0], role=hint[1],
                            evalue=_num(f.get("i_eval", f.get("hit_i_eval"))),
                            query_coverage=_num(f.get("hit_seq_cov")),
                            model_coverage=_num(f.get("hit_profile_cov")), source_refs=_refs([row])))
        elif ns == "PADLOC":
            # domtbl supports both hmmscan and hmmsearch orientations.
            a, b = _component(ns, f.get("query_name")), _component(ns, f.get("target_name"))
            if a and not b:
                model, hint, plen, qlen = f["query_name"], a, _num(f.get("query_length")), _num(f.get("target_length"))
            elif b and not a:
                model, hint, plen, qlen = f["target_name"], b, _num(f.get("target_length")), _num(f.get("query_length"))
            elif row.get("native_model") and _component(ns, row["native_model"]):
                hint = _component(ns, row["native_model"])
                out.append(dict(model=row["native_model"], family=hint[0], role=hint[1],
                                evalue=_num(f.get("full.seq.E.value")), query_coverage=_num(f.get("target.coverage")),
                                model_coverage=_num(f.get("hmm.coverage")), source_refs=_refs([row])))
                continue
            else:
                continue
            if not qlen or not plen or length and qlen != length:
                continue
            coords = [_num(f.get(n)) for n in ("ali_from", "ali_to", "hmm_from", "hmm_to")]
            if None in coords:
                continue
            if max(coords[:2]) > qlen or max(coords[2:]) > plen:
                continue
            grouped[(ns, model, hint, qlen, plen)].append((coords, row))
    for (_, model, hint, qlen, plen), domains in grouped.items():
        q = _union([(int(c[0]), int(c[1])) for c, _ in domains])
        m = _union([(int(c[2]), int(c[3])) for c, _ in domains])
        evs = [_num(r["native_fields"].get("evalue")) for _, r in domains]
        out.append(dict(model=model, family=hint[0], role=hint[1], evalue=min((v for v in evs if v is not None), default=None),
                        query_coverage=q / qlen if q is not None else None, model_coverage=m / plen if m is not None else None,
                        source_refs=_refs([r for _, r in domains])))
    return [r for r in out if all(r[k] is None or 0 <= r[k] <= 1 for k in ("query_coverage", "model_coverage"))]


def _qualifies(p):
    return (p["evalue"] is not None and p["evalue"] <= Policy().hmm_evalue_max
            and p["query_coverage"] is not None and p["query_coverage"] >= Policy().hmm_query_coverage_min
            and p["model_coverage"] is not None and p["model_coverage"] >= Policy().hmm_model_coverage_min)


def _ordered(positions, names, span):
    def walk(i, chosen):
        if i == len(names):
            return chosen[-1] - chosen[0] <= span
        return any(walk(i + 1, [*chosen, p]) for p in positions.get(names[i], []) if not chosen or p > chosen[-1])
    return walk(0, [])


def _motifs(rows, sequence, text, scan_complete):
    positions = defaultdict(set)
    accepted = []
    if scan_complete and sequence:
        # Python and native R gregexpr both use nonoverlapping matches. The
        # original architecture evaluator deliberately does not apply plotting
        # position windows from motif export labels.
        for name, pattern in MOTIF_PATTERNS.items():
            positions[name].update(m.start() + 1 for m in re.finditer(pattern, sequence))
    for row in rows:
        f = row.get("native_fields", {})
        start, stop = _num(f.get("start_aa")), _num(f.get("end_aa"))
        match = _value(f.get("match"))
        if not sequence or start is None or stop is None or not match:
            continue
        start, stop = int(start), int(stop)
        if start < 1 or stop > len(sequence) or sequence[start - 1:stop] != str(match).upper():
            continue
        name = str(f.get("motif", ""))
        if not scan_complete:
            positions[name].add(start)
        accepted.append(row)
    norm = {"P-loop": {"P-loop", "ResIII-WA"}, "Walker-B": {"HsdR-WB", "ResIII-WB"},
            "Motif-III": {"HsdR-MIII", "ResIII-MIII"}, "PD-ExK": {"HsdR-PD", "PD-ExK", "ResIII-PD"}}
    for key, aliases in norm.items():
        positions[key] = set().union(*(positions[a] for a in aliases))
    motor = _ordered(positions, ["P-loop", "Walker-B", "Motif-III"], 750)
    mmei = _ordered(positions, ["MmeI-X", "MmeI-I", "Amino-IV"], 650)
    rease = bool(re.search(r"restriction[- ]?(modification )?(endonuclease|enzyme)|restriction subunit|atp-dependent endonuclease|\bhsdr\b|\bres subunit\b|modification-dependent restriction|methylated[- ]dna.*restriction", text, re.I))
    rease &= not bool(re.search(r"crispr|\bcas[0-9a-z]*\b|transposase|integrase|recombinase|rna-guided|homing endonuclease|endonuclease (iii|iv|v|viii)|exonuclease|ribonuclease|dnase i|uvr[abc]|\bmut[ls]\b|\brec[abcdgq]\b|add[ab] helicase", text, re.I))
    nuc = (bool(positions["PD-ExK"]) and rease
           or bool(positions["HNH"]) and bool(re.search(r"\bhnh\b|h-n-h|pf01844|pf13391", text, re.I))
           or bool(positions["GIY-YIG"]) and (rease or bool(re.search(r"giy[-_ ]?yig|pf01541|pf02820", text, re.I)))
           or len(positions["PLD-HKD"]) >= 2 and bool(re.search(r"pf13091|pldc_2|pld[- ]like|phospholipase d/nuclease|hkd", text, re.I))
           or bool(positions["MmeI-PDExK"]) and (mmei or bool(re.search(r"\bmmei\b|pf2046[5-7]|pf20473|ipr04681[689]", text, re.I))))
    mt = bool(sequence and 150 <= len(sequence) <= 2000 and (positions["SAM"] and positions["Amino-IV"] or _ordered(positions, ["N5C-PC", "N5C-ENV"], 250) or mmei))
    raw = sorted(name for name in ("PD-ExK", "HNH", "GIY-YIG", "PLD-HKD", "MmeI-PDExK") if positions[name])
    return dict(legacy_motor=True if motor else False if scan_complete else None,
                legacy_nuclease=True if nuc else False if scan_complete else None,
                legacy_mtase_motif=True if mt else False if scan_complete else None,
                raw_nuclease_classes=raw), _refs(accepted)


def _reference_fields(feature, gene, family, role, alignments):
    """Reconstruct exported alignment spans, never reference length from a hit.

    Native DNMB uses aligned spans when available, otherwise alignment length;
    reference metadata must supply the independent subject length. Optional
    imported alignment rows preserve their own source and metadata receipts.
    """
    f = dict(feature)
    raw_subject = _value(f.get("hit_label", f.get("blast_match")))
    curated_subject = (_value(f.get("curated_blast_resolved_subject_id")) or
                       _value(f.get("curated_blast_subject_id")))
    subject = curated_subject or raw_subject
    if not subject or str(subject).startswith(("typeI_context:", "typeII_context:", "typeIII_context:", "typeIV_context:")):
        return f, []
    matches = [r for r in alignments
               if r.get("native_fields", {}).get("subject_id") == subject
               and r.get("operation") in {"stored_exact_database_alignment", "fresh_reference_alignment"}
               and (r.get("operation") == "stored_exact_database_alignment" or r.get("command_exit_code") == 0)
               and r.get("protein_sha256") == gene.get("protein_sha256", gene.get("sequence_sha256"))
               and re.fullmatch(r"[a-f0-9]{64}", str(r.get("reference_sequence_sha256", "")))
               and _refs([r]) and _refs([{"source_refs": r.get("reference_metadata_source_refs", [])}])]
    sources = []
    if len(matches) == 1:
        row = matches[0]
        n = row["native_fields"]
        f["reference_search_scope_limited"] = row.get("operation") == "fresh_reference_alignment" and row.get("search_is_exhaustive") is not True
        f["reference_measurement_diagnostic"] = {"subject_id": subject, "native_fields": dict(n),
            "alignment_status": row.get("alignment_status"), "search_scope": row.get("search_scope", "saved_nominated_subjects"),
            "search_is_exhaustive": row.get("search_is_exhaustive") is True, "source_refs": _refs([row])}
        if row.get("operation") == "fresh_reference_alignment":
            # A current alignment is a new observation. Never borrow an old
            # coverage, compatibility flag or precomputed quality from it.
            for field in ("blast_identity", "blast_evalue", "blast_length", "blast_qstart", "blast_qend", "blast_sstart", "blast_send", "blast_query_coverage", "blast_reference_coverage", "blast_alignment_quality", "blast_role_compatible", "blast_family_compatible", "factor_strong", "factor_moderate", "rebase_reference_aa_len", "reference_aa_len", "reference_family", "reference_role"):
                f.pop(field, None)
        # Scalar metadata is admissible only with a separate hashed receipt.
        metadata_refs = row.get("reference_metadata_source_refs", [])
        mapping = {"identity": "blast_identity", "evalue": "blast_evalue", "qstart": "blast_qstart", "qend": "blast_qend", "sstart": "blast_sstart", "send": "blast_send", "length": "blast_length", "query_coverage": "blast_query_coverage", "reference_coverage": "blast_reference_coverage"}
        for a,b in mapping.items():
            if _value(n.get(a)) is not None:
                f[b] = n[a]
        if _refs([{"source_refs": metadata_refs}]):
            for a,b in (("reference_length", "rebase_reference_aa_len"), ("reference_family", "reference_family"), ("reference_role", "reference_role")):
                if _value(n.get(a)) is not None:
                    f[b] = n[a]
            sources += _refs([{"source_refs": metadata_refs}])
        sources += _refs([row])
        if row.get("alignment_status") == "no_alignment" and row.get("command_exit_code") == 0:
            for field in ("blast_identity", "blast_evalue", "blast_length", "blast_qstart", "blast_qend", "blast_sstart", "blast_send", "blast_query_coverage", "blast_reference_coverage"):
                f.pop(field, None)
            if not f["reference_search_scope_limited"]:
                f["factor_strong"] = False
                f["factor_moderate"] = False
                f["reference_alignment_absent"] = True
    length = _num(f.get("blast_length"))
    qlen = _num(gene.get("protein_length"))
    rlen = _num(f.get("rebase_reference_aa_len", f.get("reference_aa_len")))
    for dest, startkey, endkey, whole in (("blast_query_coverage", "blast_qstart", "blast_qend", qlen),
                                         ("blast_reference_coverage", "blast_sstart", "blast_send", rlen)):
        start, end = _num(f.get(startkey)), _num(f.get(endkey))
        span = abs(end - start) + 1 if start is not None and end is not None else length
        if _value(f.get(dest)) is None and whole and span is not None:
            f[dest] = min(1., span / whole)
    # Candidate labels are the hypothesis being tested, never independent
    # reference metadata. A curated reference is usable only for this subject.
    rf, rr = _value(f.get("reference_family")), _value(f.get("reference_role"))
    fresh = len(matches) == 1 and matches[0].get("operation") == "fresh_reference_alignment"
    measured_subject = subject if len(matches) == 1 else raw_subject
    if not fresh and _value(f.get("curated_blast_subject_id")) == measured_subject:
        rf = rf if rf is not None else _value(f.get("curated_blast_family"))
        rr = rr if rr is not None else _value(f.get("curated_blast_role"))
    if _value(f.get("blast_family_compatible")) is None and rf in RM_TYPES and family in RM_TYPES:
        f["blast_family_compatible"] = rf == family
    if _value(f.get("blast_role_compatible")) is None and rr in ROLES and role in ROLES:
        f["blast_role_compatible"] = rr == role or role == "RM" and rr in {"R", "M"} or role in {"R", "M"} and rr == "RM"
    if f.get("reference_search_scope_limited"):
        values = [_num(f.get(k)) for k in ("blast_identity", "blast_query_coverage", "blast_reference_coverage", "blast_evalue")]
        # A nominated-subject search can establish positive homology support.
        # It cannot rule out a stronger compatible reference elsewhere.
        positive = (all(v is not None for v in values) and values[0] >= .3 and values[1] >= .7
                    and values[2] >= .7 and values[3] <= 1e-3 and _bool(f.get("blast_role_compatible")) is True
                    and _bool(f.get("blast_family_compatible")) is not False)
        if not positive:
            for field in ("blast_identity", "blast_evalue", "blast_length", "blast_qstart", "blast_qend", "blast_sstart", "blast_send", "blast_query_coverage", "blast_reference_coverage", "blast_alignment_quality", "blast_role_compatible", "blast_family_compatible", "factor_strong", "factor_moderate"):
                f.pop(field, None)
    return f, sources


def compile_rm_evidence(data: dict) -> dict:
    """Return auditable evidence inputs without mutating imported observations."""
    if not isinstance(data, dict) or not isinstance(data.get("genes"), list):
        raise ValueError("Automatic R-M evaluation requires a canonical DNMB suite bundle")
    genes = {(r["sample_id"], r["gene_id"]): r for r in data["genes"]}
    suite = data.get("suite", {})
    indexed = {}
    for table in ("features", "annotations", "profiles", "motifs", "contexts", "structures", "reference_alignments"):
        ix = defaultdict(list)
        for row in suite.get(table, []):
            if row.get("sample_id") and row.get("gene_id"):
                ix[(row["sample_id"], row["gene_id"])].append(row)
        indexed[table] = ix
    calls = defaultdict(list)
    for row in data.get("calls", []):
        calls[(row.get("sample_id"), row.get("gene_id"))].append(row)
    decisions = {(r["sample_id"], r["gene_id"]): r for r in data.get("saved_decisions", [])}
    keys = set(decisions) | set(indexed["features"])
    keys.update(k for k, rows in calls.items() if any(r.get("rm_type") in RM_TYPES and r.get("evidence_kind") == "final_member" for r in rows))
    records = []
    for key in sorted(keys):
        sid, gid = key
        gene = genes.get(key)
        if gene is None:
            raise ValueError(f"Automatic R-M candidate has no canonical gene: {sid}/{gid}")
        sha = gene.get("protein_sha256", gene.get("sequence_sha256"))
        rows = indexed["features"][key]
        annotations = indexed["annotations"][key]
        feature = dict(rows[0].get("native_fields", {})) if rows else {}
        all_rows = [gene, *rows, *annotations, *calls[key]]
        src = _refs(all_rows)
        reasons, axes = [], {}
        evidence = {f.name: None for f in fields(Evidence)}
        for field in ("profiles", "raw_nuclease_classes", "reference_nuclease_anchor_profile_indices"):
            evidence[field] = []
        def set_axis(name, value, sources, rule):
            evidence[name] = value
            axes[name] = {"rule": rule, "source_refs": sources, "status": "unknown" if value is None else "evaluated"}
        # Exact full-protein binding is mandatory before combining observations.
        sequences = []
        for row in [*annotations, *rows]:
            seq = _value(row.get("native_fields", {}).get("translation"))
            if seq:
                seq = "".join(str(seq).split()).upper()
                if hashlib.sha256(seq.encode()).hexdigest() != sha:
                    reasons.append("annotation_protein_hash_mismatch")
                    continue
                sequences.append(seq)
        raw_sequence = sequences[0] if sequences else None
        # Bind the unmodified imported protein first. Native motif coordinates
        # then refer to its separately recorded A-Z-normalized sequence.
        sequence = re.sub(r"[^A-Za-z]", "", raw_sequence).upper() if raw_sequence else None
        valid = bool(sequence) if raw_sequence else None
        set_axis("valid_sequence", valid, _refs([gene, *annotations]), "exact_whole_protein_hash_then_native_nonletter_normalization")
        textvalues, motiftextvalues, productvalues, signaturevalues = [], [], [], []
        bound_annotations = []
        for row in annotations:
            raw_sequence = _value(row.get("native_fields", {}).get("translation"))
            declared_hash = row.get("protein_sha256", row.get("sequence_sha256"))
            if declared_hash is not None and declared_hash != sha:
                reasons.append("annotation_protein_hash_mismatch")
                continue
            if raw_sequence and hashlib.sha256("".join(str(raw_sequence).split()).upper().encode()).hexdigest() != sha:
                continue
            bound_annotations.append(row)
        annotations = bound_annotations
        for row in [*annotations, *rows]:
            for name, value in row.get("native_fields", {}).items():
                if MOTIF_ANNOTATION_COLUMNS.search(name) and _value(value) is not None:
                    motiftextvalues.append(str(value))
                # Native profile_support reads every Signature-prefixed
                # field, independently of the narrower curation selector.
                if re.match(r"^Signature[ ._]", name) and _value(value) is not None:
                    signaturevalues.append(str(value))
                if ANNOTATION_COLUMNS.search(name) and _value(value) is not None:
                    textvalues.append(str(value))
                    if name == "product":
                        productvalues.append(str(value))
        text, product, signature = " | ".join(textvalues), " | ".join(productvalues), " | ".join(signaturevalues)
        motiftext = " | ".join(motiftextvalues)
        flags, pflags, sflags = _flags(text), _flags(product), _flags(signature)
        profiles = _profiles([*indexed["profiles"][key], *calls[key]], gene)
        candidates = []
        saved = decisions.get(key, {})
        family = _value(feature.get("family_id", feature.get("rm_type")))
        role = _value(feature.get("enzyme_role", feature.get("subunit")))
        if family is None:
            family = _value(saved.get("family", saved.get("rm_type")))
        if role is None:
            role = _value(saved.get("role"))
        original_family, original_role = family, role
        if role is None:
            role = ("RM" if flags["dna_mtase"] and flags["rm_rease"] else
                    "S" if flags["rm_specificity"] else "M" if flags["dna_mtase"] else
                    "R" if flags["rm_rease"] else None)
        if family in RM_TYPES and role in ROLES:
            candidates.append((family, role, "saved_candidate_assignment" if original_role is not None else "native_missing_role_annotation"))
        final_hints = sorted({_component(r.get("caller"), r.get("native_model")) for r in calls[key] if r.get("evidence_kind") == "final_member"} - {None})
        comparator_hints = sorted({(*hint, "exact_comparator_component") for hint in final_hints}
                                  | {(p["family"], p["role"], "qualifying_comparator_profile") for p in profiles if _qualifies(p)})
        # User-authorized union nomination is separate from DNMB annotation.
        # Comparator observations cannot overrule a native type or role.
        if original_family is None and original_role is None:
            candidates.extend((a, b, why) for a, b, why in comparator_hints if role is None or b == role)
        pairs = {(a, b) for a, b, _ in candidates}
        if len(pairs) == 1:
            family, role = next(iter(pairs))
            status = "supported"
        elif len(pairs) > 1:
            family, role, status = "Unresolved", "unknown", "unresolved"
            reasons.append("conflicting_family_or_role_assignments")
        else:
            family = family if family is not None else "Unresolved"
            role = role if role is not None else "unknown"
            status = "supported" if family in RM_TYPES and role in ROLES else "unresolved"
        if comparator_hints and any((a, b) != (family, role) for a,b,_ in comparator_hints):
            reasons.append("comparator_assignment_disagreement_diagnostic")
        if any(flags[k] for k in ("explicit_non_dna_mtase", "repair_nuclease", "other_defense", "mobile_nuclease", "unrelated_metabolic")):
            reasons.append("alternative_annotation_diagnostic_not_scope_exclusion")
        evidence.update(family=family, role=role)
        axes["family"] = axes["role"] = {"rule": "native_assignment_precedes_union_nomination_hints", "source_refs": _refs([*rows, *calls[key]]), "status": status,
            "comparator_hints": [{"family": a, "role": b, "rule": why} for a,b,why in comparator_hints]}
        # Only compatible, same-profile coherent measurements enter Profile.
        # Comparator HMMs are useful exact component observations but are not
        # silently substituted for the historically selected DNMB profile
        # library. Preserve diagnostics independently from numerical credit.
        evidence["profiles"] = []
        axes["profiles"] = {"rule": "dnmb_profile_library_not_in_native_export", "source_refs": [], "status": "unknown", "comparator_profile_details": profiles}
        set_axis("hmm_role_ambiguous", None, [], "dnmb_profile_library_not_in_native_export")
        ann = {"R": flags["rm_rease"], "M": flags["dna_mtase"], "S": flags["rm_specificity"], "RM": flags["rm_rease"] and flags["dna_mtase"]}.get(role)
        set_axis("annotation", ann if text else None, _refs([*annotations, *rows]), "native_role_specific_annotation_flags")
        # old_profile is the native Signature-column factor, not an HMM label.
        sigann = {"R": sflags["rm_rease"], "M": sflags["dna_mtase"], "S": sflags["rm_specificity"], "RM": sflags["rm_rease"] and sflags["dna_mtase"]}.get(role)
        set_axis("old_profile", sigann if any(any(re.match(r"^Signature[ ._]", n) for n in row.get("native_fields", {})) for row in annotations + rows) else None, _refs(annotations), "native_signature_role_factor")
        motifrows = indexed["motifs"][key]
        # Fresh deterministic evaluation of the original motif/architecture
        # rules on the exact stored protein; no detector or external search.
        motif, mrefs = _motifs(motifrows, sequence, motiftext, valid is True)
        if valid is True:
            mrefs = _refs([gene, *annotations]) + mrefs
        for name, value in motif.items():
            set_axis(name, value, mrefs, "frozen_native_sequence_motif_architecture:" + MOTIF_CATALOG_SHA256)
        if role == "S":
            set_axis("old_motif", False, _refs(rows), "specificity_subunit_motif_not_applicable")
        else:
            mt, nuc, motor = (evidence[n] for n in ("legacy_mtase_motif", "legacy_nuclease", "legacy_motor"))
            m = mt if role == "M" else _and(nuc, motor) if role == "R" else _and(mt, nuc) if role == "RM" else None
            # The selected frozen source passes score_motif_support (lines
            # 970–978), which intentionally differs from the annotation-gated
            # diagnostic motif_support used by its conflict predicates.
            set_axis("old_motif", m, mrefs + _refs(annotations), "selected_native_score_motif_support")
        feature, reference_refs = _reference_fields(feature, gene, family, role, indexed["reference_alignments"][key])
        src += [ref for ref in reference_refs if ref not in src]
        for prefix in ("", "curated_"):
            mapping = {"identity": "blast_identity", "query_coverage": "blast_query_coverage", "reference_coverage": "blast_reference_coverage", "evalue": "blast_evalue", "role_compatible": "blast_role_compatible", "family_compatible": "blast_family_compatible"}
            for dest, native in mapping.items():
                name = prefix + dest
                if name not in evidence:
                    continue
                raw = feature.get(prefix + native)
                v = _bool(raw) if "compatible" in dest else _num(raw)
                if dest.endswith("coverage") or dest == "identity":
                    v = v if v is None or v <= 1 else None
                set_axis(name, v, _refs(rows) + reference_refs, "coherent_exported_reference_measurements")
        # Never infer subject length from the aligned segment. A successful hit
        # or high identity does not establish full reference coverage.
        quality = _value(feature.get("blast_alignment_quality"))
        limited_reference_search = feature.get("reference_search_scope_limited") is True
        if feature.get("reference_alignment_absent") is True:
            set_axis("old_strong", False, _refs(rows) + reference_refs, "completed_exact_reference_alignment_no_hit")
            set_axis("old_moderate", False, _refs(rows) + reference_refs, "completed_exact_reference_alignment_no_hit")
        compatible = evidence["role_compatible"] is True and evidence["family_compatible"] is not False
        if quality in {"strong", "moderate", "weak"} and evidence["role_compatible"] is not None:
            for name in ("strong", "moderate"):
                set_axis("old_" + name, compatible and quality == name, _refs(rows) + reference_refs, "native_alignment_quality_and_compatibility")
        elif all(evidence[n] is not None for n in ("identity", "query_coverage", "reference_coverage", "evalue", "role_compatible", "family_compatible")):
            i, q, r, ev = (evidence[n] for n in ("identity", "query_coverage", "reference_coverage", "evalue"))
            strong = compatible and i >= .5 and q >= .8 and r >= .8 and ev <= 1e-5
            moderate = compatible and not strong and i >= .3 and q >= .7 and r >= .7 and ev <= 1e-3
            set_axis("old_strong", True if strong else None if limited_reference_search else False, _refs(rows) + reference_refs, "frozen_native_reference_strong_thresholds")
            set_axis("old_moderate", True if moderate else None if limited_reference_search else False, _refs(rows) + reference_refs, "frozen_native_reference_moderate_thresholds")
        # A measured failing necessary threshold legitimately rules out the
        # corresponding present-evidence credit even if another axis is absent.
        for name, cuts in (("old_strong", (.5, .8, .8, 1e-5)), ("old_moderate", (.3, .7, .7, 1e-3))):
            values = [evidence[k] for k in ("identity", "query_coverage", "reference_coverage", "evalue")]
            fails = any(v is not None and (v > cutoff if i == 3 else v < cutoff) for i,(v,cutoff) in enumerate(zip(values,cuts)))
            fails = fails or evidence["role_compatible"] is False or evidence["family_compatible"] is False
            if evidence[name] is None and fails and not limited_reference_search:
                set_axis(name, False, _refs(rows) + reference_refs, "measured_reference_threshold_not_met")
        context_values, contextrefs = {}, []
        for ctxfamily, column in (("Type I", "typei_operon_supported"), ("Type II", "typeii_operon_supported"), ("Type III", "typeiii_operon_supported")):
            value = _bool(feature.get(column))
            if value is not None:
                context_values[ctxfamily] = value
                contextrefs += _refs(rows)
        for row in indexed["contexts"][key]:
            if row.get("family") in {"Type I", "Type II", "Type III"}:
                value = _bool(row.get("native_fields", {}).get("operon_supported"))
                if value is not None:
                    ctxfamily = row["family"]
                    # Native evidence is an OR across all exported contexts,
                    # independent of the candidate's current family label.
                    context_values[ctxfamily] = context_values.get(ctxfamily, False) or value
                    contextrefs += _refs([row])
        set_axis("old_context", _or(*(context_values.get(f) for f in ("Type I", "Type II", "Type III"))), contextrefs, "native_any_family_operon_OR_missing_exports_unknown")
        if family in {"Type I", "Type III", "Type IV"}:
            set_axis("divergent_context", False, _refs(rows), "divergent_cognate_extension_is_type_ii_only")
        rawpartial = _value(feature.get("partial_status"))
        same_partial_hypothesis = (original_family == family and original_role == role)
        if rawpartial in {"full_length", "partial_or_short"} and same_partial_hypothesis:
            partial = rawpartial == "partial_or_short"
            explicit = bool(re.search(r"annotation suggests|pseudogene|frameshift|internal stop|contig edge", str(feature.get("partial_reason", "")), re.I))
            qcov, rcov = evidence["query_coverage"], evidence["reference_coverage"]
            rescued = _and(evidence["old_strong"], None if qcov is None else qcov >= .9,
                           None if rcov is None else rcov >= .9)
            set_axis("partial", _and(partial, _not(_and(not explicit, rescued))), _refs(rows), "native_partial_status_with_full_length_reference_rescue")
        elif sequence and any(re.search(r"pseudogene|frameshift|internal stop", v, re.I) for v in textvalues):
            set_axis("partial", True, _refs(annotations), "explicit_partial_annotation")
        elif rawpartial in {"full_length", "partial_or_short"} and not same_partial_hypothesis:
            reasons.append("native_partial_status_requires_reassessment_for_inferred_assignment")
        # Explicit frozen feature factors win over reconstructed values; no
        # default false is used for columns absent from legacy exports.
        aliases = {"legacy_nuclease": "rease_nuclease_verified", "legacy_motor": "rease_motor_verified", "legacy_mtase_motif": "mtase_motif_verified", "old_strong": "factor_strong", "old_moderate": "factor_moderate", "old_profile": "factor_profile", "annotation": "factor_annotation", "old_motif": "factor_motif", "old_context": "factor_context", "strict_structure": "numeric_structure_admissible_structural_evidence", "alternative_annotation": "alternative_function_annotation", "conflict": "evidence_conflict", "divergent_context": "divergent_context"}
        for name, column in aliases.items():
            value = _bool(feature.get(column))
            if value is not None:
                set_axis(name, value, _refs(rows) + reference_refs, "explicit_native_evidence_field:" + column)
        if _bool(feature.get("factor_motif")) is None:
            # Exported legacy motif flags supersede the rescan before the
            # selected native score_motif_support predicate is evaluated.
            mt, nuc, motor = (evidence[n] for n in ("legacy_mtase_motif", "legacy_nuclease", "legacy_motor"))
            motif_score = (False if role == "S" else mt if role == "M" else
                           _and(nuc, motor) if role == "R" else _and(mt, nuc) if role == "RM" else None)
            set_axis("old_motif", motif_score, mrefs + _refs(rows), "selected_native_score_motif_support_after_exported_flags")
        if text:
            # The original rescue and conflict expressions are evaluated with
            # three-valued logic. Missing structure/homology cannot silently
            # make a diagnostic alternative penalty true.
            homology_rescue = _and(evidence["old_strong"], ann)
            context_rescue = _and(evidence["old_context"], flags["rm_rease"])
            diagnostic_motif = _and(evidence["old_motif"], ann)
            axes["diagnostic_motif_support"] = {"rule": "annotation_gated_native_diagnostic_motif_support",
                                                "value": diagnostic_motif, "source_refs": mrefs + _refs(annotations)}
            positive = _or(evidence["old_strong"], evidence["old_moderate"], ann,
                           evidence["old_context"], diagnostic_motif, evidence["strict_structure"])
            defense_values = [str(v) for row in annotations for n,v in row.get("native_fields", {}).items()
                              if n.startswith(("DefenseFinder_", "PADLOC_")) and _value(v) is not None]
            dflags = _flags(" | ".join(defense_values))
            explicit_other = pflags["other_defense"] or dflags["other_defense"]
            other_conflict = _and(flags["other_defense"], not explicit_other, positive)
            other_call = _and(flags["other_defense"], _not(other_conflict))
            mixed_dna = pflags["dna_mtase"] and sflags["dna_mtase"]
            non_dna = pflags["explicit_non_dna_mtase"] or (flags["explicit_non_dna_mtase"] and not mixed_dna)
            non_dna &= not (role == "S" and pflags["rm_specificity"] and flags["rm_specificity"])
            r_candidate = role in {"R", "RM"}
            unrelated = _and(r_candidate, flags["unrelated_helicase"], _not(_or(homology_rescue, context_rescue, evidence["strict_structure"])))
            repair = _and(flags["repair_nuclease"], _or(r_candidate, _and(not flags["dna_mtase"], _not(ann))),
                          _or(pflags["repair_nuclease"], _not(_and(homology_rescue, flags["rm_rease"]))))
            toxin = _and(r_candidate, flags["toxin_atpase"], _not(_and(homology_rescue, flags["rm_rease"])))
            mobile = _and(r_candidate, flags["mobile_nuclease"], _not(_and(homology_rescue, flags["rm_rease"])))
            fold = _and(r_candidate, flags["non_rm_restriction_fold"], _not(_or(homology_rescue, evidence["strict_structure"])))
            alternative = _or(non_dna, repair, unrelated, toxin, mobile, fold, flags["unrelated_metabolic"], other_call)
            if evidence["alternative_annotation"] is None:
                set_axis("alternative_annotation", alternative, _refs(annotations), "native_alternative_function_with_rescue_predicates")
            annotation_conflict = ((role in {"M", "RM"} and flags["rm_rease"] and not flags["dna_mtase"])
                                   or (role == "R" and flags["dna_mtase"] and not flags["rm_rease"])
                                   or (role == "S" and (flags["dna_mtase"] or flags["rm_rease"]) and not flags["rm_specificity"]))
            non_dna_conflict = flags["explicit_non_dna_mtase"] and mixed_dna and not pflags["explicit_non_dna_mtase"]
            reference_conflict = evidence["role_compatible"] is False or evidence["family_compatible"] is False
            conflict = _or(annotation_conflict, non_dna_conflict, other_conflict,
                           _and(reference_conflict, _not(_and(evidence["old_context"], ann))),
                           _and(evidence["alternative_annotation"], positive))
            if evidence["conflict"] is None:
                set_axis("conflict", conflict, _refs([*rows, *annotations]), "native_role_annotation_reference_and_alternative_conflicts")
        # Extension anchors require an explicit registered geometry/profile
        # receipt; this compiler does not manufacture that optional credit.
        set_axis("reference_nuclease_anchor_support", False, [], "optional_anchor_extension_not_run")
        record = dict(sample_id=sid, gene_id=gid, protein_sha256=sha, family=family, role=role,
                      assignment_status=status, assignment_support=[{"family": a, "role": b, "rule": c} for a,b,c in sorted(set(candidates))],
                      evidence=evidence, unknown_axes=sorted({k for k,v in evidence.items() if v is None} | {"profiles"} | ({"raw_nuclease_classes"} if valid is not True else set())),
                      axis_sources=axes, source_refs=src, reason_codes=reasons,
                      reference_measurement_diagnostic=feature.get("reference_measurement_diagnostic"),
                      source_scope=[name for name, yes in (("saved_candidates", key in decisions or bool(rows)), ("original_tool_rm_members", any(r.get("rm_type") in RM_TYPES and r.get("evidence_kind") == "final_member" for r in calls[key]))) if yes],
                      motif_sequence_normalization={"rule": "native_strip_nonletters_then_uppercase", "input_protein_sha256": sha,
                          "normalized_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest() if sequence else None,
                          "normalized_length": len(sequence) if sequence else None,
                          "changed": bool(raw_sequence and sequence != raw_sequence)},
                      biological_activity="not_inferred", protein_binding_status="verified" if sequence and not reasons.count("annotation_protein_hash_mismatch") else "unassessed")
        if not re.fullmatch(r"[a-f0-9]{64}", str(sha or "")) or not src or "annotation_protein_hash_mismatch" in reasons:
            record["assignment_status"] = "unresolved"
            record["reason_codes"].append("missing_exact_protein_or_hashed_source_binding")
        records.append(record)
    scopes = []
    for sample in data.get("samples", []):
        sid = sample["sample_id"]
        rs = [r for r in records if r["sample_id"] == sid]
        unmapped = [r for r in suite.get("unmapped_records", []) if r.get("sample_id") == sid and r.get("module") == "rebasefinder"]
        scopes.append(dict(sample_id=sid, status="partial" if unmapped else "complete", gene_ids=[r["gene_id"] for r in rs],
                           source_refs=_refs([{"source_refs": r["source_refs"]} for r in rs]),
                           unbound_native_candidate_count=len(unmapped), description="All imported saved candidates plus original R-M final members; no genome-wide truth claim"))
    return dict(schema_version=1, rule_version=RULE_VERSION, rules=RULE_PROVENANCE,
                records=records, candidate_scope=scopes,
                sources=[s for s in suite.get("sources", []) if s.get("module") in {"rebasefinder", "defensefinder", "padloc", "interproscan", "suite"}],
                engine_runs=0, motif_scans=sum(r["evidence"]["valid_sequence"] is True for r in records),
                motif_catalog_sha256=MOTIF_CATALOG_SHA256,
                analysis_kind="fresh_evaluation_of_stored_observations_and_exact_protein_features")
