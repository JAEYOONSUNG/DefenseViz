"""Fresh, source-bound alignment of saved REBASE nominees to configured assets.

Historical calls remain unchanged. Current reference sequences and measurements
are never represented as the historical database or as biological activity.
"""
from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
import csv
import hashlib
from importlib.resources import files
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

MISSING = {"", "NA", "N/A", "None", "null", "unknown", "Unresolved", "-"}
ROLES = {"R", "M", "S", "RM"}
TYPES = {"Type I", "Type II", "Type III", "Type IV"}
BLAST_FIELDS = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
                "sstart", "send", "evalue", "bitscore", "qlen", "slen")
HSP_MERGE_SOURCE = "DNMB/R/rebasefinder_metadata.R:.dnmb_rebasefinder_collapse_blast_hsps"
HSP_MERGE_SOURCE_SHA256 = "831de92dd3e2cff2f0a86c7caf88548622f8e565f3ad0ab31ada9a9aa442be93"


def _sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _digest(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def _sequence(value):
    if value is None or str(value).strip() in MISSING:
        return None
    sequence = "".join(str(value).split()).upper()
    return sequence if sequence and re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWYUBZOX]+", sequence) else None


def _reference_sequence(value):
    """Remove an explicit terminal REBASE record separator, never internal text."""
    raw = "".join(str(value).split()).upper() if value is not None else ""
    normalized = raw[:-2] if raw.endswith("<>") else raw
    return _sequence(normalized), (["terminal_REBASE_record_separator_removed"] if normalized != raw else [])


def _value(value):
    return None if value is None or str(value).strip() in MISSING else str(value).strip()


def _interval_union(intervals):
    intervals = sorted((min(a,b), max(a,b)) for a,b in intervals)
    if not intervals:
        return None
    left,right=intervals[0]; total=0
    for a,b in intervals[1:]:
        if a<=right+1:
            right=max(right,b)
        else:
            total+=right-left+1;left,right=a,b
    return total+right-left+1


def _collapse_hsps(hsps):
    """Exact pinned native per-query/per-subject HSP summary.

    Identity weights every HSP's alignment length, including overlaps; coverage
    counts the interval union once. Missing intervals fall back to the highest
    bitscore HSP exactly as the native rule does. Raw HSP rows stay in the TSV.
    """
    if not hsps:
        return None
    best=min(hsps,key=lambda h:(-h['bitscore'] if h.get('bitscore') is not None else math.inf,
                               h['evalue'] if h.get('evalue') is not None else math.inf))
    out=dict(best)
    if len(hsps)==1:
        out['dnmb_aligned_query_length']=(abs(out['qend']-out['qstart'])+1
            if out.get('qstart') is not None and out.get('qend') is not None else out.get('length'))
        out['dnmb_aligned_subject_length']=(abs(out['send']-out['sstart'])+1
            if out.get('sstart') is not None and out.get('send') is not None else out.get('length'))
        out['dnmb_hsp_count']=1
        return out
    if not all(h.get(k) is not None for h in hsps for k in ('qstart','qend','sstart','send')):
        out.update(dnmb_aligned_query_length=out.get('length'),dnmb_aligned_subject_length=out.get('length'),dnmb_hsp_count=1)
        return out
    q_union=_interval_union([(h['qstart'],h['qend']) for h in hsps])
    s_union=_interval_union([(h['sstart'],h['send']) for h in hsps])
    weighted=[(h['pident'], h['length'] if h.get('length') is not None and h['length']>0 else 1)
              for h in hsps if h.get('pident') is not None]
    out['pident']=sum(v*w for v,w in weighted)/sum(w for _,w in weighted) if weighted else None
    scores=[h['bitscore'] for h in hsps if h.get('bitscore') is not None]
    evalues=[h['evalue'] for h in hsps if h.get('evalue') is not None]
    out.update(bitscore=sum(scores) if scores else None,evalue=min(evalues) if evalues else None,
        length=min(q_union,s_union),qstart=min(min(h['qstart'],h['qend']) for h in hsps),
        qend=max(max(h['qstart'],h['qend']) for h in hsps),sstart=min(min(h['sstart'],h['send']) for h in hsps),
        send=max(max(h['sstart'],h['send']) for h in hsps),dnmb_aligned_query_length=q_union,
        dnmb_aligned_subject_length=s_union,dnmb_hsp_count=len(hsps))
    return out


def _family(value):
    value = _value(value)
    if value is None:
        return None
    value = value.replace("_", " ")
    if value in TYPES:
        return value
    match = re.fullmatch(r"Type\s+(II[A-Z]+)", value, re.I)
    return "Type II" if match else None


def _metadata_assignment(row):
    """Parse explicit reference descriptions before historical heuristic fields.

    Older DNMB metadata tests Type II before III and can therefore store III as
    II. Boundary-aware description parsing fixes that projection; the original
    fields and disagreements remain in the receipt.
    """
    text = _value(row.get("enz_type")) or ""
    tokens = re.findall(r"\bType[ _-]*(IV|III|II[A-Z]*|I)(?![A-Z0-9])", text, re.I)
    families = {"Type " + ("IV" if t.upper().startswith("IV") else "III" if t.upper().startswith("III")
                           else "II" if t.upper().startswith("II") else "I") for t in tokens}
    native_family = _family(row.get("reference_family", row.get("rm_type")))
    family = next(iter(families)) if len(families) == 1 else None if families else native_family
    low = text.lower()
    negative_restriction = bool(re.search(r"\b(?:anti[- ]?|non[- ]?)restriction\b|\bnot (?:an? )?restriction\b|\black(?:s|ing)? restriction\b", low))
    fused = bool(re.search(r"restriction enzyme/methyltransferase|restriction[- ]modification enzyme", low))
    restriction_component = bool(re.search(
        r"\brestriction (?:endo)?nuclease\b|\brestriction enzyme\b|\brestriction(?:[- ]modification)?(?: system)? (?:subunit r|r subunit)\b|\b(?:hsdr|res) subunit\b", low))
    inferred_role = (None if negative_restriction else "RM" if fused else
                     "S" if "specificity" in low else "M" if "methyltransferase" in low
                     else "R" if restriction_component else None)
    native_role = _value(row.get("reference_role", row.get("subunit")))
    native_role = native_role if native_role in ROLES else None
    role = None if negative_restriction else inferred_role or native_role
    return dict(reference_family=family, reference_role=role,
        native_family=native_family, native_role=native_role, reference_description=text,
        family_basis="explicit_reference_description" if families else "explicit_reference_field" if family else "unassigned",
        role_basis="explicit_reference_description" if inferred_role else "explicit_reference_field" if role else "unassigned",
        family_metadata_disagreement=bool(native_family and family and native_family != family),
        role_metadata_disagreement=bool(native_role and ((role and native_role != role) or negative_restriction)),
        role_description_conflict=negative_restriction,
        ambiguous_description=len(families) > 1)


def _run(command, commands, *, stdout_path=None):
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    receipt = dict(command=list(map(str, command)), exit_code=result.returncode,
                   elapsed_seconds=time.monotonic() - started, stdout=result.stdout, stderr=result.stderr)
    commands.append(receipt)
    if stdout_path is not None:
        Path(stdout_path).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise ValueError(f"Reference evaluation command failed ({result.returncode}): {command[0]}: {result.stderr.strip()}")
    return result.stdout.strip()


def _load_metadata(cache, requested, work, rscript, commands):
    json_path = cache / "reference_metadata.json"
    rds_path = cache / "rebase_data.rds"
    if json_path.is_file():
        payload = json.loads(json_path.read_text())
        if not isinstance(payload, dict) or payload.get("schema_version") != 1 or not isinstance(payload.get("records"), list):
            raise ValueError("reference_metadata.json requires schema_version 1 and records")
        rows = [dict(row, fasta_id=row.get("fasta_id", row.get("subject_id")),
                     enzyme_name=row.get("enzyme_name", row.get("subject_id"))) for row in payload["records"]]
        return rows, json_path, []
    if not rds_path.is_file():
        raise ValueError("Reference cache requires rebase_data.rds or reference_metadata.json")
    executable = rscript or shutil.which("Rscript")
    if not executable:
        raise ValueError("Rscript is required to read native REBASE RDS metadata")
    requests = work / "requested_subjects.txt"
    requests.write_text("".join(s + "\n" for s in sorted(requested)))
    projection = work / "reference_metadata.tsv"
    reader = files("defenseviz").joinpath("data/reference_metadata.R")
    _run([str(executable), "--version"], commands)
    _run([str(executable), "--vanilla", str(reader), str(rds_path), str(requests), str(projection)], commands)
    with projection.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    return rows, rds_path, [requests, projection]


def _read_fasta(path, wanted):
    """Stream the full asset for database length; retain only requested records."""
    selected = defaultdict(list)
    current, chunks = None, []
    residues = sequences = 0
    def flush():
        if current in wanted:
            selected[current].append("".join(chunks).upper())
    with Path(path).open() as stream:
        for raw in stream:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current, chunks = line[1:].split()[0], []
                sequences += 1
            else:
                if current is None:
                    raise ValueError("Reference FASTA sequence precedes its identifier")
                value = "".join(line.split())
                # REBASE parser exports can retain the terminal record marker
                # <>; punctuation is not an amino-acid residue in BLAST.
                residues += sum(c.isalpha() and c.isascii() for c in value)
                if current in wanted:
                    chunks.append(value)
    flush()
    if not residues or not sequences:
        raise ValueError("Reference FASTA contains no sequences")
    return selected, residues, sequences


def _nominees(data):
    genes = {(r["sample_id"], r["gene_id"]): r for r in data.get("genes", [])}
    annotations = defaultdict(list)
    for row in data.get("suite", {}).get("annotations", []) + data.get("suite", {}).get("features", []):
        annotations[row.get("sample_id"), row.get("gene_id")].append(row)
    matches = defaultdict(list)
    for row in data.get("suite", {}).get("reference_matches", []):
        subject = _value(row.get("reference_id"))
        if subject and not subject.startswith(("typeI_context:", "typeII_context:", "typeIII_context:", "typeIV_context:")):
            if any(c in subject for c in "\n\r\t"):
                raise ValueError("Reference subject identifier contains a control character")
            matches[row.get("sample_id"), row.get("gene_id"), subject].append(row)
    records, omitted = [], []
    for (sid, gid, subject), source in sorted(matches.items()):
        key = sid, gid
        gene = genes.get(key, {})
        digest = gene.get("protein_sha256", gene.get("sequence_sha256"))
        sequences, conflict = set(), False
        for row in annotations[key]:
            value = row.get("native_fields", {}).get("translation")
            declared = row.get("protein_sha256", row.get("sequence_sha256"))
            if declared and declared != digest:
                conflict = True
            if _value(value) is not None:
                raw = "".join(str(value).split()).upper()
                # Canonical DNMB hashes retain terminal stop markers. Verify
                # that exact input before native A-Z normalization for BLAST.
                seq = re.sub(r"[^A-Za-z]", "", raw).upper()
                if not seq or not raw.isascii() or _digest(raw) != digest:
                    conflict = True
                else:
                    sequences.add((raw, seq))
        if len(sequences) != 1 or conflict:
            omitted.append(dict(sample_id=sid, gene_id=gid, subject_id=subject,
                reason="missing_or_conflicting_exact_protein_binding", source_refs=[r.get("source_ref") for r in source]))
            continue
        raw, seq = next(iter(sequences))
        records.append(dict(sample_id=sid, gene_id=gid, subject_id=subject, protein_sha256=digest,
                            sequence=seq, query_source_sequence_sha256=digest,
                            query_sequence_sha256=_digest(seq), query_sequence_length=len(seq),
                            query_sequence_normalizations=["native_nonletter_removal"] if seq != raw else [],
                            source_refs=[r.get("source_ref") for r in source]))
    return records, omitted


def enrich_reference_evidence(data, cache_path, work_dir, blastp="blastp", rscript=None):
    """Return new measurements from a configured current cache, without mutations.

    The optional portable metadata JSON has schema_version=1 and records with
    subject_id, sequence or sequence_sha256, and explicit reference_family/role
    (or native rm_type/subunit/enz_type). Native DNMB RDS is supported directly.
    """
    started = time.monotonic()
    cache, work = Path(cache_path).resolve(), Path(work_dir).resolve()
    if not cache.is_dir():
        raise ValueError("Configured reference cache is not a directory")
    if work == cache or cache in work.parents:
        raise ValueError("Reference outputs must be separate from the read-only cache")
    if work.exists() and any(work.iterdir()):
        raise ValueError("Reference evidence workspace must be empty")
    work.mkdir(parents=True, exist_ok=True)
    fasta = cache / "rebase_db.fasta"
    if not fasta.is_file():
        raise ValueError("Reference cache requires rebase_db.fasta")
    commands, artifacts, records, omitted = [], [], [], []
    resource_hashes = {str(fasta): _sha(fasta)}
    expected_metadata = cache / ("reference_metadata.json" if (cache / "reference_metadata.json").is_file() else "rebase_data.rds")
    if expected_metadata.is_file():
        resource_hashes[str(expected_metadata)] = _sha(expected_metadata)
    receipt_path = work / "reference_evaluation.json"
    try:
        nominees, omitted = _nominees(data)
        wanted = {n["subject_id"] for n in nominees}
        metadata, metadata_path, metadata_artifacts = _load_metadata(cache, wanted, work, rscript, commands)
        artifacts += metadata_artifacts
        if str(metadata_path) not in resource_hashes or _sha(metadata_path) != resource_hashes[str(metadata_path)]:
            raise ValueError("Read-only reference metadata changed while loading")
        by_id = defaultdict(list)
        for row in metadata:
            for ident in {_value(row.get("fasta_id")), _value(row.get("enzyme_name"))} - {None}:
                by_id[ident].append(row)
        selected_ids = {str(row["fasta_id"]) for subject in wanted for row in by_id[subject] if row.get("fasta_id")}
        sequences, db_residues, db_sequences = _read_fasta(fasta, selected_ids)
        resolved, reference_audit = {}, []
        for subject in sorted(wanted):
            rows = by_id[subject]
            reason = None
            if len(rows) != 1:
                reason = "reference_metadata_missing" if not rows else "ambiguous_reference_name"
            else:
                row = rows[0]; candidates = sequences.get(str(row.get("fasta_id")), [])
                raw_sequence = candidates[0] if len(candidates) == 1 else ""
                sequence, normalizations = _reference_sequence(raw_sequence)
                expected_raw = "".join(str(row.get("sequence") or "").split()).upper()
                expected_hash = row.get("sequence_sha256") or (_digest(expected_raw) if expected_raw else None)
                if len(candidates) != 1:
                    reason = "reference_fasta_missing_or_duplicated_identifier"
                elif sequence is None or expected_hash is None or _digest(raw_sequence) != expected_hash:
                    reason = "reference_metadata_fasta_sequence_mismatch"
                elif row.get("seq_length") not in {None, ""} and float(row["seq_length"]) != len(raw_sequence):
                    reason = "reference_metadata_length_mismatch"
                else:
                    assignment = _metadata_assignment(row)
                    resolved[subject] = dict(sequence=sequence, reference_sequence_sha256=_digest(sequence),
                        reference_source_sequence_sha256=_digest(raw_sequence), sequence_normalizations=normalizations,
                        reference_length=len(sequence), fasta_id=row["fasta_id"], metadata_row=row.get("metadata_row"), **assignment)
            reference_audit.append(dict(subject_id=subject, status="resolved" if reason is None else "unresolved", reason=reason,
                metadata_candidate_count=len(rows), **({k:v for k,v in resolved[subject].items() if k!="sequence"} if subject in resolved else {})))
        usable = []
        for n in nominees:
            if n["subject_id"] not in resolved:
                omitted.append({k:v for k,v in n.items() if k!="sequence"} | {"reason":"unresolved_current_reference"})
            else:
                usable.append(n)
        metadata_receipt = work / "reference_mapping.json"
        metadata_receipt.write_text(json.dumps(reference_audit, indent=2) + "\n")
        artifacts.append(metadata_receipt)
        result_file = work / "reference_alignments.tsv"
        query_file, subject_file = work / "queries.faa", work / "reference_subset.faa"
        query_ids, subject_ids = {}, {}
        for n in usable:
            query_ids.setdefault((n["sample_id"],n["gene_id"]), ("q"+str(len(query_ids)), n["sequence"]))
            if n["subject_id"] not in subject_ids:
                subject_ids[n["subject_id"]] = ("r"+str(len(subject_ids)), resolved[n["subject_id"]]["sequence"])
        query_file.write_text("".join(f">{ident}\n{seq}\n" for ident,seq in query_ids.values()))
        subject_file.write_text("".join(f">{ident}\n{seq}\n" for ident,seq in subject_ids.values()))
        artifacts += [query_file,subject_file]
        version = None
        if usable:
            executable = shutil.which(str(blastp)) or str(blastp)
            version = _run([executable,"-version"],commands)
            command = [executable,"-query",str(query_file),"-subject",str(subject_file),"-out",str(result_file),
                "-outfmt","6 "+" ".join(BLAST_FIELDS),"-evalue","0.001",
                "-max_target_seqs",str(max(5,len(subject_ids))),"-dbsize",str(db_residues)]
            _run(command,commands,stdout_path=work/"blast.log")
            artifacts.append(work/"blast.log")
        else:
            result_file.write_text("")
        artifacts.append(result_file)
        aligned = defaultdict(list)
        with result_file.open(newline="") as stream:
            for line,row in enumerate(csv.reader(stream,delimiter="\t"),1):
                if len(row)!=len(BLAST_FIELDS):
                    raise ValueError("Unexpected BLAST output column count")
                values=dict(zip(BLAST_FIELDS,row)); aligned[values["qseqid"],values["sseqid"]].append((line,values))
        result_ref=dict(path=str(result_file),sha256=_sha(result_file))
        metadata_refs=[dict(path=str(metadata_path),sha256=resource_hashes[str(metadata_path)]),
                       dict(path=str(fasta),sha256=resource_hashes[str(fasta)]),
                       dict(path=str(metadata_receipt),sha256=_sha(metadata_receipt))]
        for n in usable:
            ref=resolved[n["subject_id"]]; qid=query_ids[n["sample_id"],n["gene_id"]][0]; rid=subject_ids[n["subject_id"]][0]
            hits=aligned[qid,rid]
            native=dict(subject_id=n["subject_id"],reference_length=ref["reference_length"],
                        reference_family=ref["reference_family"],reference_role=ref["reference_role"])
            source=dict(result_ref)
            measurements=[]
            for line,h in hits:
                nums={k:float(h[k]) for k in BLAST_FIELDS[2:]}
                if not all(math.isfinite(x) and x>=0 for x in nums.values()) or not 0<=nums["pident"]<=100:
                    raise ValueError("BLAST returned invalid numeric values")
                if nums["qlen"]!=len(n["sequence"]) or nums["slen"]!=ref["reference_length"]:
                    raise ValueError("BLAST query/reference length binding mismatch")
                for a,b,length in (("qstart","qend",nums["qlen"]),("sstart","send",nums["slen"])):
                    if not 1<=min(nums[a],nums[b])<=max(nums[a],nums[b])<=length:
                        raise ValueError("BLAST alignment coordinates outside bound protein")
                measurements.append(nums)
            if measurements:
                nums=_collapse_hsps(measurements)
                native.update({k:nums[k] for k in ("evalue","qstart","qend","sstart","send","length","bitscore")})
                native["identity"]=nums["pident"]/100.
                native.update({k:nums[k] for k in ("dnmb_aligned_query_length","dnmb_aligned_subject_length","dnmb_hsp_count")})
                native.update(query_length=len(n["sequence"]),query_coverage=min(1.,nums['dnmb_aligned_query_length']/len(n["sequence"])),
                              reference_coverage=min(1.,nums['dnmb_aligned_subject_length']/ref['reference_length']))
                source["rows"]=[line for line,_ in hits]
                if len(hits)==1:
                    source["row"]=hits[0][0]
            records.append(dict(sample_id=n["sample_id"],gene_id=n["gene_id"],operation="fresh_reference_alignment",
                protein_sha256=n["protein_sha256"],reference_sequence_sha256=ref["reference_sequence_sha256"],
                query_source_sequence_sha256=n["query_source_sequence_sha256"],query_sequence_sha256=n["query_sequence_sha256"],
                query_sequence_length=n["query_sequence_length"],query_sequence_normalizations=n["query_sequence_normalizations"],
                query_sequence_source_ref=dict(path=str(query_file),sha256=_sha(query_file)),
                source_ref=source,reference_metadata_source_refs=metadata_refs,nomination_source_refs=n["source_refs"],
                native_fields=native,alignment_status="aligned" if hits else "no_alignment",command_exit_code=0,
                search_scope="saved_nominated_subjects",search_is_exhaustive=False,
                hsp_merge_rule=HSP_MERGE_SOURCE,hsp_merge_source_sha256=HSP_MERGE_SOURCE_SHA256,
                database_relation="configured_current_reference_not_assumed_historical",query_id=qid,reference_id=rid))
        if any(_sha(path)!=digest for path,digest in resource_hashes.items()):
            raise ValueError("Read-only reference asset changed during evaluation")
        receipt=dict(schema_version=1,operation="fresh_reference_alignment",status="complete",activity_inferred=False,
            detector_runs=0,homology_search_runs=int(bool(usable)),commands=commands,blast_version=version,
            reference_resource_sha256=resource_hashes,reference_sequences=db_sequences,reference_residues=db_residues,
            original_nominee_pairs=len(nominees),aligned_pairs=sum(r["alignment_status"]=="aligned" for r in records),
            no_alignment_pairs=sum(r["alignment_status"]=="no_alignment" for r in records),unresolved=omitted,
            reference_resolution_status="partial" if omitted else "complete",
            code_sha256={"automatic_reference.py":_sha(__file__),"reference_metadata.R":_sha(files("defenseviz").joinpath("data/reference_metadata.R"))},
            metadata_disagreements=[r for r in reference_audit if r.get("family_metadata_disagreement") or r.get("role_metadata_disagreement")],
            methods=dict(search_evalue=.001,hsps_per_pair="all_reported_significant_hsps",hsp_merge_rule=HSP_MERGE_SOURCE,
                         hsp_merge_source_sha256=HSP_MERGE_SOURCE_SHA256,search_scope="saved_nominated_subjects",search_is_exhaustive=False,
                         search_database_length="full configured FASTA residues",
                         native_calls_changed=False,metadata_type_parser="explicit_description_boundary_aware",
                         source="https://www.ncbi.nlm.nih.gov/books/NBK279684/"),
            elapsed_seconds=time.monotonic()-started,artifacts=[str(p) for p in artifacts]+[str(receipt_path)])
        receipt_path.write_text(json.dumps(receipt,indent=2)+"\n")
        result=deepcopy(data)
        result.setdefault("suite",{}).setdefault("reference_alignments",[]).extend(records)
        result["suite"]["automatic_reference"]=receipt
        inputs=result.setdefault("inputs",[])
        for path,digest in resource_hashes.items():
            if path not in inputs:inputs.append(path)
            result.setdefault("input_sha256",{})[path]=digest
            roles=result.setdefault("input_roles",{}).setdefault(path,[])
            role=dict(module="rebasefinder",kind="configured_current_reference_database",native_result=False)
            if role not in roles:roles.append(role)
        return result
    except Exception as exc:
        receipt_path.write_text(json.dumps(dict(schema_version=1,status="failed",error=str(exc),commands=commands,
            reference_resource_sha256=resource_hashes,elapsed_seconds=time.monotonic()-started),indent=2)+"\n")
        raise
