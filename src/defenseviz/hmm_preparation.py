#!/usr/bin/env python3
"""Prepare an exact HMM/role-metadata pair from authorized source files.

Only explicit MicrobeMod metadata supplies role/type hints. This script does not
download resources, run searches, modify GA cutoffs or infer query function.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys
import time

ROLE_MAP = {"MT": ("M",), "RE": ("R",), "SP": ("S",), "IIG": ("M", "R")}
SYSTEM_MAP = {prefix + name: name for prefix in ("RM_Type_", "Type_")
              for name in ("I", "II", "III", "IV", "IIG", "IIM")}
FIELDS = ["profile_name", "roles", "type_hints", "off_target", "profile_length", "library_sha256",
          "native_enzyme_type", "native_system", "metadata_lines", "role_mapping_basis",
          "GA_sequence", "GA_domain", "profile_record_sha256"]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def hmm_profiles(path):
    """Read exact record headers; preserve the complete library byte-for-byte."""
    profiles, record = [], []
    for line in Path(path).read_bytes().splitlines(keepends=True):
        if line.startswith(b"HMMER3/"):
            if record:
                raise ValueError("HMM record lacks // terminator")
            record = [line]
        elif record:
            record.append(line)
            if line.strip() == b"//":
                raw = b"".join(record)
                headers = {}
                for header in raw.decode("utf-8").splitlines()[1:]:
                    if header.startswith("HMM "):
                        break
                    parts = header.split()
                    if parts and parts[0] in {"NAME", "LENG", "GA"}:
                        if parts[0] in headers:
                            raise ValueError("Duplicate HMM header: " + parts[0])
                        headers[parts[0]] = parts[1:]
                name = headers.get("NAME", [])
                length = headers.get("LENG", [])
                ga = headers.get("GA", [])
                if len(name) != 1 or len(length) != 1 or not length[0].isdigit() or int(length[0]) <= 0:
                    raise ValueError("Every HMM needs one NAME and positive integer LENG")
                try:
                    valid_ga = len(ga) == 2 and all(math.isfinite(float(x.rstrip(";"))) for x in ga)
                except ValueError:
                    valid_ga = False
                if not valid_ga:
                    raise ValueError("HMM lacks valid native GA sequence/domain cutoffs: " + name[0])
                profiles.append(dict(profile_name=name[0], profile_length=length[0],
                                     GA_sequence=ga[0].rstrip(";"), GA_domain=ga[1].rstrip(";"),
                                     profile_record_sha256=hashlib.sha256(raw).hexdigest()))
                record = []
        elif line.strip():
            raise ValueError("Content outside a HMMER3 record")
    if record:
        raise ValueError("Final HMM record lacks // terminator")
    names = [row["profile_name"] for row in profiles]
    if not names or len(names) != len(set(names)):
        raise ValueError("HMM profiles must have nonempty unique names")
    return profiles


def source_metadata(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames or []
        if not {"Name", "Enzyme_type", "System"} <= set(fields) or len(fields) != len(set(fields)):
            raise ValueError("Source CSV needs unique Name, Enzyme_type and System columns")
        rows = list(reader)
    by_name = {}
    for line, row in enumerate(rows, 2):
        if None in row or None in row.values() or not row["Name"].strip():
            raise ValueError("Malformed source metadata row: " + str(line))
        by_name.setdefault(row["Name"].strip(), []).append((line, {k: v.strip() for k, v in row.items()}))
    return by_name


def prepare(hmm, metadata, off_target, out, *, argv=None):
    started = time.monotonic()
    source_paths = {key: Path(value).expanduser().resolve() for key, value in
                    (("hmm", hmm), ("metadata", metadata), ("off_target", off_target))}
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError("Output must be a new directory: " + str(out))
    inputs = {key: dict(path=str(path), sha256=digest(path)) for key, path in source_paths.items()}
    profiles = hmm_profiles(source_paths["hmm"])
    by_name = source_metadata(source_paths["metadata"])
    excluded = {line.strip() for line in source_paths["off_target"].read_text(encoding="utf-8").splitlines()
                if line.strip()}
    rows = []
    for profile in profiles:
        name = profile["profile_name"]
        matches = by_name.get(name, [])
        enzymes = sorted({row["Enzyme_type"] for _, row in matches})
        systems = sorted({row["System"] for _, row in matches})
        is_off_target = name in excluded
        if is_off_target:
            roles, types = [], []
            basis = "explicit source off-target list; no R/M/S assignment"
        else:
            if len(enzymes) != 1 or enzymes[0] not in ROLE_MAP:
                raise ValueError("Missing, unknown or ambiguous source Enzyme_type for " + name)
            if len(systems) != 1 or systems[0] not in SYSTEM_MAP:
                raise ValueError("Missing, unknown or ambiguous source System for " + name)
            roles, types = list(ROLE_MAP[enzymes[0]]), [SYSTEM_MAP[systems[0]]]
            basis = "explicit source Enzyme_type/System; family hints only"
        rows.append(dict(**profile, roles=";".join(roles), type_hints=";".join(types),
                         off_target=str(is_off_target).lower(), library_sha256=inputs["hmm"]["sha256"],
                         native_enzyme_type=";".join(enzymes), native_system=";".join(systems),
                         metadata_lines=";".join(str(line) for line, _ in matches), role_mapping_basis=basis))
    # Validate everything before creating output; never replace a previous set.
    for key, path in source_paths.items():
        if digest(path) != inputs[key]["sha256"]:
            raise ValueError("Source changed while preparing: " + str(path))
    out.mkdir(parents=True)
    try:
        shutil.copyfile(source_paths["hmm"], out / "RM_HMMs.hmm")
        if digest(out / "RM_HMMs.hmm") != inputs["hmm"]["sha256"]:
            raise ValueError("Copied HMM library differs from source")
        with (out / "profile_roles.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        for key, path in source_paths.items():
            if digest(path) != inputs[key]["sha256"]:
                raise ValueError("Source changed during preparation: " + str(path))
        names = {row["profile_name"] for row in rows}
        receipt = dict(schema_version=1, operation="source_explicit_HMM_resource_preparation", passed=True,
                       created_at=datetime.now(timezone.utc).isoformat(), argv=list(argv or sys.argv),
                       exit_code=0, elapsed_seconds=time.monotonic()-started, inputs=inputs,
                       script_sha256=digest(__file__), profiles=len(rows), off_target_profiles=len(names & excluded),
                       metadata_profiles_not_in_library=sorted(set(by_name)-names),
                       off_target_names_not_in_library=sorted(excluded-names),
                       role_mapping={k: list(v) for k, v in ROLE_MAP.items()}, system_mapping=SYSTEM_MAP,
                       inputs_unchanged=True, original_GA_cutoffs_preserved=True,
                       detector_runs=0, biological_function_inferred=False,
                       outputs={name: digest(out/name) for name in ("RM_HMMs.hmm", "profile_roles.tsv")})
        (out / "preparation.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        return receipt
    except BaseException as error:
        (out / "failure.json").write_text(json.dumps(dict(passed=False, error=str(error), inputs=inputs), indent=2)+"\n")
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hmm", required=True, help="Authorized source HMMER3 library with native GA cutoffs")
    parser.add_argument("--metadata", required=True, help="Matching MicrobeMod restriction_metadata.csv")
    parser.add_argument("--off-target", required=True, help="Matching source off_target.txt")
    parser.add_argument("--out", required=True, help="New directory for the prepared HMM/metadata pair")
    args = parser.parse_args(argv)
    try:
        receipt = prepare(args.hmm, args.metadata, args.off_target, args.out,
                          argv=[parser.prog, *(argv if argv is not None else sys.argv[1:])])
    except (ValueError, OSError) as error:
        parser.exit(2, "prepare_hmm_resources: " + str(error) + "\n")
    print(json.dumps(dict(output=str(Path(args.out).absolute()), profiles=receipt["profiles"],
                          off_target_profiles=receipt["off_target_profiles"], passed=True)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
