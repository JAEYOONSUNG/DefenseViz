"""Acquire public reference snapshots; never run detectors or change models.

REBASE publishes its protein collection over anonymous FTP. That transport is
not authenticated: the recorded SHA-256 pins the bytes received, not their
authenticity or a historical release. MicrobeMod v3 is pinned by the provider's
file ID, byte size and MD5, then additionally recorded with SHA-256.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import tarfile
import time
import urllib.parse
import urllib.request
import uuid

REBASE_URL = "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt"
REBASE_SOURCE = "https://rebase.neb.com/rebase/rebase.seqs.html"
REBASE_TERMS = "https://www.neb.com/en-us/terms-of-use"
MICROBEMOD_URL = "https://api.figshare.com/v2/file/download/55209713"
MICROBEMOD_METADATA_URL = "https://api.figshare.com/v2/articles/24456481/versions/3"
MICROBEMOD_DOI = "10.6084/m9.figshare.24456481.v3"
MICROBEMOD_SIZE = 195132093
MICROBEMOD_MD5 = "279f19f6d204fc4e1b750aea9828151d"
MEMBERS = {
    "hmm": "db/HMMs/RM_HMMs.hmm",
    "hmm_metadata": "db/restriction_metadata.csv",
    "off_target": "db/HMMs/off_target.txt",
}
_FILES = {"rebase_source": "REBASE_protein_seqs.txt", **MEMBERS}
_MAX_DOWNLOAD = 2 * 1024**3
_MAX_UNPACKED = 4 * 1024**3
_BLOCK = 1024**2


def _now():
    return datetime.now(timezone.utc).isoformat()


def _digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm, usedforsecurity=algorithm != "md5")
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(_BLOCK), b""):
            h.update(block)
    return h.hexdigest()


def _atomic_json(path, value):
    temp = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    try:
        with temp.open("x", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, ensure_ascii=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


class _HTTPSRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise ValueError("Provider HTTPS redirect attempted an insecure transport")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open(url):
    # FTP is allowed only for the one publicly documented REBASE resource.
    if urllib.parse.urlsplit(url).scheme != "https" and url != REBASE_URL:
        raise ValueError("Unsupported provider URL")
    request = urllib.request.Request(url, headers={
        "User-Agent": "DefenseViz/1.0 reference-setup", "Accept-Encoding": "identity"})
    return urllib.request.build_opener(_HTTPSRedirect()).open(request, timeout=60)


def _download(url, destination, *, size=None, md5=None, maximum=_MAX_DOWNLOAD, require_length=False):
    """Stream into a fresh partial file; publish only after transfer integrity."""
    started, tick = _now(), time.monotonic()
    partial = destination.with_name(destination.name + ".partial")
    sha, legacy, received = hashlib.sha256(), hashlib.md5(usedforsecurity=False), 0
    with _open(url) as response, partial.open("xb") as target:
        status = getattr(response, "status", None)
        if status not in (None, 200):
            raise ValueError(f"Provider returned HTTP {status}; expected a complete file")
        headers = {k.lower(): v for k, v in response.headers.items()}
        if "content-range" in headers:
            raise ValueError("Unexpected partial provider response")
        if headers.get("content-encoding", "identity").lower() != "identity":
            raise ValueError("Unexpected HTTP content encoding")
        announced = headers.get("content-length")
        announced = int(announced) if announced is not None else None
        if require_length and announced is None:
            raise ValueError("Provider did not declare the complete resource size")
        if announced is not None and (announced <= 0 or announced > maximum):
            raise ValueError("Invalid provider Content-Length")
        if size is not None and announced is not None and announced != size:
            raise ValueError("Provider Content-Length differs from pinned size")
        while block := response.read(_BLOCK):
            received += len(block)
            if received > maximum or (size is not None and received > size):
                raise ValueError("Provider download exceeds its size limit")
            target.write(block)
            sha.update(block)
            legacy.update(block)
        if not received or (announced is not None and received != announced):
            raise ValueError("Incomplete provider download")
        if size is not None and received != size:
            raise ValueError("Provider download differs from pinned size")
        if md5 is not None and legacy.hexdigest() != md5:
            raise ValueError("Provider archive MD5 mismatch")
        target.flush()
        os.fsync(target.fileno())
        final = urllib.parse.urlsplit(response.geturl())
        # Do not persist short-lived cloud-storage signed query credentials.
        final_url = urllib.parse.urlunsplit((final.scheme, final.netloc, final.path, "", ""))
    os.replace(partial, destination)
    return dict(url=url, final_url_without_query=final_url, transport=urllib.parse.urlsplit(url).scheme,
                tls_verified=urllib.parse.urlsplit(url).scheme == "https", started_at=started,
                finished_at=_now(), elapsed_seconds=time.monotonic()-tick, size=received,
                sha256=sha.hexdigest(), md5=legacy.hexdigest(),
                response_headers={k: headers[k] for k in ("content-length", "content-type", "etag", "last-modified") if k in headers})


def _validate_rebase(path):
    """Check the provider record framing; downstream native parsing is separate."""
    records, residues, current, ended = 0, 0, False, False
    inline, standalone, empty = 0, 0, 0
    declared_zero = False
    with path.open(encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if not current and not records and line == "REBASE protein sequences":
                continue
            if line.startswith(">REBASE:"):
                if current and not ended:
                    raise ValueError(f"Incomplete REBASE record before line {number}")
                if not re.match(r">REBASE:\S+\s+.*EnzType:", line):
                    raise ValueError(f"Malformed REBASE record header at line {number}")
                lengths = re.findall(r"(?:^|\s)SeqLength:([^\s]+)", line)
                declared_zero = lengths == ["0"]
                records += 1
                residues, current, ended = 0, True, False
            elif line == "<>":
                if not current or ended or (not residues and not declared_zero):
                    raise ValueError(f"Malformed REBASE record terminator at line {number}")
                empty += not residues
                standalone += 1
                ended = True
            elif current and not ended:
                # Official records may attach <> to the final sequence line.
                # Only a terminal marker is valid; embedded/repeated markers
                # still fail the sequence grammar. Preserve source bytes.
                terminal = line.endswith("<>")
                sequence = line[:-2].rstrip() if terminal else line
                if declared_zero or not re.fullmatch(r"[A-Za-z* .-]+", sequence):
                    raise ValueError(f"Unexpected content in REBASE protein collection at line {number}")
                residues += len("".join(sequence.split()))
                ended = terminal
                inline += terminal
            else:
                raise ValueError(f"Unexpected content in REBASE protein collection at line {number}")
    if not records or not ended or records == empty:
        raise ValueError("Empty or truncated REBASE protein collection")
    return {"records": records, "records_with_sequence": records-empty,
            "explicit_empty_records": empty, "inline_terminators": inline,
            "standalone_terminators": standalone, "framing_valid": True,
            "semantic_validation": "native prepare-resources required"}


def _extract(archive, out):
    """Write exactly the three regular reference files, never tar.extract()."""
    found, seen, total = set(), set(), 0
    with gzip.open(archive, "rb") as compressed:
        with tarfile.open(fileobj=compressed, mode="r|") as bundle:
            for member in bundle:
                name = member.name
                if name.startswith("./"):
                    name = name[2:]
                path = PurePosixPath(name)
                if (not name or path.is_absolute() or ".." in path.parts or "\\" in name
                        or ":" in name or str(path) in seen or str(path) != name.rstrip("/")):
                    raise ValueError("Unsafe or duplicate archive member: " + member.name)
                seen.add(str(path))
                if not (member.isdir() or member.isfile()) or member.issparse():
                    raise ValueError("Archive contains a link, device or unsupported member: " + name)
                total += member.size
                if member.size < 0 or total > _MAX_UNPACKED:
                    raise ValueError("Archive exceeds unpacked size limit")
                if name not in MEMBERS.values():
                    continue
                if not member.isfile() or not 0 < member.size <= 256 * _BLOCK:
                    raise ValueError("Invalid required provider file: " + name)
                dest = out.joinpath(*path.parts)
                dest.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                if source is None:
                    raise ValueError("Unreadable archive member: " + name)
                received = 0
                with source, dest.open("xb") as target:
                    while block := source.read(_BLOCK):
                        received += len(block)
                        target.write(block)
                if received != member.size:
                    raise ValueError("Truncated archive member: " + name)
                found.add(name)
        # tar EOF precedes gzip EOF; drain to verify CRC and truncated trailers.
        while compressed.read(_BLOCK):
            pass
    if found != set(MEMBERS.values()):
        raise ValueError("Missing required MicrobeMod files: " + ", ".join(sorted(set(MEMBERS.values())-found)))
    with (out / MEMBERS["hmm"]).open("rb") as stream:
        if not stream.read(16).startswith(b"HMMER3/"):
            raise ValueError("Invalid HMM library signature")
    with (out / MEMBERS["hmm_metadata"]).open(encoding="utf-8-sig", newline="") as stream:
        if not {"Name", "Enzyme_type", "System"} <= set(csv.DictReader(stream).fieldnames or []):
            raise ValueError("Invalid MicrobeMod metadata schema")
    (out / MEMBERS["off_target"]).read_text(encoding="utf-8")


def _regular(root, relative):
    parts = PurePosixPath(relative)
    if parts.is_absolute() or ".." in parts.parts or str(parts) != relative or "\\" in relative:
        raise ValueError("Unsafe cache receipt path")
    path = root
    for part in parts.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError("Provider cache contains a symlink")
    if not path.is_file():
        raise ValueError("Missing provider cache file: " + str(path))
    return path


def _reuse(root):
    pointer = json.loads(_regular(root, "ready.json").read_text())
    snapshot = pointer["snapshot"]
    if not re.fullmatch(r"snapshots/[0-9a-f]{32}", snapshot):
        raise ValueError("Invalid provider snapshot pointer")
    receipt_path = _regular(root, snapshot + "/receipt.json")
    if _digest(receipt_path) != pointer["receipt_sha256"]:
        raise ValueError("Provider receipt checksum mismatch")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("ready") is not True or receipt.get("schema_version") != 1:
        raise ValueError("Provider snapshot is not ready")
    expected = set(_FILES.values()) | {"microbemod_db_v3.tar.gz", "microbemod_metadata.json"}
    if set(receipt["files"]) != expected:
        raise ValueError("Provider receipt has an unexpected file set")
    base = root / snapshot
    for name, item in receipt["files"].items():
        path = _regular(base, name)
        if path.stat().st_size != item["size"] or _digest(path) != item["sha256"]:
            raise ValueError("Corrupt provider cache: " + str(path) + "; use refresh=True for a new snapshot")
    archive = base / "microbemod_db_v3.tar.gz"
    if archive.stat().st_size != MICROBEMOD_SIZE or _digest(archive, "md5") != MICROBEMOD_MD5:
        raise ValueError("Cached MicrobeMod archive differs from pinned provider release")
    return {**{key: str(base / name) for key, name in _FILES.items()},
            "receipt": str(receipt_path), "source_label": receipt["source_label"]}


def _acquire_lock(root):
    """The kernel releases this lock even if a downloader is killed abruptly."""
    try:
        import fcntl
    except ImportError as error:
        raise RuntimeError("Managed resource setup requires Linux, macOS or Windows through WSL") from error
    path = root / ".setup.lock"
    if path.is_symlink():
        raise ValueError("Provider setup lock must not be a symlink")
    descriptor = os.open(path, os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0), 0o600)
    stream = os.fdopen(descriptor, "r+")
    try:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        stream.close()
        raise RuntimeError("Provider setup is already active in " + str(root)) from error
    except BaseException:
        stream.close()
        raise
    # Retain the lock inode; unlinking it would let another process lock a
    # different inode while the first process is still using this one.
    return stream


def fetch_provider_sources(cache_root, refresh=False):
    """Return six string fields pointing to a verified persistent source set.

    Valid reuse is entirely offline. Corrupt ready caches fail rather than
    silently acquiring a different REBASE snapshot. Explicit refresh creates a
    new snapshot and leaves every previous source file unchanged. These are raw
    provider resources; analysis readiness still requires prepare-resources.
    """
    if type(refresh) is not bool:
        raise TypeError("refresh must be a bool")
    started, tick = _now(), time.monotonic()
    root = Path(cache_root).expanduser().resolve() / "providers"
    root.mkdir(parents=True, exist_ok=True)
    if root.is_symlink():
        raise ValueError("Provider cache root must not be a symlink")
    lock = _acquire_lock(root)
    stage = None
    try:
        if (root / "ready.json").exists() or (root / "ready.json").is_symlink():
            if not refresh:
                return _reuse(root)
        identifier = uuid.uuid4().hex
        stage = root / (".partial-" + identifier)
        stage.mkdir()
        metadata_transfer = _download(MICROBEMOD_METADATA_URL, stage / "microbemod_metadata.json", maximum=2*_BLOCK)
        metadata = json.loads((stage / "microbemod_metadata.json").read_text())
        records = [x for x in metadata.get("files", []) if x.get("id") == 55209713]
        if (metadata.get("doi") != MICROBEMOD_DOI or len(records) != 1
                or records[0].get("size") != MICROBEMOD_SIZE
                or records[0].get("computed_md5") != MICROBEMOD_MD5):
            raise ValueError("Figshare metadata differs from pinned MicrobeMod v3")
        archive = stage / "microbemod_db_v3.tar.gz"
        microbe = _download(MICROBEMOD_URL, archive, size=MICROBEMOD_SIZE, md5=MICROBEMOD_MD5)
        _extract(archive, stage)
        rebase = _download(REBASE_URL, stage / _FILES["rebase_source"], require_length=True)
        framing = _validate_rebase(stage / _FILES["rebase_source"])
        label = "REBASE unknown release; retrieved " + rebase["finished_at"] + "; SHA256 " + rebase["sha256"] + "; MicrobeMod " + MICROBEMOD_DOI
        names = set(_FILES.values()) | {archive.name, "microbemod_metadata.json"}
        receipt = dict(schema_version=1, ready=True, operation="public_provider_source_acquisition",
            created_at=_now(), started_at=started, elapsed_seconds=time.monotonic()-tick,
            acquisition_code_sha256=_digest(__file__), source_label=label,
            files={name: dict(size=(stage/name).stat().st_size, sha256=_digest(stage/name)) for name in sorted(names)},
            rebase=dict(provider="REBASE / New England Biolabs", version="unknown", source_page=REBASE_SOURCE,
                        terms_url=REBASE_TERMS, redistribution_license="not established; retain provider terms",
                        transport_authenticity="unauthenticated FTP; SHA256 records acquired bytes only",
                        transfer=rebase, format_check=framing),
            microbemod=dict(provider="Alex Crits-Christoph / MicrobeMod", doi=MICROBEMOD_DOI,
                           provider_license=metadata.get("license"), file_id=55209713,
                           expected_size=MICROBEMOD_SIZE, expected_md5=MICROBEMOD_MD5,
                           transfer=microbe, metadata_transfer=metadata_transfer),
            detector_runs=0, analysis_ready=False)
        _atomic_json(stage / "receipt.json", receipt)
        snapshots = root / "snapshots"
        snapshots.mkdir(exist_ok=True)
        if snapshots.is_symlink():
            raise ValueError("Snapshot directory must not be a symlink")
        dest = snapshots / identifier
        stage.rename(dest)
        stage = None
        _atomic_json(root / "ready.json", dict(snapshot="snapshots/"+identifier,
                                              receipt_sha256=_digest(dest/"receipt.json")))
        return _reuse(root)
    except BaseException as error:
        if stage is not None and stage.is_dir():
            _atomic_json(stage / "failure.json", dict(ready=False, failed_at=_now(), error=str(error)))
        raise
    finally:
        lock.close()
