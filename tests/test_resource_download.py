"""Synthetic HTTP/FTP streams exercise the real downloader and tar reader."""
import builtins
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tarfile
import urllib.request

import pytest

from defenseviz import resource_download as download

requires_posix = pytest.mark.skipif(
    sys.platform == "win32", reason="Managed provider cache requires POSIX flock; use Windows through WSL")


REBASE = b"\nREBASE protein sequences\n\n>REBASE:R.Test\tEnzType:Type II restriction enzyme\tSeqLength:4\nM A A A\n<>\n"
SOURCES = {
    "db/HMMs/RM_HMMs.hmm": b"HMMER3/f\nNAME test\nLENG 1\nGA 20 20;\n//\n",
    "db/restriction_metadata.csv": b"Name,Enzyme_type,System\ntest,RE,RM_Type_II\n",
    "db/HMMs/off_target.txt": b"unrelated_profile\n",
    "db/unneeded_reference.fasta": b">not-extracted\nMAAA\n",
}


def archive(files=None, extra=None):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as tar:
        for name, data in (files if files is not None else SOURCES).items():
            item = tarfile.TarInfo(name)
            item.size = len(data)
            tar.addfile(item, io.BytesIO(data))
        if extra is not None:
            tar.addfile(extra, io.BytesIO(b"x" * extra.size) if extra.isfile() else None)
    return stream.getvalue()


class Response(io.BytesIO):
    def __init__(self, data, url, status=200, headers=None):
        super().__init__(data)
        self.status = status
        self.headers = {"Content-Length": str(len(data)), **(headers or {})}
        self.url = url

    def geturl(self):
        return self.url


@pytest.fixture
def provider(monkeypatch):
    payload = archive()
    monkeypatch.setattr(download, "MICROBEMOD_SIZE", len(payload))
    monkeypatch.setattr(download, "MICROBEMOD_MD5", hashlib.md5(payload).hexdigest())
    state = {"archive": payload, "rebase": REBASE, "calls": [], "responses": {}}

    def metadata():
        return json.dumps(dict(doi=download.MICROBEMOD_DOI, license={"name": "MIT"}, files=[dict(
            id=55209713, size=download.MICROBEMOD_SIZE, computed_md5=download.MICROBEMOD_MD5)])).encode()

    def open_url(url):
        state["calls"].append(url)
        if url in state["responses"]:
            value = state["responses"][url]
            if isinstance(value, Exception):
                raise value
            return value()
        data = {download.MICROBEMOD_METADATA_URL: metadata(), download.MICROBEMOD_URL: state["archive"],
                download.REBASE_URL: state["rebase"]}[url]
        return Response(data, url, status=None if url == download.REBASE_URL else 200)

    monkeypatch.setattr(download, "_open", open_url)
    return state


def ready_pointer(cache):
    return cache / "providers/ready.json"


def pin_archive(monkeypatch, provider, payload):
    provider["archive"] = payload
    monkeypatch.setattr(download, "MICROBEMOD_SIZE", len(payload))
    monkeypatch.setattr(download, "MICROBEMOD_MD5", hashlib.md5(payload).hexdigest())


def test_public_provider_pins():
    assert download.MICROBEMOD_SIZE == 195132093
    assert download.MICROBEMOD_MD5 == "279f19f6d204fc4e1b750aea9828151d"
    assert download.MICROBEMOD_DOI == "10.6084/m9.figshare.24456481.v3"
    assert download.REBASE_URL == "ftp://ftp.neb.com/pub/rebase/protein_seqs.txt"


@requires_posix
def test_fresh_sources_integrity_and_entirely_offline_reuse(tmp_path, provider):
    result = download.fetch_provider_sources(tmp_path)
    assert set(result) == {"rebase_source", "hmm", "hmm_metadata", "off_target", "receipt", "source_label"}
    assert all(isinstance(value, str) for value in result.values())
    assert Path(result["rebase_source"]).read_bytes() == REBASE
    for key, relative in download.MEMBERS.items():
        assert Path(result[key]).read_bytes() == SOURCES[relative]
    base = Path(result["receipt"]).parent
    assert sorted(p.relative_to(base).as_posix() for p in (base / "db").rglob("*") if p.is_file()) == sorted(download.MEMBERS.values())
    receipt = json.loads(Path(result["receipt"]).read_text())
    assert receipt["ready"] is True and receipt["analysis_ready"] is False
    assert receipt["rebase"]["version"] == "unknown"
    assert receipt["rebase"]["transfer"]["tls_verified"] is False
    assert receipt["microbemod"]["transfer"]["tls_verified"] is True
    assert "unknown release" in result["source_label"]
    files_before = {p: p.read_bytes() for p in base.rglob("*") if p.is_file()}
    calls = list(provider["calls"])
    provider["responses"] = {u: OSError("offline") for u in calls}
    assert download.fetch_provider_sources(tmp_path) == result
    assert provider["calls"] == calls
    assert all(p.read_bytes() == value for p, value in files_before.items())


@pytest.mark.parametrize("key", ["hmm", "hmm_metadata", "off_target", "rebase_source", "receipt"])
@requires_posix
def test_corrupt_cache_never_silently_redownloads(tmp_path, provider, key):
    result = download.fetch_provider_sources(tmp_path)
    path = Path(result[key])
    path.write_bytes(path.read_bytes() + b"corruption")
    before = list(provider["calls"])
    with pytest.raises(ValueError, match="[Cc]hecksum|[Cc]orrupt"):
        download.fetch_provider_sources(tmp_path)
    assert provider["calls"] == before


@requires_posix
def test_refresh_preserves_old_snapshot_and_failed_refresh_keeps_ready_pointer(tmp_path, provider):
    old = download.fetch_provider_sources(tmp_path)
    old_bytes = {p: p.read_bytes() for p in Path(old["receipt"]).parent.rglob("*") if p.is_file()}
    provider["rebase"] = REBASE.replace(b"R.Test", b"R.New")
    new = download.fetch_provider_sources(tmp_path, refresh=True)
    assert new["receipt"] != old["receipt"]
    assert all(p.read_bytes() == value for p, value in old_bytes.items())
    pointer = ready_pointer(tmp_path).read_bytes()
    provider["responses"][download.REBASE_URL] = OSError("server unavailable")
    with pytest.raises(OSError, match="unavailable"):
        download.fetch_provider_sources(tmp_path, refresh=True)
    assert ready_pointer(tmp_path).read_bytes() == pointer
    assert download.fetch_provider_sources(tmp_path) == new


@pytest.mark.parametrize("bad", [b"<html>provider challenge</html>", b"", REBASE[:-4], REBASE.replace(b"M A A A", b"<script>")])
@requires_posix
def test_bad_rebase_never_becomes_ready(tmp_path, provider, bad):
    provider["rebase"] = bad
    with pytest.raises(ValueError):
        download.fetch_provider_sources(tmp_path)
    assert not ready_pointer(tmp_path).exists()
    assert list((tmp_path / "providers").glob(".partial-*/failure.json"))


@pytest.mark.parametrize("kind", ["md5", "short", "content_length", "partial_status", "challenge", "truncated_gzip"])
@requires_posix
def test_transfer_failures_never_publish(tmp_path, provider, monkeypatch, kind):
    if kind == "md5":
        provider["archive"] = provider["archive"][:-1] + bytes([provider["archive"][-1] ^ 1])
    elif kind == "short":
        provider["archive"] = provider["archive"][:-1]
    elif kind == "content_length":
        provider["responses"][download.MICROBEMOD_URL] = lambda: Response(provider["archive"][:-5], download.MICROBEMOD_URL,
                                                                             headers={"Content-Length": str(download.MICROBEMOD_SIZE)})
    elif kind == "partial_status":
        provider["responses"][download.MICROBEMOD_URL] = lambda: Response(provider["archive"], download.MICROBEMOD_URL, status=206)
    elif kind == "challenge":
        provider["responses"][download.MICROBEMOD_URL] = lambda: Response(b"", download.MICROBEMOD_URL, status=202)
    else:
        # Even a checksum-pinned truncated gzip must fail its format/CRC checks.
        pin_archive(monkeypatch, provider, provider["archive"][:-8])
    with pytest.raises((ValueError, EOFError, tarfile.TarError, OSError)):
        download.fetch_provider_sources(tmp_path)
    assert not ready_pointer(tmp_path).exists()


@pytest.mark.parametrize("kind", ["traversal", "absolute", "symlink", "hardlink", "device", "duplicate", "backslash"])
def test_unsafe_tar_rejected_including_unselected_members(tmp_path, kind):
    item = tarfile.TarInfo({"traversal": "db/../../escape", "absolute": "/escape", "backslash": "db\\escape",
                          "duplicate": "db/HMMs/RM_HMMs.hmm"}.get(kind, "db/unused"))
    item.size = 1
    if kind in {"symlink", "hardlink", "device"}:
        item.type = {"symlink": tarfile.SYMTYPE, "hardlink": tarfile.LNKTYPE, "device": tarfile.CHRTYPE}[kind]
        item.linkname = "../../escape"
        item.size = 0
    compressed = tmp_path / "provider.tar.gz"
    compressed.write_bytes(archive(extra=item))
    with pytest.raises(ValueError, match="[Aa]rchive"):
        download._extract(compressed, tmp_path / "extracted")
    assert not (tmp_path / "escape").exists()


@requires_posix
def test_missing_required_member_rejected(tmp_path, provider, monkeypatch):
    pin_archive(monkeypatch, provider, archive({k: v for k, v in SOURCES.items() if k != download.MEMBERS["hmm_metadata"]}))
    with pytest.raises(ValueError, match="Missing required"):
        download.fetch_provider_sources(tmp_path)
    assert not ready_pointer(tmp_path).exists()


@requires_posix
def test_cache_symlink_and_pointer_traversal_rejected_offline(tmp_path, provider):
    result = download.fetch_provider_sources(tmp_path)
    path = Path(result["hmm"])
    target = tmp_path / "same-bytes.hmm"
    target.write_bytes(path.read_bytes())
    path.unlink()
    path.symlink_to(target)
    calls = list(provider["calls"])
    with pytest.raises(ValueError, match="symlink"):
        download.fetch_provider_sources(tmp_path)
    ready_pointer(tmp_path).write_text(json.dumps({"snapshot": "../outside", "receipt_sha256": "0"*64}))
    with pytest.raises(ValueError, match="pointer"):
        download.fetch_provider_sources(tmp_path)
    assert provider["calls"] == calls


@requires_posix
def test_provider_metadata_must_match_pinned_release(tmp_path, provider):
    provider["responses"][download.MICROBEMOD_METADATA_URL] = lambda: Response(b'{"doi":"other","files":[]}', download.MICROBEMOD_METADATA_URL)
    with pytest.raises(ValueError, match="metadata differs"):
        download.fetch_provider_sources(tmp_path)
    assert provider["calls"] == [download.MICROBEMOD_METADATA_URL]


@requires_posix
def test_interrupted_stream_keeps_no_ready_file(tmp_path, provider):
    class Broken(Response):
        def read(self, n=-1):
            raise OSError("connection reset")
    provider["responses"][download.REBASE_URL] = lambda: Broken(REBASE, download.REBASE_URL)
    with pytest.raises(OSError, match="connection reset"):
        download.fetch_provider_sources(tmp_path)
    assert not ready_pointer(tmp_path).exists()
    assert not list((tmp_path / "providers").glob(".partial-*/REBASE_protein_seqs.txt"))


@requires_posix
def test_concurrent_setup_fails_without_network(tmp_path, provider):
    root = tmp_path / "providers"
    root.mkdir()
    with download._acquire_lock(root):
        with pytest.raises(RuntimeError, match="already active"):
            download.fetch_provider_sources(tmp_path)
    assert provider["calls"] == []


@requires_posix
def test_killed_downloader_lock_releases_and_partial_is_preserved(tmp_path, provider):
    script = """
from pathlib import Path
import sys, time
from defenseviz.resource_download import _acquire_lock
root=Path(sys.argv[1])/'providers'
root.mkdir()
lock=_acquire_lock(root)
partial=root/'.partial-interrupted'
partial.mkdir()
(partial/'REBASE_protein_seqs.txt.partial').write_bytes(b'unfinished transfer')
print('LOCKED',flush=True)
time.sleep(30)
"""
    process = subprocess.Popen([sys.executable, "-B", "-c", script, str(tmp_path)], stdout=subprocess.PIPE, text=True)
    try:
        assert process.stdout.readline().strip() == "LOCKED"
        with pytest.raises(RuntimeError, match="already active"):
            download.fetch_provider_sources(tmp_path)
        assert provider["calls"] == []
        process.kill()
        process.wait(timeout=5)
        result = download.fetch_provider_sources(tmp_path)
        assert Path(result["receipt"]).is_file()
        assert (tmp_path / "providers/.partial-interrupted/REBASE_protein_seqs.txt.partial").read_bytes() == b'unfinished transfer'
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        process.stdout.close()


@requires_posix
def test_symlink_lock_is_not_followed(tmp_path, provider):
    root = tmp_path / "providers"
    root.mkdir()
    target = tmp_path / "unrelated"
    target.write_text("untouched")
    (root / ".setup.lock").symlink_to(target)
    with pytest.raises(ValueError, match="lock must not be a symlink"):
        download.fetch_provider_sources(tmp_path)
    assert target.read_text() == "untouched"
    assert provider["calls"] == []


def test_redirect_never_downgrades_https():
    request = urllib.request.Request(download.MICROBEMOD_URL)
    with pytest.raises(ValueError, match="insecure"):
        download._HTTPSRedirect().redirect_request(request, None, 302, "Found", {}, "http://example.com/file")


def test_size_guard_and_signed_url_redaction(tmp_path, provider):
    url = download.MICROBEMOD_URL
    provider["responses"][url] = lambda: Response(b"abc", "https://provider.example/file?signature=do-not-log")
    with pytest.raises(ValueError, match="Content-Length"):
        download._download(url, tmp_path / "too-large", maximum=2)
    receipt = download._download(url, tmp_path / "small", maximum=3)
    assert receipt["final_url_without_query"] == "https://provider.example/file"
    assert "do-not-log" not in json.dumps(receipt)


@requires_posix
def test_rebase_requires_a_declared_complete_transfer_size(tmp_path, provider):
    def no_length():
        response = Response(REBASE, download.REBASE_URL, status=None)
        response.headers = {}
        return response
    provider["responses"][download.REBASE_URL] = no_length
    with pytest.raises(ValueError, match="complete resource size"):
        download.fetch_provider_sources(tmp_path)
    assert not ready_pointer(tmp_path).exists()


def test_duplicate_directory_alias_is_rejected(tmp_path):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w:gz") as tar:
        for name in ["db/", "./db"]:
            item = tarfile.TarInfo(name)
            item.type = tarfile.DIRTYPE
            tar.addfile(item)
    compressed = tmp_path / "provider.tar.gz"
    compressed.write_bytes(stream.getvalue())
    with pytest.raises(ValueError, match="duplicate"):
        download._extract(compressed, tmp_path / "extracted")


@requires_posix
def test_missing_cached_archive_is_an_offline_error(tmp_path, provider):
    result = download.fetch_provider_sources(tmp_path)
    (Path(result["receipt"]).parent / "microbemod_db_v3.tar.gz").unlink()
    calls = list(provider["calls"])
    with pytest.raises(ValueError, match="Missing provider cache"):
        download.fetch_provider_sources(tmp_path)
    assert provider["calls"] == calls


def test_missing_fcntl_gives_actionable_platform_error(tmp_path, monkeypatch):
    # On Windows this exercises the actual missing-module path. On POSIX,
    # simulate only the unavailable import, never alternate lock semantics.
    if sys.platform != "win32":
        original_import = builtins.__import__

        def without_fcntl(name, *args, **kwargs):
            if name == "fcntl":
                raise ImportError("fcntl unavailable on this platform")
            return original_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", without_fcntl)
    with pytest.raises(RuntimeError, match="Windows through WSL"):
        download._acquire_lock(tmp_path)
    assert not (tmp_path / ".setup.lock").exists()


@pytest.mark.parametrize("broken", [False, True])
def test_portable_extraction_checks_gzip_crc_and_exact_file_set(tmp_path, broken):
    compressed = tmp_path / "provider.tar.gz"
    data = archive()
    compressed.write_bytes(data[:-8] if broken else data)
    out = tmp_path / "extracted"
    if broken:
        with pytest.raises((EOFError, OSError, tarfile.TarError)):
            download._extract(compressed, out)
    else:
        download._extract(compressed, out)
        assert {p.relative_to(out).as_posix(): p.read_bytes() for p in out.rglob("*") if p.is_file()} == {
            name: SOURCES[name] for name in download.MEMBERS.values()}


@pytest.mark.parametrize("failure", ["md5", "short", "challenge"])
def test_portable_transfer_rejects_integrity_failures(tmp_path, provider, failure):
    target = tmp_path / "archive"
    expected = provider["archive"]
    if failure == "md5":
        provider["archive"] = expected[:-1] + bytes([expected[-1] ^ 1])
    elif failure == "short":
        provider["archive"] = expected[:-1]
    else:
        provider["responses"][download.MICROBEMOD_URL] = lambda: Response(b"", download.MICROBEMOD_URL, status=202)
    with pytest.raises(ValueError):
        download._download(download.MICROBEMOD_URL, target, size=len(expected), md5=hashlib.md5(expected).hexdigest())
    assert not target.exists()
