"""Tests must never dispatch into a workstation's installed native runtime."""
import pytest


@pytest.fixture(autouse=True)
def isolated_managed_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('DEFENSEVIZ_CACHE', str(tmp_path / 'unconfigured-managed-cache'))
