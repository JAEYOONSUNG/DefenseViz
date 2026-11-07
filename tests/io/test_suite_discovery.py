"""Project declarations are discovery boundaries, not extra primary samples."""
import pytest

from defenseviz.io.suite import load_suite
from defenseviz.io.suite_discovery import discover_suite
from tests.test_suite_discovery import native


@pytest.mark.parametrize("project_has_modules", [False, True])
def test_nested_project_is_skipped_before_native_modules(tmp_path, project_has_modules):
    keep = native(tmp_path / "primary")
    project = tmp_path / "independent-case"
    child = native(project if project_has_modules else project / "source")
    declaration = project / "defenseviz_project.json"
    declaration.write_text('{"schema_version":1}', encoding="utf-8")
    inventory = discover_suite(tmp_path)
    assert [row["directory"] for row in inventory["samples"]] == [str(keep)]
    assert {"path": str(project), "reason": "nested_project"} in inventory["skipped"]
    # Explicitly selecting the project itself still discovers intended outputs.
    assert [row["directory"] for row in discover_suite(project)["samples"]] == [str(child)]
    data = load_suite(tmp_path)
    assert [row["sample_id"] for row in data["samples"]] == [keep.name]
    assert all(str(project) not in path for path in data["input_sha256"])


def test_root_project_declaration_does_not_disable_discovery(tmp_path):
    (tmp_path / "defenseviz_project.json").write_text('{}', encoding="utf-8")
    expected = {str(native(tmp_path / name)) for name in ("first", "second")}
    assert {row["directory"] for row in discover_suite(tmp_path)["samples"]} == expected


def test_malformed_nested_declaration_still_prevents_accidental_intake(tmp_path):
    native(tmp_path / "kept")
    project = tmp_path / "unfinished-project"
    native(project / "source")
    (project / "defenseviz_project.json").write_text('not JSON', encoding="utf-8")
    inventory = discover_suite(tmp_path)
    assert [row["spec"]["sample_id"] for row in inventory["samples"]] == ["kept"]
    assert any(row["reason"] == "nested_project" for row in inventory["skipped"])


@pytest.mark.parametrize("outside", [False, True])
@pytest.mark.parametrize("link_descendant", [False, True])
def test_symlink_cannot_bypass_project_boundary(tmp_path, outside, link_descendant):
    root = tmp_path / "cohort"
    native(root / "kept")
    project = (tmp_path if outside else root) / "independent"
    child = native(project / "source")
    (project / "defenseviz_project.json").write_text('{}', encoding="utf-8")
    alias = root / "a-alias"
    try:
        alias.symlink_to(child if link_descendant else project, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks unavailable for this user/platform")
    inventory = discover_suite(root)
    assert [row["spec"]["sample_id"] for row in inventory["samples"]] == ["kept"]
    reason = "nested_project_descendant" if link_descendant else "nested_project"
    assert {"path": str(alias), "reason": reason} in inventory["skipped"]


def test_explicitly_selected_root_project_symlink_remains_eligible(tmp_path):
    project = tmp_path / "selected"
    child = native(project / "source")
    (project / "defenseviz_project.json").write_text('{}', encoding="utf-8")
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(project, target_is_directory=True)
    except OSError:
        pytest.skip("Directory symlinks unavailable for this user/platform")
    assert [row["directory"] for row in discover_suite(alias)["samples"]] == [str(child)]
