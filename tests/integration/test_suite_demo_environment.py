"""Synthetic demo must not depend on configured genome analysis resources."""
import json
import shlex


def test_suite_demo_uses_stored_examples_with_engine_configured(tmp_path, monkeypatch, capsys):
    from defenseviz.cli import main
    engine_config = tmp_path / "engine.json"
    engine_config.write_text("invalid engine configuration must not be read")
    monkeypatch.setenv("DEFENSEVIZ_ENGINE_CONFIG", str(engine_config))
    output = tmp_path / "demo"
    code = main(["demo", "--suite", "--out", str(output)])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    response = json.loads(captured.out)
    assert response["samples"] == 3
    receipt = json.loads((output / "results/receipt.json").read_text())
    assert receipt["options"]["evaluation"] == "stored"
    assert receipt["evaluation_execution"]["rm_mode"] == "stored_replay"
    assert receipt["evaluation_execution"]["fresh_full_execution"] is False
    assert receipt["detector_runs"] == receipt["hmm_search_runs"] == receipt["homology_search_runs"] == 0
    assert receipt["operation"] == "saved_suite_result_import_aggregation_and_render"
    assert (output / "results/figures/publication.pdf").stat().st_size > 1000
    assert list((output / "results/figures").glob("*.svg"))

    # Follow the generated instructions under the same configured environment.
    # This checks the user-visible command, not a separately constructed replay.
    readme = (output / "inputs/README.txt").read_text()
    command = shlex.split(readme.split("reproduce with: ", 1)[1])
    assert command[0] == "defenseviz"
    monkeypatch.chdir(output / "inputs")
    assert main(command[1:]) == 0
    rerun = json.loads((output / "rerun/receipt.json").read_text())
    assert rerun["options"]["evaluation"] == "stored"
    assert rerun["detector_runs"] == rerun["hmm_search_runs"] == rerun["homology_search_runs"] == 0
    original = json.loads((output / "results/analysis.json").read_text())
    repeated = json.loads((output / "rerun/analysis.json").read_text())
    assert repeated["global_summary"] == original["global_summary"]
