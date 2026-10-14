from __future__ import annotations

import json
from pathlib import Path

import pytest

import cellpilot.cli as cli
from cellpilot.provenance import verify_manifest


def test_interpret_command_updates_status_history_and_manifest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = tmp_path / "project"
    run_dir = root / "results" / "existing-run"
    run_dir.mkdir(parents=True)
    (root / "pyproject.toml").write_text("[project]\nname='fixture'\n", encoding="utf-8")
    config_path = root / "config.yaml"
    config_path.write_text("project_name: fixture\n", encoding="utf-8")

    def fake_interpret(_root: Path, _config: object, selected_run: Path) -> None:
        target = selected_run / "interpretation" / "interpretation.json"
        target.parent.mkdir(parents=True)
        target.write_text('{"metadata":{"external_api_called":false}}\n', encoding="utf-8")

    monkeypatch.setattr(cli, "_run_interpretation", fake_interpret)
    exit_code = cli.main(
        [
            "--project-root",
            str(root),
            "interpret",
            "--config",
            str(config_path),
            "--run-dir",
            str(run_dir),
        ]
    )

    assert exit_code == 0
    status = json.loads((run_dir / "stage_status.json").read_text(encoding="utf-8"))
    assert status["stages"]["interpretation"]["status"] == "completed"
    assert status["stages"]["manifest"]["status"] == "completed"
    history = json.loads((run_dir / "execution_history.json").read_text(encoding="utf-8"))
    assert history["events"][-1]["operation"] == "interpret"
    report = verify_manifest(run_dir)
    assert report.valid
    assert report.checked_inputs == 1


def test_changed_historical_config_prevents_manifest_rewrite(tmp_path: Path) -> None:
    root = tmp_path / "project"
    root.mkdir()
    config_path = root / "config.yaml"
    config_path.write_text("project_name: first\n", encoding="utf-8")
    run_dir = root / "results" / "run"
    run_dir.mkdir(parents=True)
    cli._record_execution(
        run_dir,
        root,
        operation="all",
        config_path=config_path,
    )

    config_path.write_text("project_name: changed\n", encoding="utf-8")

    with pytest.raises(ValueError, match="earlier run operation has changed"):
        cli._collect_run_inputs(root, run_dir)
