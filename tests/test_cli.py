from pathlib import Path
import os
import subprocess
import sys

import pytest

import gwaspeek.cli as cli_mod
from gwaspeek.cli import build_parser, main


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"
SRC = Path(__file__).resolve().parents[1] / "src"


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    current = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = str(SRC) if not current else f"{SRC}{os.pathsep}{current}"
    return env


def test_cli_manhattan_smoke() -> None:
    cmd = [
        sys.executable,
        "-m",
        "gwaspeek.cli",
        "-s",
        str(FIXTURE),
        "--width",
        "60",
        "--height",
        "18",
        "--ascii",
    ]
    out = subprocess.check_output(cmd, text=True, env=_subprocess_env())
    assert "Manhattan" in out
    assert ":" in out


def test_cli_interactive_positional_same_as_dash_i() -> None:
    cmd = [
        sys.executable,
        "-m",
        "gwaspeek.cli",
        str(FIXTURE),
        "--width",
        "60",
        "--height",
        "18",
        "--ascii",
    ]
    out = subprocess.check_output(cmd, text=True, input="q\n", env=_subprocess_env())
    assert "h help" in out and "v vars" in out and "t 37|38|off" in out


def test_cli_manhattan_interactive_quit() -> None:
    cmd = [
        sys.executable,
        "-m",
        "gwaspeek.cli",
        "-i",
        str(FIXTURE),
        "--width",
        "60",
        "--height",
        "18",
        "--ascii",
    ]
    out = subprocess.check_output(cmd, text=True, input="q\n", env=_subprocess_env())
    assert "h help" in out and "v vars" in out and "t 37|38|off" in out


def test_cli_interactive_prints_loading_message(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli_mod, "run_interactive_manhattan", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_mod, "_app_version", lambda: "0.4.0")
    main(["-i", str(FIXTURE), "--width", "60", "--height", "18", "--ascii"])
    captured = capsys.readouterr()
    assert "[gwaspeek v0.4.0] Loading" in captured.err


def test_cli_positional_default_interactive_prints_loading_message(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cli_mod, "run_interactive_manhattan", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_mod, "_app_version", lambda: "0.4.0")
    main([str(FIXTURE), "--width", "60", "--height", "18", "--ascii"])
    captured = capsys.readouterr()
    assert "[gwaspeek v0.4.0] Loading" in captured.err


def test_cli_version_long_flag(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert "gwaspeek" in capsys.readouterr().out


def test_cli_version_short_flag(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as exc:
        main(["-v"])
    assert exc.value.code == 0
    assert "gwaspeek" in capsys.readouterr().out


def test_cli_default_skip_is_3() -> None:
    parser = build_parser()
    args = parser.parse_args([str(FIXTURE)])
    assert args.skip == 3.0
    assert args.build == "37"
    assert args.no_color is False
    assert args.chr_col is None and args.pos_col is None and args.p_col is None


def test_cli_rejects_multiple_inputs() -> None:
    with pytest.raises(SystemExit):
        main(["-s", str(FIXTURE), str(FIXTURE)])
