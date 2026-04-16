from pathlib import Path
import subprocess
import sys

import pytest

import gwaspeek.cli as cli_mod
from gwaspeek.cli import build_parser, main


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


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
    out = subprocess.check_output(cmd, text=True)
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
    out = subprocess.check_output(cmd, text=True, input="q\n")
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
    out = subprocess.check_output(cmd, text=True, input="q\n")
    assert "h help" in out and "v vars" in out and "t 37|38|off" in out


def test_cli_interactive_prints_loading_message(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(cli_mod, "run_interactive_manhattan", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_mod, "_app_version", lambda: "0.1.0")
    main(["-i", str(FIXTURE), "--width", "60", "--height", "18", "--ascii"])
    captured = capsys.readouterr()
    assert "[gwaspeek v0.1.0] Loading" in captured.err


def test_cli_positional_default_interactive_prints_loading_message(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cli_mod, "run_interactive_manhattan", lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_mod, "_app_version", lambda: "0.1.0")
    main([str(FIXTURE), "--width", "60", "--height", "18", "--ascii"])
    captured = capsys.readouterr()
    assert "[gwaspeek v0.1.0] Loading" in captured.err


def test_cli_default_skip_is_5() -> None:
    parser = build_parser()
    args = parser.parse_args([str(FIXTURE)])
    assert args.skip == 5.0
    assert args.chr_col is None and args.pos_col is None and args.p_col is None


def test_cli_rejects_multiple_inputs() -> None:
    with pytest.raises(SystemExit):
        main(["-s", str(FIXTURE), str(FIXTURE)])
