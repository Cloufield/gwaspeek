from pathlib import Path
import subprocess
import sys


FIXTURE = Path(__file__).parent / "fixtures" / "sumstats_small.tsv"


def test_cli_manhattan_smoke() -> None:
    cmd = [
        sys.executable,
        "-m",
        "gwast.cli",
        "manhattan",
        "--input",
        str(FIXTURE),
        "--width",
        "60",
        "--height",
        "18",
        "--ascii",
    ]
    out = subprocess.check_output(cmd, text=True)
    assert "Manhattan plot" in out


def test_cli_qq_smoke() -> None:
    cmd = [
        sys.executable,
        "-m",
        "gwast.cli",
        "qq",
        "--input",
        str(FIXTURE),
        "--width",
        "60",
        "--height",
        "18",
        "--ascii",
    ]
    out = subprocess.check_output(cmd, text=True)
    assert "QQ plot" in out
