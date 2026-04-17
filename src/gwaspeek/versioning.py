from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


def package_version() -> str:
    try:
        return version("gwaspeek")
    except PackageNotFoundError:
        return "dev"
