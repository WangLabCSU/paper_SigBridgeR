#!/usr/bin/env python3

import tarfile
from pathlib import Path

ROOT: Path = Path(".").cwd()
OUTPUT: Path = ROOT / "paper_SigBridgeR.tar.gz"
MAX_SIZE: int = 50 * 1024 * 1024  # 50 MB
OPTION: list = ["script", "result"]


def main(option: str | None):
    with tarfile.open(OUTPUT, mode="w:gz") as archive:
        if option == "script":
            for path in ROOT.rglob("*"):
                if not path.is_file():
                    continue

                if path.resolve() == OUTPUT:
                    continue

                if path.suffix not in {
                    ".R",
                    ".sh",
                    ".py",
                    ".Rproj",
                    ".md",
                    ".toml",
                    ".here",
                }:
                    continue

                if path.stat().st_size > MAX_SIZE:
                    continue

                archive.add(path, arcname=path.relative_to(ROOT))

                print(f"已添加: {path.relative_to(ROOT)}")
        elif option == "result":
            for path in ROOT.rglob("*"):
                if not path.is_file():
                    continue

                if path.resolve() == OUTPUT:
                    continue

                if path.suffix not in {".png", ".sh", ".R", ".csv"}:
                    continue

                if path.stat().st_size > MAX_SIZE:
                    continue

                archive.add(path, arcname=path.relative_to(ROOT))

                print(f"已添加: {path.relative_to(ROOT)}")
        else:
            raise ValueError(f"Invalid option: {option}")

    print(f"\n完成，压缩包已生成：{OUTPUT}")


if __name__ == "__main__":
    main("script")
