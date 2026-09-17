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
                # 只处理普通文件
                if not path.is_file():
                    continue

                # 排除生成的压缩包本身
                if path.resolve() == OUTPUT:
                    continue

                # 只匹配 .R 和 .sh 文件
                if path.suffix not in {".R", ".sh"}:
                    continue

                # 排除大于 50 MB 的文件
                if path.stat().st_size > MAX_SIZE:
                    continue

                # 使用相对路径写入压缩包
                archive.add(path, arcname=path.relative_to(ROOT))

                print(f"已添加: {path.relative_to(ROOT)}")
        elif option == "result":
            for path in ROOT.rglob("*"):
                # 只处理普通文件
                if not path.is_file():
                    continue

                # 排除生成的压缩包本身
                if path.resolve() == OUTPUT:
                    continue

                # file target
                if path.suffix not in {".png", ".sh", ".R", ".csv"}:
                    continue

                # 排除大于 50 MB 的文件
                if path.stat().st_size > MAX_SIZE:
                    continue

                # 使用相对路径写入压缩包
                archive.add(path, arcname=path.relative_to(ROOT))

                print(f"已添加: {path.relative_to(ROOT)}")
        else:
            raise ValueError(f"Invalid option: {option}")

    print(f"\n完成，压缩包已生成：{OUTPUT}")


if __name__ == "__main__":
    main()
