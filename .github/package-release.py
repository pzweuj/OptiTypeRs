"""Bundle a release binary with its project and dependency notices."""

from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import zipfile


def package(platform: str, executable: str, version: str) -> None:
    vendor = Path("target/release-vendor")
    release = Path("release")
    release.mkdir(exist_ok=True)
    name = f"optitype-rs-v{version}-{platform}"
    with tempfile.TemporaryDirectory() as temporary:
        bundle = Path(temporary) / name
        bundle.mkdir()
        shutil.copy2(Path("target/release") / executable, bundle / executable)
        for license_file in ("LICENSE", "LICENSE-OptiType"):
            shutil.copy2(license_file, bundle / license_file)

        notices = sorted(
            file for file in vendor.rglob("*")
            if file.is_file() and file.name.upper().startswith(
                ("LICENSE", "LICENCE", "COPYING", "COPYRIGHT", "NOTICE", "UNLICENSE")
            ) and file.suffix.lower() not in (".rs", ".c", ".h", ".cpp")
        )
        if not notices:
            raise RuntimeError("No dependency license notices found")
        for file in notices:
            destination = bundle / "dependency-notices" / file.relative_to(vendor)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(file, destination)

        if platform.startswith("windows"):
            with zipfile.ZipFile(release / f"{name}.zip", "w", zipfile.ZIP_DEFLATED) as archive:
                for file in sorted(bundle.rglob("*")):
                    if file.is_file():
                        archive.write(file, file.relative_to(bundle.parent))
        else:
            with tarfile.open(release / f"{name}.tar.gz", "w:gz") as archive:
                archive.add(bundle, arcname=name)


if __name__ == "__main__":
    package(*sys.argv[1:])
