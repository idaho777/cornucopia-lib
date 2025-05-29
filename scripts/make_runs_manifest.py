# scripts/make_runs_manifest.py
import json
from pathlib import Path


def make_runs_manifest(runs_dir: Path, out_file: Path) -> None:
    """
    Scan runs_dir for all subdirectories named 'fits-*'
    and write their names to out_file as a JSON list.
    """
    print(f"Scanning {runs_dir} for runs...", end="")

    runs: list[str] = sorted(
        [
            d.name
            for d in runs_dir.iterdir()
            if d.is_dir() and d.name.startswith("fits-")
        ]
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(runs, indent=2), encoding="utf-8")

    print("Complete")


if __name__ == "__main__":
    make_runs_manifest(runs_dir=Path("./runs"), out_file=Path("./runs/runs.json"))
