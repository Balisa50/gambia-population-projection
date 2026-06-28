"""
Fetch the UN World Population Prospects 2024 inputs for The Gambia.

Source: the UN Population Division's official `wpp2024` R data package
(github.com/PPgp/wpp2024), whose `data-raw/` directory holds the published WPP
2024 series as plain tab-separated text. We pin to a specific commit SHA so the
download is exactly reproducible, stream each (sometimes ~50 MB) file, and keep
only The Gambia's rows (country_code 270). The full multi-country files are
never written to disk - only the small Gambia extracts - which also keeps us
within tight disk limits.

Run:  python src/fetch_wpp.py
Out:  data/raw/wpp2024/*.tsv  +  data/raw/wpp2024/_manifest.json
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests

# --- reproducibility: pin to a commit, not a moving branch -------------------
REPO = "PPgp/wpp2024"
PINNED_SHA = "2da7768ae64fc74105d3f9e98f9a74d37b62f99a"  # main @ 2025-06-24
RAW = f"https://raw.githubusercontent.com/{REPO}/{PINNED_SHA}/data-raw/"

GAMBIA_CODE = "270"  # UN M49 / country_code for The Gambia

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "raw" / "wpp2024"

# Curated set of source files, grouped by role in the analysis.
FILES = {
    # Age-specific central death rates m(x,t), single age, 1950-2023 - the
    # Lee-Carter fitting series.
    "mortality": ["mxF.txt", "mxM.txt", "mxB.txt"],
    # Population by single age - exposures for the Poisson model + projection base.
    "population": ["popF.txt", "popM.txt"],
    # Life expectancy at birth (annual estimates) - headline + sanity check.
    "e0_estimates": ["e0B.txt", "e0F.txt", "e0M.txt"],
    # UN's official probabilistic e0 projection (median + 80/95% bounds) -
    # the benchmark our independent forecast is compared against (RQ4).
    "e0_projection": [
        "e0Bproj.txt", "e0Fproj.txt", "e0Mproj.txt",
        "e0Bproj80l.txt", "e0Bproj80u.txt", "e0Bproj95l.txt", "e0Bproj95u.txt",
    ],
    # Fertility inputs + projection (for the cohort-component model).
    "fertility": [
        "tfr.txt", "tfrprojMed.txt", "percentASFR.txt",
        "tfrproj80l.txt", "tfrproj80u.txt", "tfrproj95l.txt", "tfrproj95u.txt",
    ],
    # Sex ratio at birth.
    "srb": ["sexRatio.txt"],
    # Net migration (estimates + medium projection) - explicit scenarios later.
    "migration": ["migration.txt", "migrationprojMed.txt"],
    # UN projected population by single age & sex (medium) - benchmark for our
    # own projection.
    "population_projection": ["popFprojMed.txt", "popMprojMed.txt"],
    # Counts / crude rates for cross-checks.
    "counts": ["deaths.txt", "births.txt", "cdr.txt", "cbr.txt"],
}


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def fetch_one(session: requests.Session, fname: str) -> dict:
    """Stream one data-raw file and keep header + The Gambia's rows."""
    url = RAW + fname
    header = None
    kept: list[str] = []
    with session.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        for i, line in enumerate(r.iter_lines(decode_unicode=True)):
            if line is None:
                continue
            if i == 0:
                header = line
                continue
            # country_code is the first tab-delimited field
            if line.split("\t", 1)[0] == GAMBIA_CODE:
                kept.append(line)

    if header is None:
        raise RuntimeError(f"empty/invalid file: {fname}")

    out_text = "\n".join([header] + kept) + "\n"
    out_path = OUT_DIR / fname.replace(".txt", "_GMB.tsv")
    out_path.write_text(out_text, encoding="utf-8")

    n_cols = len(header.split("\t"))
    return {
        "file": fname,
        "url": url,
        "saved_as": out_path.name,
        "n_data_rows": len(kept),
        "n_columns": n_cols,
        "sha256": _sha256(out_text),
        "bytes": len(out_text.encode("utf-8")),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "gambia-pop-research/0.1"})

    manifest = {
        "source_repo": REPO,
        "pinned_sha": PINNED_SHA,
        "raw_base": RAW,
        "country_code": GAMBIA_CODE,
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "files": [],
    }

    print(f"Fetching WPP 2024 (Gambia, code {GAMBIA_CODE}) @ {PINNED_SHA[:10]}")
    print(f"{'file':<22}{'rows':>6}{'cols':>6}  status")
    print("-" * 50)
    for group, files in FILES.items():
        for fname in files:
            try:
                rec = fetch_one(session, fname)
                rec["group"] = group
                manifest["files"].append(rec)
                print(f"{fname:<22}{rec['n_data_rows']:>6}{rec['n_columns']:>6}  ok")
            except Exception as e:  # noqa: BLE001 - report and continue
                print(f"{fname:<22}{'':>6}{'':>6}  FAIL: {e!r}")
                manifest["files"].append({"file": fname, "group": group, "error": repr(e)})

    (OUT_DIR / "_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    ok = sum(1 for f in manifest["files"] if "error" not in f)
    print("-" * 50)
    print(f"Done: {ok}/{len(manifest['files'])} files -> {OUT_DIR}")


if __name__ == "__main__":
    main()
