"""
Fetch WPP 2024 single-age mortality m(x,t) for a West-African REFERENCE GROUP,
used by the Li-Lee coherent model (Phase 6) to borrow strength for The Gambia.

Same pinned, reproducible source as fetch_wpp.py. One streamed pass over each
mx file keeps every reference-country row, so the download cost is one pass.
"""
from __future__ import annotations
import json
from pathlib import Path
import requests

PINNED_SHA = "2da7768ae64fc74105d3f9e98f9a74d37b62f99a"
RAW = f"https://raw.githubusercontent.com/PPgp/wpp2024/{PINNED_SHA}/data-raw/"

# Region: The Gambia + neighbours / Sahel-coastal West Africa with a comparable
# mortality regime. (UN M49 country codes.)
REF = {
    270: "Gambia", 686: "Senegal", 624: "Guinea-Bissau", 324: "Guinea",
    466: "Mali", 694: "Sierra Leone", 478: "Mauritania", 854: "Burkina Faso",
}
CODES = {str(c) for c in REF}

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "raw" / "wpp2024_refgroup"


def fetch(session, fname):
    header, kept = None, []
    with session.get(RAW + fname, stream=True, timeout=300) as r:
        r.raise_for_status()
        for i, line in enumerate(r.iter_lines(decode_unicode=True)):
            if line is None:
                continue
            if i == 0:
                header = line
                continue
            if line.split("\t", 1)[0] in CODES:
                kept.append(line)
    out = OUT / fname.replace(".txt", "_ref.tsv")
    out.write_text("\n".join([header] + kept) + "\n", encoding="utf-8")
    return len(kept)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    s = requests.Session(); s.headers.update({"User-Agent": "gambia-pop/0.1"})
    print(f"Reference group ({len(REF)}): {', '.join(REF.values())}")
    for f in ["mxB.txt", "mxF.txt", "mxM.txt"]:
        n = fetch(s, f)
        print(f"  {f}: {n} rows ({n // 101} countries x 101 ages)")
    (OUT / "_refgroup.json").write_text(
        json.dumps({"pinned_sha": PINNED_SHA, "countries": REF}, indent=2), encoding="utf-8")
    print(f"-> {OUT}")


if __name__ == "__main__":
    main()
