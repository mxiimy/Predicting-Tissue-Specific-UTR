import argparse
import re
import pandas as pd
from pathlib import Path

def parse_rnafold_out(filepath: Path) -> dict:
    mfe_map = {}
    if not filepath.exists():
        return mfe_map
    with open(filepath) as fh:
        lines = [l.rstrip() for l in fh if l.strip()]

    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i][1:]
            transcript_id = header.split("|")[0]
            i += 2 # Skip header and sequence lines
            if i < len(lines):
                match = re.search(r'\(\s*(-?\d+\.?\d*)\s*\)', lines[i])
                mfe_map[transcript_id] = float(match.group(1)) if match else None
            i += 1
        else:
            i += 1
    return mfe_map

def main():
    parser = argparse.ArgumentParser(description="Parse RNAfold output → CSV")
    parser.add_argument("--input",   required=True, help="Original input CSV")
    parser.add_argument("--rnafold", required=True, help="Dir with 5_UTR.out, CDS.out, 3_UTR.out")
    parser.add_argument("--output",  required=True, help="Path for output CSV")
    args = parser.parse_args()

    # Load original CSV matching the 8-column header structure
    print(f"[INFO] Loading {args.input}")
    df = pd.read_csv(
        args.input,
        header=0,
        names=["gene_id", "transcript_id", "5_UTR", "CDS", "3_UTR", "5_len", "CDS_len", "3_len"]
    )

    rnafold_dir = Path(args.rnafold)
    for region in ["5_UTR", "CDS", "3_UTR"]:
        out_file = rnafold_dir / f"{region}.out"
        print(f"[INFO] Parsing {out_file}")
        mfe_map = parse_rnafold_out(out_file)
        
        df[f"{region}_MFE"] = df["transcript_id"].map(mfe_map)
        print(f"  → {df[f'{region}_MFE'].notna().sum()} MFE values mapped")

    df.to_csv(args.output, index=False)
    print(f"[DONE] Saved → {args.output}")

if __name__ == "__main__":
    main()