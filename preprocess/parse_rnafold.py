import argparse
import re
import pandas as pd
from pathlib import Path

def parse_rnafold_out(filepath: Path) -> dict:
    mfe_map = {}
    if not filepath.exists():
        return mfe_map
        
    print(f"  Reading {filepath.name}...")
    with open(filepath) as fh:
        # Process line by line to save memory and allow progress tracking
        current_id = None
        count = 0
        for line in fh:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                current_id = line[1:].split("|")[0]
            elif "(" in line and ")" in line:
                match = re.search(r'\(\s*(-?\d+\.?\d*)\s*\)', line)
                if match and current_id:
                    mfe_map[current_id] = float(match.group(1))
                    count += 1
                    if count % 1000 == 0:
                        print(f"    ... parsed {count} MFE values")
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