import argparse
import re
import pandas as pd
from pathlib import Path

def convert_u_to_t(sequence):
    if isinstance(sequence, str):
        return sequence.replace("U", "T")
    return sequence

def infer_region(path: Path) -> str | None:
    match = re.search(r"(5_UTR|CDS|3_UTR)", path.as_posix())
    return match.group(1) if match else None

def parse_rnafold_out(filepath: Path) -> dict:
    mfe_map = {}
    if not filepath.exists():
        return mfe_map
        
    print(f"  Reading {filepath.name}...")
    with open(filepath) as fh:
        current_id = None
        count = 0
        for line in fh:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                # Extracts ID before the first pipe symbol
                current_id = line[1:].split("|")[0]
            elif "(" in line and ")" in line:
                match = re.search(r'\(\s*(-?\d+\.?\d*)\s*\)', line)
                if match and current_id:
                    mfe_map[current_id] = float(match.group(1))
                    count += 1
                    if count % 1000 == 0:
                        print(f"    ... parsed {count} MFE values")
    return mfe_map

def parse_rnafold_root(root_dir: Path) -> dict:
    region_maps = {"5_UTR": {}, "CDS": {}, "3_UTR": {}}
    resolved_path = root_dir.resolve()
    
    if not resolved_path.exists():
        print(f"[ERROR] RNAfold results directory not found: {resolved_path}")
        return region_maps

    print(f"[INFO] Searching for .out files in: {resolved_path}")
    for out_file in sorted(resolved_path.rglob("*.out")):
        region = infer_region(out_file)
        if region is None:
            continue
        print(f"[INFO] Parsing {out_file.name} (Region: {region})")
        region_maps[region].update(parse_rnafold_out(out_file))

    return region_maps

def main():
    parser = argparse.ArgumentParser(description="Map RNAfold MFE values to Sequence Database")
    
    parser.add_argument("-i", "--input", 
                        default="./data/initialDB_sequences.csv", 
                        help="Path to input CSV (default: ./data/initialDB_sequences.csv)")
    
    parser.add_argument("-o", "--output", 
                        default="./data/mfe_results.csv", 
                        help="Path to output CSV (default: ./data/mfe_results.csv)")
    
    parser.add_argument("-r", "--rnafold", 
                        default="./results", 
                        help="RNAfold results directory (default: ./results)")
    
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}")
        return

    print(f"[INFO] Loading {input_path}")
    df = pd.read_csv(
        input_path,
        header=0,
        names=["gene_id", "transcript_id", "5_UTR", "CDS", "3_UTR", "5_len", "CDS_len", "3_len"]
    )

    # RNAfold needed T to be U so convert back for mapping
    for region in ["5_UTR", "CDS", "3_UTR"]:
        df[region] = df[region].map(convert_u_to_t)

    region_maps = parse_rnafold_root(Path(args.rnafold))

    for region in ["5_UTR", "CDS", "3_UTR"]:
        mfe_map = region_maps[region]
        df[f"{region}_MFE"] = df["transcript_id"].map(mfe_map)
        print(f"  → {df[f'{region}_MFE'].notna().sum()} MFE values mapped for {region}")

    output_path = Path(args.output).resolve()
    df.to_csv(output_path, index=False)
    print(f"\n[DONE] Successfully saved results to: {output_path}")

if __name__ == "__main__":
    main()