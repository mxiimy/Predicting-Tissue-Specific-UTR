import argparse
import os
import numpy as np
import pandas as pd
import RNA                  
from multiprocessing import Pool
from pathlib import Path


# ── helpers ──────────────────────────────────────────────────────────────────

def calculate_mfe(seq: str):
    """Return (structure, mfe_kcal_mol) for a single RNA sequence string.
    Returns (None, NaN) for empty / non-string input."""
    if not seq or not isinstance(seq, str) or seq.strip() == "" or seq.upper() == "NAN":
        return None, np.nan
    seq_rna = seq.strip().upper().replace("T", "U")
    fc = RNA.fold_compound(seq_rna)
    ss, mfe = fc.mfe()
    return ss, mfe


def mfe_only(seq: str) -> float:
    """Wrapper returning only the MFE float (for use with multiprocessing.Pool)."""
    return calculate_mfe(seq)[1]


def make_junction_seq(utr5: str, cds: str, utr5_tail: int = 30, cds_head: int = 10) -> str:
    """Return the last `utr5_tail` nt of 5'UTR + first `cds_head` nt of CDS.
    Returns empty string if both parts are missing."""
    utr5 = utr5 if isinstance(utr5, str) else ""
    cds  = cds  if isinstance(cds,  str) else ""
    if utr5.upper() == "NAN": utr5 = ""
    if cds.upper() == "NAN": cds = ""
    return utr5[-utr5_tail:] + cds[:cds_head]


def parallel_mfe(series: pd.Series, n_cpus: int) -> pd.Series:
    """Apply mfe_only across a pandas Series using a multiprocessing Pool."""
    seqs = series.fillna("").tolist()
    with Pool(processes=n_cpus) as pool:
        results = pool.map(mfe_only, seqs)
    return pd.Series(results, index=series.index, dtype=float)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute MFE for 5'UTR–CDS junction only if missing")
    parser.add_argument("--input",  required=True,  help="Path to input CSV")
    parser.add_argument("--output", required=True,  help="Path to output CSV")
    parser.add_argument("--cpus",   type=int, default=1,
                        help="Number of CPUs for parallel folding (match --cpus-per-task)")
    parser.add_argument("--nrows",  type=int, default=None,
                        help="Only process first N rows (useful for testing)")
    args = parser.parse_args()

    # 1. Check for pre-calculated junction scores in master database
    master_path = Path("../data/finalDB.csv")
    existing_data = None
    
    if master_path.exists():
        print(f"[INFO] Found existing database at {master_path}. Checking for pre-calculated junction MFE...")
        try:
            existing_data = pd.read_csv(master_path, usecols=['transcript_id', '5UTR_CDS_junction_MFE'])
            existing_data['transcript_id'] = existing_data['transcript_id'].astype(str)
        except Exception as e:
            print(f"[WARNING] Could not read existing database columns ({e}). Running complete recalculation.")
            existing_data = None

    # ── load input data ───────────────────────────────────────────────────────
    print(f"[INFO] Loading input file: {args.input}")
    df = pd.read_csv(args.input, nrows=args.nrows)
    
    # Standardize expected columns dynamically to avoid structural mismatches
    rename_dict = {}
    col_5utr = None
    col_cds = None
    
    for col in df.columns:
        c_low = col.lower()
        if col.lower() in ['gene_id', 'geneid']: 
            rename_dict[col] = 'gene_id'
        elif col.lower() in ['transcript_id', 'transcriptid']: 
            rename_dict[col] = 'transcript_id'
        elif '5' in c_low and 'utr' in c_low and 'mfe' not in c_low and 'len' not in c_low:
            col_5utr = col
            rename_dict[col] = '5_UTR'
        elif c_low == 'cds':
            col_cds = col
            rename_dict[col] = 'CDS'
    
    # Fallback to positional values if naming structures vary completely
    if not col_5utr: rename_dict[df.columns[2]] = '5_UTR'
    if not col_cds:  rename_dict[df.columns[3]] = 'CDS'
    
    df = df.rename(columns=rename_dict)
    df['transcript_id'] = df['transcript_id'].astype(str)
    print(f"[INFO] {len(df)} rows loaded.")

    # ── merge existing scores if found ────────────────────────────────────────
    if existing_data is not None:
        if '5UTR_CDS_junction_MFE' in df.columns:
            df = df.drop(columns=['5UTR_CDS_junction_MFE'])
        df = df.merge(existing_data, on='transcript_id', how='left')

    # Ensure sequence data columns are clean strings
    df['5_UTR'] = df['5_UTR'].fillna("").astype(str).replace(["nan", "NaN", "NAN"], "")
    df['CDS'] = df['CDS'].fillna("").astype(str).replace(["nan", "NaN", "NAN"], "")

    # ── MFE: 5'UTR–CDS junction (last 30 nt of 5'UTR + first 10 nt of CDS) ───
    junction_mfe_col = "5UTR_CDS_junction_MFE"
    if junction_mfe_col not in df.columns:
        df[junction_mfe_col] = np.nan

    print("[INFO] Constructing junction sequences...")
    # Safe 1D slicing isolation handles duplicate columns without returning DataFrames
    utr5_array = df[["5_UTR"]].iloc[:, 0].astype(str).tolist()
    cds_array = df[["CDS"]].iloc[:, 0].astype(str).tolist()

    junction_seqs = pd.Series(
        [make_junction_seq(u, c) for u, c in zip(utr5_array, cds_array)],
        index=df.index,
    ).fillna("").astype(str)
    
    df["5UTR_CDS_junction_length"] = junction_seqs.str.len()
    
    # Determine exact calculation masking criteria
    needs_junction_mask = (junction_seqs != "") & (junction_seqs != "nan") & (df[junction_mfe_col].isna())
    skipped_junctions = len(df) - needs_junction_mask.sum()
    
    if needs_junction_mask.sum() == 0:
        print(f"[INFO] Skipping 5UTR_CDS_junction: All {skipped_junctions} sequences already calculated.")
    else:
        print(f"[INFO] Calculating MFE for 5UTR_CDS_junction ({needs_junction_mask.sum()} remaining, skipping {skipped_junctions})...")
        calculated_junction_scores = parallel_mfe(junction_seqs[needs_junction_mask], n_cpus=args.cpus)
        df.loc[needs_junction_mask, junction_mfe_col] = calculated_junction_scores
        
    non_null_junc = df[junction_mfe_col].notna().sum()
    print(f"[INFO]   → Total {non_null_junc}/{len(df)} junction sequences filled successfully")

    # ── save ──────────────────────────────────────────────────────────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"[INFO] Saved → {args.output}")


if __name__ == "__main__":
    main()