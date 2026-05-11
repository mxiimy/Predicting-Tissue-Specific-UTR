import argparse
import os
import numpy as np
import pandas as pd
import RNA                      # ViennaRNA Python bindings (module load viennarna)
from multiprocessing import Pool
from pathlib import Path


# ── helpers ──────────────────────────────────────────────────────────────────

def calculate_mfe(seq: str):
    """Return (structure, mfe_kcal_mol) for a single RNA sequence string.
    Returns (None, NaN) for empty / non-string input."""
    if not seq or not isinstance(seq, str) or seq.strip() == "":
        return None, np.nan
    # ViennaRNA expects RNA alphabet; convert T→U just in case input is DNA
    seq_rna = seq.strip().upper().replace("T", "U")
    fc = RNA.fold_compound(seq_rna)
    ss, mfe = fc.mfe()
    return ss, mfe


def mfe_only(seq: str) -> float:
    """Wrapper returning only the MFE float (for use with multiprocessing.Pool)."""
    return calculate_mfe(seq)[1]


def parallel_mfe(series: pd.Series, n_cpus: int) -> pd.Series:
    """Apply mfe_only across a pandas Series using a multiprocessing Pool."""
    seqs = series.fillna("").tolist()
    with Pool(processes=n_cpus) as pool:
        results = pool.map(mfe_only, seqs)
    return pd.Series(results, index=series.index, dtype=float)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Compute MFE for 5'UTR / CDS / 3'UTR")
    parser.add_argument("--input",  required=True,  help="Path to input CSV")
    parser.add_argument("--output", required=True,  help="Path to output CSV")
    parser.add_argument("--cpus",   type=int, default=1,
                        help="Number of CPUs for parallel folding (match --cpus-per-task)")
    parser.add_argument("--nrows",  type=int, default=None,
                        help="Only process first N rows (useful for testing)")
    args = parser.parse_args()

    # ── load ──────────────────────────────────────────────────────────────────
    print(f"[INFO] Loading {args.input}")
    df = pd.read_csv(
        args.input,
        sep=",",
        header=0,
        comment="#",
        names=["gene_id", "transcript_id", "5_UTR", "CDS", "3_UTR"],
        nrows=args.nrows,
    )
    print(f"[INFO] {len(df)} rows loaded")

    # ── lengths ───────────────────────────────────────────────────────────────
    for col in ["5_UTR", "CDS", "3_UTR"]:
        df[col] = df[col].fillna("")
        df[f"{col}_length"] = df[col].str.len()

    # ── MFE calculation ───────────────────────────────────────────────────────
    for col in ["5_UTR", "CDS", "3_UTR"]:
        print(f"[INFO] Calculating MFE for {col} using {args.cpus} CPU(s)...")
        df[f"{col}_MFE"] = parallel_mfe(df[col], n_cpus=args.cpus)
        non_null = df[f"{col}_MFE"].notna().sum()
        print(f"[INFO]   → {non_null}/{len(df)} sequences folded successfully")

    # ── save ──────────────────────────────────────────────────────────────────
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"[INFO] Saved → {args.output}")


if __name__ == "__main__":
    main()