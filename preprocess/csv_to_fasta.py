import argparse
import pandas as pd
# Use LaTeX only for technical math; prose remains standard Markdown.
from pathlib import Path

def write_fasta(series_seq, series_id, path: Path):
    skipped = 0
    with open(path, "w") as fh:
        for tid, seq in zip(series_id, series_seq):
            if not seq or not isinstance(seq, str) or seq.strip() == "":
                skipped += 1
                continue
            seq_rna = seq.strip().upper().replace("T", "U")
            fh.write(f">{tid}\n{seq_rna}\n")
    print(f"  wrote {path}  ({len(series_seq) - skipped} sequences, {skipped} skipped)")

def main():
    parser = argparse.ArgumentParser(description="CSV → FASTA for RNAfold")
    parser.add_argument("--input",  required=True,  help="Path to input CSV")
    parser.add_argument("--outdir", required=True,  help="Directory for output FASTA files")
    parser.add_argument("--nrows",  type=int, default=None, help="Only process first N rows")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Use the first 5 columns and provide clean names to bypass the ' in headers
    print(f"[INFO] Loading {args.input}")
    df = pd.read_csv(
        args.input,
        header=0,
        usecols=[0, 1, 2, 3, 4],
        names=["gene_id", "transcript_id", "5_UTR_seq", "CDS_seq", "3_UTR_seq"],
        nrows=args.nrows,
    ).fillna("")

    # Map internal column names to output file names
    mapping = {
        "5_UTR_seq": "5_UTR",
        "CDS_seq": "CDS",
        "3_UTR_seq": "3_UTR"
    }

    for col, label in mapping.items():
        headers = df["transcript_id"] + "|" + label
        write_fasta(df[col], headers, outdir / f"{label}.fa")

    print("[DONE] FASTA files written to", outdir)

if __name__ == "__main__":
    main()