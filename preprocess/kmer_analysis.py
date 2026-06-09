import argparse
from collections import Counter
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "data" / "finalDB.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "full_kmer_analysis"


def append_csv(df: pd.DataFrame, path: Path):
    if df.empty:
        return
    write_header = not path.exists()
    df.to_csv(path, index=False, mode="a", header=write_header)

def normalize_sequence(value) -> str:
    if pd.isna(value):
        return ""
    sequence = str(value).strip().upper().replace("T", "U")
    return sequence


def pick_column(df: pd.DataFrame, candidates: list[str]) -> str:
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    raise KeyError(f"None of these columns were found: {', '.join(candidates)}")


def iter_kmers(sequence: str, k: int):
    limit = len(sequence) - k + 1
    for index in range(limit):
        yield sequence[index : index + k]


def summarize_region(sequences: pd.Series, region: str, min_k: int, max_k: int, top_n: int):
    summary_rows = []
    top_rows = []
    full_count_rows = []

    cleaned_sequences = [normalize_sequence(value) for value in sequences.tolist()]
    cleaned_sequences = [sequence for sequence in cleaned_sequences if sequence]

    print(f"[INFO] {region}: {len(cleaned_sequences)} non-empty sequences")

    for k in range(min_k, max_k + 1):
        counter = Counter()
        sequences_with_windows = 0
        total_windows = 0

        for sequence in cleaned_sequences:
            if len(sequence) < k:
                continue

            sequences_with_windows += 1
            window_count = len(sequence) - k + 1
            total_windows += window_count
            counter.update(iter_kmers(sequence, k))

        unique_kmers = len(counter)
        if counter:
            top_kmer, top_count = counter.most_common(1)[0]
        else:
            top_kmer, top_count = "", 0

        summary_rows.append(
            {
                "region": region,
                "k": k,
                "sequences_with_windows": sequences_with_windows,
                "total_windows": total_windows,
                "unique_kmers": unique_kmers,
                "top_kmer": top_kmer,
                "top_kmer_count": top_count,
            }
        )

        for rank, (kmer, count) in enumerate(counter.most_common(top_n), start=1):
            top_rows.append(
                {
                    "region": region,
                    "k": k,
                    "rank": rank,
                    "kmer": kmer,
                    "count": count,
                }
            )

        for kmer, count in counter.items():
            full_count_rows.append(
                {
                    "region": region,
                    "k": k,
                    "kmer": kmer,
                    "count": count,
                }
            )

        print(
            f"[INFO] {region} k={k}: {unique_kmers} unique kmers, "
            f"top={top_kmer or 'NA'} ({top_count})"
        )

    return summary_rows, top_rows, full_count_rows


def checkpoint_name(prefix: str, chunk_index: int) -> str:
    return f"{prefix}_chunk_{chunk_index:05d}.csv"


def chunk_index_from_path(path: Path) -> int:
    return int(path.stem.rsplit("_chunk_", 1)[-1])


def sorted_checkpoint_paths(checkpoint_dir: Path, prefix: str) -> list[Path]:
    paths = list(checkpoint_dir.glob(f"{prefix}_chunk_*.csv"))
    return sorted(paths, key=chunk_index_from_path)


def main():
    parser = argparse.ArgumentParser(
        description="Run k-mer analysis for 5' and 3' UTRs across a k range."
    )
    parser.add_argument(
        "--input",
        default=str(DEFAULT_INPUT),
        help="Path to the input CSV (default: data/finalDB.csv)",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Directory where output tables will be written",
    )
    parser.add_argument("--min-k", type=int, default=3, help="Smallest k to analyze")
    parser.add_argument("--max-k", type=int, default=494, help="Largest k to analyze")
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top kmers to save per k and region",
    )
    parser.add_argument(
        "--write-full-counts",
        action="store_true",
        help="Write every observed kmer count to disk. This can be large.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1000,
        help="Number of sequences to process per checkpoint chunk",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from existing checkpoint chunks instead of starting over",
    )
    parser.add_argument(
        "--nrows",
        type=int,
        default=None,
        help="Only process the first N rows, which is useful for testing",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Loading {input_path}")
    df = pd.read_csv(input_path, nrows=args.nrows)
    print(f"[INFO] {len(df)} rows loaded")

    five_utr_col = pick_column(df, ["5_UTR", "5' UTR"])
    three_utr_col = pick_column(df, ["3_UTR", "3' UTR"])

    summary_path = output_dir / "kmer_summary.csv" 
    top_path = output_dir / "kmer_top_hits.csv"
    counts_path = output_dir / "kmer_counts.csv"

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    if not args.resume:
        for path in [summary_path, top_path, counts_path]:
            if path.exists():
                path.unlink()
        for checkpoint_path in checkpoint_dir.glob("*.csv"):
            checkpoint_path.unlink()
    else:
        summary_checkpoints = sorted_checkpoint_paths(checkpoint_dir, "kmer_summary")
        top_checkpoints = sorted_checkpoint_paths(checkpoint_dir, "kmer_top_hits")
        count_checkpoints = sorted_checkpoint_paths(checkpoint_dir, "kmer_counts")

        if summary_checkpoints:
            print(f"[INFO] Resuming from {len(summary_checkpoints)} checkpoint chunk(s)")
            if summary_path.exists():
                summary_path.unlink()
            if top_path.exists():
                top_path.unlink()
            if counts_path.exists():
                counts_path.unlink()

            for checkpoint_path in summary_checkpoints:
                append_csv(pd.read_csv(checkpoint_path), summary_path)
            for checkpoint_path in top_checkpoints:
                append_csv(pd.read_csv(checkpoint_path), top_path)
            if args.write_full_counts:
                for checkpoint_path in count_checkpoints:
                    append_csv(pd.read_csv(checkpoint_path), counts_path)

            completed_chunks = [chunk_index_from_path(path) for path in summary_checkpoints]
            last_completed_chunk = max(completed_chunks)
        else:
            print("[INFO] --resume requested, but no checkpoints were found; starting fresh")
            last_completed_chunk = 0

    if not args.resume:
        last_completed_chunk = 0

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be a positive integer")

    chunk_index = 0
    n_rows = len(df)

    for start in range(0, n_rows, args.chunk_size):
        chunk_index += 1
        if chunk_index <= last_completed_chunk:
            continue
        end = min(start + args.chunk_size, n_rows)
        print(f"[INFO] Processing rows {start + 1}-{end} (chunk {chunk_index})")

        summary_rows = []
        top_rows = []
        full_count_rows = []

        for region, column_name in [("5_UTR", five_utr_col), ("3_UTR", three_utr_col)]:
            region_summary, region_top, region_full = summarize_region(
                df[column_name].iloc[start:end],
                region=region,
                min_k=args.min_k,
                max_k=args.max_k,
                top_n=args.top_n,
            )
            summary_rows.extend(region_summary)
            top_rows.extend(region_top)
            if args.write_full_counts:
                full_count_rows.extend(region_full)

        chunk_summary_df = pd.DataFrame(summary_rows)
        chunk_top_df = pd.DataFrame(top_rows)

        chunk_summary_path = checkpoint_dir / checkpoint_name("kmer_summary", chunk_index)
        chunk_top_path = checkpoint_dir / checkpoint_name("kmer_top_hits", chunk_index)

        chunk_summary_df.to_csv(chunk_summary_path, index=False)
        chunk_top_df.to_csv(chunk_top_path, index=False)

        append_csv(chunk_summary_df, summary_path)
        append_csv(chunk_top_df, top_path)

        if args.write_full_counts:
            chunk_counts_df = pd.DataFrame(full_count_rows)
            chunk_counts_path = checkpoint_dir / checkpoint_name("kmer_counts", chunk_index)
            chunk_counts_df.to_csv(chunk_counts_path, index=False)
            append_csv(chunk_counts_df, counts_path)

        print(f"[INFO] Saved checkpoint files for chunk {chunk_index}")

    print(f"[INFO] Saved summary -> {summary_path}")
    print(f"[INFO] Saved top hits -> {top_path}")

    if args.write_full_counts:
        print(f"[INFO] Saved full counts -> {counts_path}")

    print("[DONE] k-mer analysis complete")


if __name__ == "__main__":
    main()
