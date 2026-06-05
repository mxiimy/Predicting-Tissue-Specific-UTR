import pandas as pd
from gtAI import Run_gtAI
import os
import sys

# Configuration
INPUT_FILE = '/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB.csv'
OUTPUT_FILE = '/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB_with_tai.csv'
CHUNK_SIZE = 5000

anticodon_counts = {
    'GCA': 29, 'AGC': 26, 'GTT': 25, 'CAT': 20, 'CTT': 15, 'AAT': 15,
    'GCC': 14, 'CTG': 13, 'CAC': 13, 'GTC': 13, 'GTA': 13, 'TTT': 12,
    'GAA': 10, 'AAC': 10, 'TCC': 9, 'GTG': 9, 'CAG': 9, 'AGG': 9,
    'AAG': 9, 'AGT': 9, 'AGA': 9, 'TTC': 8, 'CTC': 8, 'GCT': 8,
    'TGC': 8, 'CAA': 7, 'TGG': 7, 'CCA': 7, 'ACG': 7, 'TCT': 6,
    'TGT': 6, 'TCG': 6, 'TTG': 6, 'CCC': 5, 'TAC': 5, 'CCT': 5,
    'CGT': 5, 'TAT': 5, 'CGG': 4, 'TGA': 4, 'TAA': 4, 'CGA': 4,
    'CCG': 4, 'CGC': 4, 'TAG': 3, 'GAT': 3, 'TCA': 1, 'ATA': 1
}

# ── Monkey-patch the library bug ──────────────────────────────────────────────
# bygaft.py does `RSCU_wai_corr["Wi"][0]` which fails when the Series index
# is not integer-positional. `.iloc[0]` is the correct fix.
try:
    import gtAI.bygaft as _bygaft
    import inspect, types

    src = inspect.getsource(_bygaft)
    if '"Wi"][0]' in src or "'Wi'][0]" in src:
        src = src.replace('"Wi"][0]', '"Wi"].iloc[0]')
        src = src.replace("'Wi'][0]", "'Wi'].iloc[0]")
        _code = compile(src, _bygaft.__file__, "exec")
        _new_mod = types.ModuleType(_bygaft.__name__)
        exec(_code, _new_mod.__dict__)
        import sys as _sys
        _sys.modules[_bygaft.__name__] = _new_mod
        # Re-import Run_gtAI so it picks up the patched bygaft
        import importlib
        import gtAI.Run_gtAI as _run_gtai
        importlib.reload(_run_gtai)
        from gtAI import Run_gtAI
        print("Monkey-patch applied successfully.")
    else:
        print("Library source looks already fixed; no patch needed.")
except Exception as patch_err:
    print(f"Warning: could not apply monkey-patch ({patch_err}). Proceeding anyway.")
# ─────────────────────────────────────────────────────────────────────────────

def is_valid_cds(seq: str) -> bool:
    """Return True only for sequences gtAI can process without crashing."""
    if len(seq) < 6:           # need at least 2 codons
        return False
    if len(seq) % 3 != 0:     # must be codon-aligned
        return False
    valid_bases = set("ACGT")
    if not set(seq).issubset(valid_bases):
        return False
    # Strip trailing stop codon if present (TAA, TAG, TGA)
    # gtAI typically expects the stop codon removed
    stop_codons = {"TAA", "TAG", "TGA"}
    if seq[-3:] in stop_codons:
        seq = seq[:-3]
    if len(seq) == 0:
        return False
    return True

def clean_cds(seq: str) -> str:
    """Strip trailing stop codon so gtAI doesn't choke on it."""
    stop_codons = {"TAA", "TAG", "TGA"}
    if len(seq) >= 3 and seq[-3:] in stop_codons:
        seq = seq[:-3]
    return seq


# ── Load & clean data ─────────────────────────────────────────────────────────
df = pd.read_csv(INPUT_FILE)
df = df.dropna(subset=['CDS', 'gene_id'])
df['CDS'] = df['CDS'].astype(str).str.upper().str.replace(r'[^ACGT]', '', regex=True)

# Apply validity filter BEFORE chunking
df['CDS'] = df['CDS'].apply(clean_cds)
valid_mask = df['CDS'].apply(is_valid_cds)
n_dropped = (~valid_mask).sum()
if n_dropped:
    print(f"Dropping {n_dropped} rows with invalid CDS sequences.")
df = df[valid_mask].reset_index(drop=True)

total_chunks = (len(df) // CHUNK_SIZE) + 1
results_list = []

for i in range(total_chunks):
    chunk = df.iloc[i * CHUNK_SIZE : (i + 1) * CHUNK_SIZE]
    if chunk.empty:
        continue

    print(f"Processing chunk {i+1}/{total_chunks} ({len(chunk)} sequences)...")

    temp_fasta = f"temp_chunk_{i}.fasta"
    try:
        with open(temp_fasta, "w") as f:
            for _, row in chunk.iterrows():
                f.write(f">{row['gene_id']}\n{row['CDS']}\n")

        df_tai, _, _ = Run_gtAI.gtai_analysis(
            main_fasta=temp_fasta,
            GtRNA=anticodon_counts,
            genetic_code_number=1,
            size_pop=60,
            generation_number=100
        )

        # Defensive: reset index so downstream merges are safe
        df_tai = df_tai.reset_index(drop=True)
        results_list.append(df_tai[['gene_id', 'tai']])

    except Exception as e:
        print(f"  Error in chunk {i+1}: {e}")
    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

# ── Merge & save ──────────────────────────────────────────────────────────────
if results_list:
    all_results = pd.concat(results_list, ignore_index=True)
    df = df.merge(all_results, on='gene_id', how='left')
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"Done! Results written to {OUTPUT_FILE}")
else:
    print("Analysis failed to produce any results.")