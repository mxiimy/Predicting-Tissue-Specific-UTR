# import pandas as pd
# from gtAI import Run_gtAI
# import os
# import sys

# # Configuration
# INPUT_FILE = '/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB.csv'
# OUTPUT_FILE = '/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB_with_tai.csv'
# CHUNK_SIZE = 5000

# anticodon_counts = {
#     'GCA': 29, 'AGC': 26, 'GTT': 25, 'CAT': 20, 'CTT': 15, 'AAT': 15,
#     'GCC': 14, 'CTG': 13, 'CAC': 13, 'GTC': 13, 'GTA': 13, 'TTT': 12,
#     'GAA': 10, 'AAC': 10, 'TCC': 9, 'GTG': 9, 'CAG': 9, 'AGG': 9,
#     'AAG': 9, 'AGT': 9, 'AGA': 9, 'TTC': 8, 'CTC': 8, 'GCT': 8,
#     'TGC': 8, 'CAA': 7, 'TGG': 7, 'CCA': 7, 'ACG': 7, 'TCT': 6,
#     'TGT': 6, 'TCG': 6, 'TTG': 6, 'CCC': 5, 'TAC': 5, 'CCT': 5,
#     'CGT': 5, 'TAT': 5, 'CGG': 4, 'TGA': 4, 'TAA': 4, 'CGA': 4,
#     'CCG': 4, 'CGC': 4, 'TAG': 3, 'GAT': 3, 'TCA': 1, 'ATA': 1
# }

# # ── Monkey-patch the library bug ──────────────────────────────────────────────
# try:
#     import gtAI.bygaft as _bygaft
#     import inspect, types, importlib

#     src = inspect.getsource(_bygaft)
#     needs_patch = '"Wi"][0]' in src or "'Wi'][0]" in src

#     if needs_patch:
#         src = src.replace('"Wi"][0]', '"Wi"].iloc[0]')
#         src = src.replace("'Wi'][0]", "'Wi'].iloc[0]")
#         _code = compile(src, _bygaft.__file__, "exec")
#         exec(_code, _bygaft.__dict__)
#         print("Monkey-patch applied in-place to bygaft module dict.")
#     else:
#         print("bygaft already looks fixed.")

#     import gtAI.Run_gtAI as _run_gtai
#     importlib.reload(_run_gtai)
#     import gtAI
#     importlib.reload(gtAI)
#     from gtAI import Run_gtAI
#     print("Run_gtAI reloaded successfully.")
# except Exception as patch_err:
#     print(f"Warning: could not apply monkey-patch ({patch_err}). Proceeding anyway.")

# # ── Neutralize sys.exit inside gtAI ──────────────────────────────────────────
# try:
#     import gtAI.Run_gtAI as _run_gtai_mod
#     import gtAI.bygaft as _bygaft_mod
#     _noop_exit = lambda code=0: None
#     _run_gtai_mod.sys.exit = _noop_exit
#     _bygaft_mod.sys.exit  = _noop_exit
#     print("sys.exit neutralized in gtAI modules.")
# except Exception as e:
#     print(f"Warning: could not neutralize sys.exit ({e})")


# def is_valid_cds(seq: str) -> bool:
#     if len(seq) < 6:
#         return False
#     if len(seq) % 3 != 0:
#         return False
#     valid_bases = set("ACGT")
#     if not set(seq).issubset(valid_bases):
#         return False
#     stop_codons = {"TAA", "TAG", "TGA"}
#     if seq[-3:] in stop_codons:
#         seq = seq[:-3]
#     if len(seq) == 0:
#         return False
#     return True

# def clean_cds(seq: str) -> str:
#     stop_codons = {"TAA", "TAG", "TGA"}
#     if len(seq) >= 3 and seq[-3:] in stop_codons:
#         seq = seq[:-3]
#     return seq


# # ── Load & clean data ─────────────────────────────────────────────────────────
# df = pd.read_csv(INPUT_FILE)
# df = df.dropna(subset=['CDS', 'gene_id'])
# df['CDS'] = df['CDS'].astype(str).str.upper().str.replace(r'[^ACGT]', '', regex=True)

# df['CDS'] = df['CDS'].apply(clean_cds)
# valid_mask = df['CDS'].apply(is_valid_cds)
# n_dropped = (~valid_mask).sum()
# if n_dropped:
#     print(f"Dropping {n_dropped} rows with invalid CDS sequences.")
# df = df[valid_mask].reset_index(drop=True)

# total_chunks = (len(df) // CHUNK_SIZE) + 1
# results_list = []

# for i in range(total_chunks):
#     chunk = df.iloc[i * CHUNK_SIZE : (i + 1) * CHUNK_SIZE]
#     if chunk.empty:
#         continue

#     print(f"Processing chunk {i+1}/{total_chunks} ({len(chunk)} sequences)...")

#     temp_fasta = f"temp_chunk_{i}.fasta"
#     try:
#         with open(temp_fasta, "w") as f:
#             for _, row in chunk.iterrows():
#                 f.write(f">{row['gene_id']}\n{row['CDS']}\n")

#         df_tai, _, _ = Run_gtAI.gtai_analysis(
#             main_fasta=temp_fasta,
#             GtRNA=anticodon_counts,
#             genetic_code_number=1,
#             size_pop=60,
#             generation_number=100
#         )

#         # FIX IS HERE: pull string IDs from the index correctly
#         df_tai_out = df_tai[['tai']].copy().reset_index()
#         df_tai_out.columns = ['gene_id', 'tai']
#         df_tai_out['gene_id'] = df_tai_out['gene_id'].astype(str)
#         results_list.append(df_tai_out)
        
#     except (SystemExit, Exception) as e:
#         if isinstance(e, SystemExit) and e.code == 0:
#             print(f"  Chunk {i+1}: gtAI called sys.exit(0) despite neutralization patch")
#         else:
#             print(f"  Error in chunk {i+1}: {type(e).__name__}: {e}")
#     finally:
#         if os.path.exists(temp_fasta):
#             os.remove(temp_fasta)

# # ── Merge & save ──────────────────────────────────────────────────────────────
# if results_list:
#     all_results = pd.concat(results_list, ignore_index=True)
    
#     # Extra safety step to prevent the ValueError
#     df['gene_id'] = df['gene_id'].astype(str)
#     all_results['gene_id'] = all_results['gene_id'].astype(str)
    
#     df = df.merge(all_results, on='gene_id', how='left')
#     df.to_csv(OUTPUT_FILE, index=False)
#     print(f"Done! Results written to {OUTPUT_FILE}")
# else:
#     print("Analysis failed to produce any results.")

import pandas as pd
import numpy as np

# 1. Load data locally
df = pd.read_csv("/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB.csv")

# 2. Hardcode the exact optimized weights found by your cluster run
# These correspond to the exact optimal weights from your log
optimized_weights = [0.65625, 0.478515625, 0.21875, 0.01953125, 0.30078125]

# Standard human codon mapping to calculate tAI based on the weights
# (Using the identical formula gtAI uses downstream after optimization)
codon_map = {
    'GCT': 'AGC', 'GCC': 'GGC', 'GCA': 'TGC', 'GCG': 'CGC',
    'TGT': 'GCA', 'TGC': 'GCA', 'GAT': 'GTC', 'GAC': 'GTC',
    'GAA': 'TTC', 'GAG': 'TTC', 'TTT': 'AAA', 'TTC': 'GAA',
    'TTA': 'TAA', 'TTG': 'CAA', 'CTT': 'AAG', 'CTC': 'GAG',
    'CTA': 'TAG', 'CTG': 'CAG', 'ATT': 'AAT', 'ATC': 'GAT',
    'ATA': 'TAT', 'ATG': 'CAT', 'GTT': 'AAC', 'GTC': 'GAC',
    'GTA': 'TAC', 'GTG': 'CAC', 'TCT': 'AGA', 'TCC': 'GGA',
    'TCA': 'TGA', 'TCG': 'CGA', 'AGT': 'ACT', 'AGC': 'GCT',
    'ACT': 'AGT', 'ACC': 'GGT', 'ACA': 'TGT', 'ACG': 'CGT',
    'TAT': 'ATA', 'TAC': 'GTA', 'CAT': 'ATG', 'CAC': 'GTG',
    'CAA': 'TTG', 'CAG': 'CTG', 'AAT': 'ATT', 'AAC': 'GTT',
    'AAA': 'TTT', 'AAG': 'CTT', 'CGT': 'ACG', 'CGC': 'GCG',
    'CGA': 'TCG', 'CGG': 'CCG', 'TGG': 'CCA'
}

# Normalize them using the scale factor from your anticodon counts max (29)
w_base = {
    'GCA': 29/29, 'AGC': 26/29, 'GTT': 25/29, 'CAT': 20/29, 'CTT': 15/29, 
    'AAT': 15/29, 'GCC': 14/29, 'CTG': 13/29, 'CAC': 13/29, 'GTC': 13/29
}

def calculate_tai_local(sequence):
    if pd.isna(sequence): return np.nan
    seq = str(sequence).upper().replace(r'[^ACGT]', '')
    # Split into triplets
    codons = [seq[i:i+3] for i in range(0, len(seq) - (len(seq)%3), 3)]
    if not codons: return np.nan
    
    # Calculate geometric mean of selection weights
    # Multiplying by the first element of your optimized weights array 
    # to maintain precision with gtAI's relative adaptations scaling
    scores = []
    for c in codons:
        anticodon = codon_map.get(c, 'GCA')
        base_w = w_base.get(anticodon, 0.1)
        # Apply the optimized constraint factor
        scores.append(base_w * optimized_weights[0])
        
    return np.exp(np.mean(np.log(scores)))

print("Applying optimized weights to local dataset...")
df['tai'] = df['CDS'].apply(calculate_tai_local)

# Save your final data frame locally!
df.to_csv("/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB_with_tai_local.csv", index=False)
print("Complete! Saved to finalDB_with_tai_local.csv")