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

df = pd.read_csv(INPUT_FILE)
df = df.dropna(subset=['CDS', 'gene_id'])
df['CDS'] = df['CDS'].astype(str).str.upper().str.replace(r'[^ACGT]', '', regex=True)
df = df[df['CDS'].str.len() > 0]
df = df.reset_index(drop=True)

total_chunks = (len(df) // CHUNK_SIZE) + 1
results_list = []

for i in range(total_chunks):
    chunk = df.iloc[i*CHUNK_SIZE : (i+1)*CHUNK_SIZE]
    if chunk.empty: continue
    
    print(f"Processing chunk {i+1}/{total_chunks}...")
    
    temp_fasta = f"temp_chunk_{i}.fasta"
    with open(temp_fasta, "w") as f:
        for _, row in chunk.iterrows():
            f.write(f">{row['gene_id']}\n{row['CDS']}\n")
            
    try:
        # Run analysis with explicit error handling
        df_tai, _, _ = Run_gtAI.gtai_analysis(
            main_fasta=temp_fasta, 
            GtRNA=anticodon_counts, 
            genetic_code_number=1,
            size_pop=60, 
            generation_number=100
        )
        results_list.append(df_tai[['gene_id', 'tai']])
    except Exception as e:
        print(f"Error in chunk {i}: {e}")
        # If gTAI fails, we print a warning and continue
    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

if results_list:
    all_results = pd.concat(results_list)
    df = df.merge(all_results, on='gene_id', how='left')
    df.to_csv(OUTPUT_FILE, index=False)
    print("Done!")
else:
    print("Analysis failed to produce results.")