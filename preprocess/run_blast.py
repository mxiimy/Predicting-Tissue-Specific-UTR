import pandas as pd
import os
import sys
import argparse

INPUT_CSV = 'data/blast.csv'
FASTA_FILE = 'queries.fasta'
BLAST_OUT = 'results.tab'

def prepare_data():
    if not os.path.exists(INPUT_CSV):
        print(f"Error: {INPUT_CSV} not found.")
        sys.exit(1)
        
    print(f"Reading {INPUT_CSV}...")
    df = pd.read_csv(INPUT_CSV)
    
    with open(FASTA_FILE, "w") as f:
        for idx, row in df.iterrows():
            f.write(f">seq_{idx}\n{row['AA Seq']}\n")
    print(f"Created {FASTA_FILE} for BLASTing.")

def process_results():
    if not os.path.exists(BLAST_OUT):
        print("Error: results.tab not found. Did blastp run successfully?")
        sys.exit(1)

    df = pd.read_csv(INPUT_CSV)
    
    results = pd.read_csv(BLAST_OUT, sep='\t', 
                         names=['qseqid', 'pident', 'length', 'qlen'])

    results['is_perfect'] = (results['pident'] == 100.0) & (results['length'] == results['qlen'])
    results['original_index'] = results['qseqid'].str.replace('seq_', '').astype(int)
    
    df['100% Match'] = False
    perfect_indices = results[results['is_perfect']]['original_index']
    df.loc[perfect_indices, '100% Match'] = True

    df.to_csv('blast_finished.csv', index=False)
    print(f"Finished! Found {df['100% Match'].sum()} perfect matches.")
    print("Output saved to: blast_finished.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--process', action='store_true')
    args = parser.parse_args()

    if args.prepare:
        prepare_data()
    elif args.process:
        process_results()
    else:
        print("Specify --prepare or --process")
        sys.exit(1)