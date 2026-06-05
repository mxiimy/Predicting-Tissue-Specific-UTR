import sys
import os
import argparse  # Add this import
import pandas as pd
import numpy as np
from scipy.signal import convolve2d

# 1. Use argparse to handle the command line arguments
parser = argparse.ArgumentParser(description="Scan sequence for PWM")
parser.add_argument('--pwm-file', required=True, help="Path to the PWM file")
args = parser.parse_args()

# Now use args.pwm_file instead of sys.argv[1]
pwm_file = args.pwm_file
rbp_name = os.path.basename(pwm_file).replace('.txt', '')

# 2. Use absolute path to the data
df = pd.read_csv("/lustre09/project/6007512/HeDS/melody/Predicting-Tissue-Specific-UTR/data/finalDB.csv")

def load_pwm(filepath):
    df_pwm = pd.read_csv(filepath, sep='\t')
    # Ensure columns A, C, G, U exist
    matrix = df_pwm[['A', 'C', 'G', 'U']].values
    return matrix

def scan_seq_fast(sequence, pwm_matrix):
    if not isinstance(sequence, str) or len(sequence) < len(pwm_matrix):
        return 0.0
    mapping = {'A': 0, 'C': 1, 'G': 2, 'U': 3, 'T': 3}
    seq_ints = [mapping.get(n, -1) for n in sequence]
    if -1 in seq_ints: return 0.0
    
    seq_one_hot = np.eye(4)[seq_ints]
    # Perform convolution to find max score
    full_scores = convolve2d(seq_one_hot, pwm_matrix[::-1, :], mode='valid').sum(axis=1)
    return np.max(full_scores) if full_scores.size > 0 else 0.0

# 3. Process just this one RBP
print(f"Scanning for {rbp_name} using file {pwm_file}...")
pwm_matrix = load_pwm(pwm_file)

# We use direct assignment to avoid potential Pandas view issues
df[f'3_UTR_{rbp_name}_score'] = df['3\' UTR'].apply(lambda seq: scan_seq_fast(seq, pwm_matrix))
df[f'5_UTR_{rbp_name}_score'] = df['5\' UTR'].apply(lambda seq: scan_seq_fast(seq, pwm_matrix))

# 4. Save only the relevant columns to a temp file
output_cols = ['transcript_id', f'3_UTR_{rbp_name}_score', f'5_UTR_{rbp_name}_score']
df[output_cols].to_csv(f"temp_{rbp_name}.csv", index=False)

print(f"Finished {rbp_name}. Results saved to temp_{rbp_name}.csv")