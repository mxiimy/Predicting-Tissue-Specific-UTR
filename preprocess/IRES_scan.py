import pandas as pd

ire_df = pd.read_csv('data/ires_library.csv')
seq_df = pd.read_csv('data/initialDB_sequences.csv')

def scan_ires(utr_sequence):
    # Ensure sequence is a string and check minimum length constraint
    if pd.isna(utr_sequence) or len(str(utr_sequence)) < 28:
        return 0, ""
    
    found_ids = []
    utr_str = str(utr_sequence).upper()
    
    # Iterate through the library to find matches
    for _, row in ire_df.iterrows():
        ire_seq = str(row['IRES sequence']).upper()
        ire_id = str(row['IRES ID'])
        
        if ire_seq in utr_str:
            found_ids.append(ire_id)
            
    return len(found_ids), "; ".join(found_ids)

seq_df[['IRES count', 'IRES IDs']] = seq_df['5\' UTR'].apply(
    lambda x: pd.Series(scan_ires(x))
)

seq_df.to_csv('data/utr_ire_scan_results.csv', index=False)

print("Scan complete. Results saved to 'data/utr_ire_scan_results.csv'.")