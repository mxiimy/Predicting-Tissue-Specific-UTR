import argparse
import pandas as pd
from pathlib import Path
import RNA
import numpy as np

def calculate_mfe(sequence):
    if not sequence or not isinstance(sequence, str) or sequence.strip() == "":
        return float("nan")
    # RNAfold works with U, not T
    sequence_rna = sequence.strip().upper().replace("T", "U")
    fold_compound = RNA.fold_compound(sequence_rna)
    _, mfe = fold_compound.mfe()
    return mfe

def main():
    input_path = Path("./data/initialDB_sequences.csv")
    output_path = Path("./data/mfe_results.csv")
    
    # 1. Load with quotechar to handle the names properly
    df = pd.read_csv(input_path, quotechar='"')
    
    # Standardize headers just in case
    df.columns = df.columns.str.strip()
    print("Columns in file:", df.columns.tolist())
    
    # 2. Map correctly using the names you provided
    # Keys = Column name in CSV, Values = Column name for MFE result
    region_map = {
        "5' UTR": "5' UTR_MFE",
        "CDS": "CDS_MFE",
        "3' UTR": "3' UTR_MFE"
    }

    # Ensure MFE columns exist in the dataframe
    for mfe_col in region_map.values():
        if mfe_col not in df.columns:
            df[mfe_col] = np.nan

    # 3. Iterate through the map
    for col_name, mfe_col in region_map.items():
        # Identify rows where MFE is missing
        missing_mask = df[mfe_col].isna()
        missing_indices = df.index[missing_mask].tolist()
        
        if not missing_indices:
            print(f"[INFO] No missing MFE values for {col_name}.")
            continue
            
        print(f"[INFO] Calculating MFE for {len(missing_indices)} missing {col_name} sequences...")
        
        # Batch processing loop
        batch_size = 1000
        for i in range(0, len(missing_indices), batch_size):
            batch_indices = missing_indices[i : i + batch_size]
            
            for idx in batch_indices:
                seq = df.at[idx, col_name]
                df.at[idx, mfe_col] = calculate_mfe(seq)
            
            # Periodically save progress
            df.to_csv(output_path, index=False)
            print(f"    Processed up to index {batch_indices[-1]} for {col_name}. Saved progress.")

    print(f"\n[DONE] Successfully saved MFE results to: {output_path}")

if __name__ == "__main__":
    main()