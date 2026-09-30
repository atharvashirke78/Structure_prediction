import pandas as pd
import os

# Your output directory
data_dir = r"C:\IndiskaAI\Data\Processed"

# The exact names of your final CSVs
csv_files = [
    'validated_results_AB_Fold_3_Final.csv',
    'validated_results_Alpha_Fold_2_Final.csv',
    'validated_results_Boltz_2_Final.csv',
    'validated_results_ESM_2_Final.csv',
    'validated_results_Ibex_Final.csv',
    'validated_results_IgFold_Final.csv',
    'validated_results_Intellifold_Final.csv'
]

dataframes = {}
pdb_sets = []

print("Scanning datasets...")
for file in csv_files:
    filepath = os.path.join(data_dir, file)
    if os.path.exists(filepath):
        # Drop duplicates just in case (fixes the IgFold 111 issue)
        df = pd.read_csv(filepath).drop_duplicates(subset=['PDB_ID'])
        dataframes[file] = df
        
        # Extract the list of successfully processed PDB IDs
        # Only count rows where Masked_Fv_RMSD didn't crash
        valid_pdbs = df[df['Masked_Fv_RMSD'] != 'ERR_CRASH']['PDB_ID'].dropna().tolist()
        pdb_sets.append(set(valid_pdbs))
        print(f"{file}: {len(valid_pdbs)} valid targets")
    else:
        print(f"MISSING FILE: {file}")

# Find the intersection of ALL sets
if pdb_sets:
    common_pdbs = set.intersection(*pdb_sets)
    print("\n" + "="*50)
    print(f"STRICT INTERSECTION: {len(common_pdbs)} targets solved by ALL models.")
    print("="*50 + "\n")

    # Filter and save the new CSVs
    for file, df in dataframes.items():
        intersected_df = df[df['PDB_ID'].isin(common_pdbs)]
        
        new_filename = file.replace("_Final.csv", "_Intersected.csv")
        new_filepath = os.path.join(data_dir, new_filename)
        
        intersected_df.to_csv(new_filepath, index=False)
        print(f"Saved {new_filename} ({len(intersected_df)} rows)")
else:
    print("No data found to intersect.")