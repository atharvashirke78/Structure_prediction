import os
import csv
import subprocess
import glob
import re
from pymol import cmd

def run_universal_pipeline():
    # --- 1. CONFIGURE YOUR RUN HERE ---
    model_name = "AB Fold 3" # Change this for Boltz 2, Alpha Fold 2, etc.
    
    # Paths dynamically map to your new structure
    model_dir = rf"C:\IndiskaAI\Models\{model_name} outputs"
    tmalign_exe = r"C:\IndiskaAI\Tools\USalignWin64\USalign\USalign.exe"
    
    input_csv = os.path.join(model_dir, "data.csv")
    output_csv = os.path.join(model_dir, f"{model_name.replace(' ', '_')}_validated_results.csv")
    temp_ref = r"C:\IndiskaAI\Pipelines\temp_ref.pdb"
    temp_pred = r"C:\IndiskaAI\Pipelines\temp_pred.pdb"

    cmd.set('suspend_undo', 1)
    
    # --- 2. DYNAMIC ERROR HANDLING: Determine data source ---
    tasks = [] 
    
    if os.path.exists(input_csv):
        print(f"Found {input_csv}. Reading targets from sheet...")
        with open(input_csv, 'r', encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                pdb_id = row.get('PDB_ID', '').strip()
                if pdb_id and pdb_id != 'PDB_ID':
                    matches = glob.glob(os.path.join(model_dir, f"*{pdb_id}*.pdb"))
                    if matches:
                        tasks.append({
                            'pdb_id': pdb_id,
                            'filepath': matches[0],
                            'sheet_tm': row.get('TM', 'N/A').strip(),
                            'sheet_h3': row.get('CDR_H3_RMSD', 'N/A').strip()
                        })
    else:
        print(f"WARNING: No data.csv found in {model_dir}!")
        print("FALLBACK TRIGGERED: Scanning folder for PDB files automatically...")
        for filepath in glob.glob(os.path.join(model_dir, "*.pdb")):
            filename = os.path.basename(filepath)
            # Extracts the 4-character PDB ID (e.g., '3QWO' from 'abfold_3QWO_final.pdb')
            match = re.search(r'[1-9][a-zA-Z0-9]{3}', filename)
            if match:
                tasks.append({
                    'pdb_id': match.group(0).upper(),
                    'filepath': filepath,
                    'sheet_tm': 'N/A',
                    'sheet_h3': 'N/A'
                })

    if not tasks:
        print("ERROR: No valid PDB files found to process. Check your folder path!")
        return

    print("\n" + "="*130)
    print(f"{'PDB_ID':<8} | {'Sheet TM':<10} | {'True Fv TM':<15} | {'PyMOL Fv RMSD':<20} | {'Sheet H3':<10} | {'PyMOL H3 RMSD':<15} | {'Status':<15}")
    print("="*130)

    # --- 3. PROCESSING LOOP ---
    with open(output_csv, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['PDB_ID', 'Sheet_TM', 'True_Fv_TM', 'PyMOL_Fv_RMSD', 'Sheet_H3', 'PyMOL_H3_RMSD', 'Status'])
        
        for task in tasks:
            pdb_id = task['pdb_id']
            filepath = task['filepath']
            cmd.reinitialize()
            
            try:
                cmd.load(filepath, "pred_full")
                cmd.fetch(pdb_id, "ref_full", type="pdb")
                
                # Surgery: Extract Fv to make it fair
                cmd.create("ref_fv", "ref_full and polymer and resi 1-120")
                cmd.create("pred_fv", "pred_full and polymer and resi 1-120")
                cmd.delete("ref_full")
                cmd.delete("pred_full")
                
                # Global Fv RMSD
                global_align = cmd.align("pred_fv and name CA", "ref_fv and name CA", cycles=0)
                calc_global = f"{global_align[0]:.4f}"

                # Framework and H3 RMSD
                chains = cmd.get_chains("pred_fv")
                h_chain = 'H' if 'H' in chains else (chains[0] if len(chains) > 0 else 'A')
                cmd.super(f"pred_fv and chain {h_chain} and ss s", "ref_fv and ss s")
                
                if cmd.count_atoms("ref_fv and resi 90-120 and name CA") < 5:
                    calc_h3_str, status = "N/A", "Loop Unresolved"
                else:
                    h3_align = cmd.align(f"pred_fv and chain {h_chain} and resi 90-120 and name CA", 
                                         "ref_fv and resi 90-120 and name CA", transform=0, cycles=0)
                    calc_h3_str, status = f"{h3_align[0]:.4f}", "Success"

                # TM-Score via USalign
                true_tm = "N/A"
                if os.path.exists(tmalign_exe):
                    cmd.save(temp_ref, "ref_fv")
                    cmd.save(temp_pred, "pred_fv")
                    result = subprocess.run([tmalign_exe, temp_pred, temp_ref], capture_output=True, text=True)
                    for line in result.stdout.split('\n'):
                        if line.startswith("TM-score="):
                            true_tm = line.split()[1]
                            break
                            
                print(f"{pdb_id:<8} | {task['sheet_tm']:<10} | {true_tm:<15} | {calc_global:<20} | {task['sheet_h3']:<10} | {calc_h3_str:<15} | {status:<15}")
                writer.writerow([pdb_id, task['sheet_tm'], true_tm, calc_global, task['sheet_h3'], calc_h3_str, status])
                outfile.flush()
                
            except Exception as e:
                print(f"{pdb_id:<8} | ERROR: Failed to align")
                writer.writerow([pdb_id, task['sheet_tm'], "ERROR", "ERROR", task['sheet_h3'], "ERROR", "Failed"])
                outfile.flush()

    if os.path.exists(temp_ref): os.remove(temp_ref)
    if os.path.exists(temp_pred): os.remove(temp_pred)
    cmd.set('suspend_undo', 0)
    print("="*130)
    print(f"PIPELINE COMPLETE. Data saved to {output_csv}")

cmd.extend("run_universal_pipeline", run_universal_pipeline)