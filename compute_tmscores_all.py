import os
import pymol
from tmalign import * 
import pandas as pd 
import re 

# Initialize PyMOL in headless mode (no GUI)
pymol.finish_launching(['pymol', '-cq'])  # '-cq' means quiet and no GUI

# Path to the reference file and the directory with mutated PDBs
cwd = os.getcwd()
sep = os.sep 
PDBS_DIR = cwd + sep + 'pdbs-alphafold'
RESULTS_DIR = cwd + sep + 'results'
SUBDIRS = ["monomer", "tetramer"]
exe_path = cwd + sep + "tmalign_exe" + sep + "TMalign_cpp"

# Function to compute TM-score for two structures
def compute_tm_score(mutant1_pdb, mutant2_pdb):
    if "wt" in mutant1_pdb:
        if "tetramer" in mutant1_pdb:
            mutant1 = "wt-tetramer"
        else:
            mutant1 = "wt"
    else:
        mutant1 = mutant1_pdb.split(".")[0].split("/")[-1].strip()
    if "wt" in mutant2_pdb:
        if "tetramer" in mutant1_pdb:
            mutant2 = "wt-tetramer"
        else:
            mutant2 = "wt"
    else:
        mutant2 = mutant2_pdb.split(".")[0].split("/")[-1].strip()
    # Align the structures and get the TM-score
    #print(exe_path)
    tm_score = tmscore(mutant1, mutant2, quiet = 1, ter = 1, exe=exe_path)
   
    return tm_score

# Function to extract numbers for sorting
def extract_number(filename):
    match = re.search(r'\d+', filename)  # Find the first number in the string
    return int(match.group()) if match else float('inf')  # Convert to int

# Function to compute TM-scores for all mutants in the pdbs directory
def compute_tm_scores(wt_pdb, pdbs_dir):
    
    files = os.listdir(pdbs_dir)
    # Sort list based on extracted number
    files_sorted = sorted(files, key=extract_number)
    if "wt" in files_sorted[-1]:
        del files_sorted[-1]

    if "tetramer" in wt_pdb:
        files_final = ["wt-tetramer.pdb"]
    else:    
        files_final = ["wt.pdb"]
    [files_final.append(f) for f in files_sorted]
    files = files_final 

    #Load all structures first
    for pdb_file in files:
        if "wt" in pdb_file:
            if "tetramer" in pdb_file:
                name = "wt-tetramer"
            else:
                name = "wt"
        else:
            name = pdb_file.split(".")[0].split("/")[-1].strip()
        pdb_path = pdbs_dir + sep + pdb_file
        cmd.load(pdb_path, name)

    tm_scores = {}
    # Loop through all files in the pdbs directory
    for pdb_file1 in files:
        
        tm_scores[pdb_file1.split(".")[0]]  = {}
        if pdb_file1.endswith('.pdb'):
            mutant1_pdb = os.path.join(pdbs_dir, pdb_file1)
        else:
            continue

        if "wt" in mutant1_pdb:
            mutant1 = "wt"
        else:
            mutant1 = mutant1_pdb.split(".")[0].split("/")[-1].strip()
            
        for pdb_file2 in files:
            if pdb_file2.endswith('.pdb'):
                mutant2_pdb = os.path.join(pdbs_dir, pdb_file2)   
            else:
                continue
            print(f"Computing TM-score between {pdb_file1} and {pdb_file2}...")
            tm_score = compute_tm_score(mutant1_pdb, mutant2_pdb)
              
            # Store the result in the dictionary
            tm_scores[pdb_file1.split(".")[0]][pdb_file2.split(".")[0]]  = tm_score
            
    return tm_scores

# Main script execution
if __name__ == "__main__":

    for SUBDIR in SUBDIRS:

        print(SUBDIR)
        PDBS_SUBDIR = PDBS_DIR + sep + SUBDIR
        if SUBDIR == "monomer":
            WT_PDB = PDBS_SUBDIR + sep + "wt.pdb"
        else:
            WT_PDB = PDBS_SUBDIR + sep + "wt-{}.pdb".format(SUBDIR)

        # Compute TM-scores for all mutants in the pdbs directory
        tm_scores = compute_tm_scores(WT_PDB, PDBS_SUBDIR)

        # Print the results
        if tm_scores:
            print("\nTM-scores for mutant PDBs:")
            for mutant1, mutant2_s in tm_scores.items():
                for mutant2, score in mutant2_s.items():
                    print(f"{mutant1} vs {mutant2}: TM-score = {score:.4f}")
        else:
            print("No TM-scores were computed.")
        
        # Save the results to a CSV file
        output_file = RESULTS_DIR + sep + "tm_scores_{}_all.csv".format(SUBDIR)
        df = pd.DataFrame.from_dict(tm_scores, orient="index2")
        df.to_csv(output_file, index=False)

    # Quit PyMOL when done
    pymol.cmd.quit()