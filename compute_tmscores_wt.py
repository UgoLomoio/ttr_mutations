import os
import pymol
from tmalign import * 
import pandas as pd 

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
def compute_tm_score(wt_pdb, mutant_pdb):
    # Load the structures into PyMOL
    cmd.load(wt_pdb, "wt")
    mutant = mutant_pdb.split(".")[0].split("/")[-1]
    cmd.load(mutant_pdb, mutant)

    # Align the structures and get the TM-score
    print(exe_path)
    tm_score = tmscore('wt', mutant, quiet = 1, ter = 1, exe=exe_path)

    # Remove the structures after alignment to avoid memory overload
    #cmd.delete("wt")
    cmd.delete(mutant)
   
    return tm_score


# Function to compute TM-scores for all mutants in the pdbs directory
def compute_tm_scores(wt_pdb, pdbs_dir):
    tm_scores = {}
    
    # Loop through all files in the pdbs directory
    for pdb_file in os.listdir(pdbs_dir):
        if pdb_file.endswith('.pdb'):
            mutant_pdb = os.path.join(pdbs_dir, pdb_file)
            print(f"Computing TM-score between {wt_pdb} and {mutant_pdb}...")
            tm_score = compute_tm_score(wt_pdb, mutant_pdb)
            
            # Store the result in the dictionary
            tm_scores[pdb_file.split(".")[0]] = tm_score
    
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
            for mutant, score in tm_scores.items():
                print(f"{mutant}: TM-score = {score:.4f}")
        else:
            print("No TM-scores were computed.")
        
        # Save the results to a CSV file
        output_file = RESULTS_DIR + sep + "tm_scores_{}.csv".format(SUBDIR)
        df = pd.DataFrame(list(tm_scores.items()), columns = ['Mutant', 'TM-score'])
        df.to_csv(output_file, index=False)

    # Quit PyMOL when done
    pymol.cmd.quit()