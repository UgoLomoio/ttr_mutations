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

real_pdb = cwd + sep + "3a4d.pdb"
alpha_pdb = PDBS_DIR + sep + "dimer" + sep + "wt-dimer.pdb"

#real_pdb = cwd + sep + "wt_real.pdb"
#alpha_pdb = PDBS_DIR + sep + "tetramer" + sep + "wt-tetramer.pdb"

exe_path = cwd + sep + "tmalign_exe" + sep + "TMalign_cpp"


# Function to compute TM-score for two structures
def compute_tm_score(real_pdb, alpha_pdb):
    # Load the structures into PyMOL
    cmd.load(real_pdb, "realwt")
    cmd.load(alpha_pdb, "alphawt")

    # Align the structures and get the TM-score
    print(exe_path)
    tm_score = tmscore('alphawt', "realwt", quiet = 1, ter = 1, exe=exe_path)

    # Remove the structures after alignment to avoid memory overload
    cmd.delete("alphawt")
    cmd.delete("realwt")
   
    return tm_score


# Main script execution
if __name__ == "__main__":

    # Compute TM-scores
    tm_score = compute_tm_score(real_pdb, alpha_pdb)
    print(f"TM-score = {tm_score:.4f}")


    # Quit PyMOL when done
    pymol.cmd.quit()