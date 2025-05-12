"""
Optimize an existing ligand to better bind mutated variants
"""

import os 
import subprocess 

cwd = os.getcwd()
sep = os.sep 

parent_dir = os.path.dirname(cwd)

dffsbdd_path = cwd + sep + "DiffSBDD"
pdb_path = parent_dir + sep + "pdbs-alphafold" + sep + "tetramer"
sdf_path = cwd + sep + "sdf-existing"
sdf_generated = cwd + sep + "sdf-optimized"


ligands = ["tafamidis", "acoramidis"]
ligand = ligands[0]

mutation_to_optimize = "y98f"
mutated_pdbs = [file for file in os.listdir(pdb_path) if ".pdb" in file]
mutated_pdb = [pdb for pdb in mutated_pdbs if mutation_to_optimize in pdb][0]

outfile = f"{ligand}_optimized_{mutation_to_optimize}"

if __name__ == "__main__":

    cmd = f"python {dffsbdd_path}/optimize.py --checkpoint {dffsbdd_path}/checkpoints/crossdocked_fullatom_cond.ckpt --pdbfile {pdb_path}/{mutated_pdb} --outfile {sdf_generated}{sep}{outfile}.sdf --ref_ligand {sdf_path}/{ligand}.sdf --objective sa --population_size 1 --evolution_steps 1 --top_k 10 --timesteps 100"
    subprocess.call(cmd, shell=True)