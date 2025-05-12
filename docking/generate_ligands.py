import subprocess
import os 

cwd = os.getcwd()
sep = os.sep 


parent_dir = os.path.join(cwd, os.pardir)
parent_dir = os.path.abspath(parent_dir)
pdb_dir = parent_dir + sep + "pdbs-alphafold" + sep + "tetramer"
mutation = "y98f"  # Specify the mutation here
mutant_path = pdb_dir + sep + f"{mutation}-tetramer.pdb"


resis = [15, 17, 54, 106, 108, 110, 117, 119] 
chains =  ["C", "D"]

residues = ""
for chain in chains:
    for resi in resis:
        residues += f"{chain}:{resi} " 
residues = residues.rstrip()

if os.path.basename(os.path.normpath(os.getcwd())) != "docking":
    os.chdir("docking")

os.makedirs(cwd + sep + "sdf-generated", exist_ok=True)

model = "crossdocked_fullatom_cond.ckpt"

# Dictionary of reference residues
residues_reference = {
    f"{mutation}_T4": residues
}


# Base command
base_cmd = "python DiffSBDD/generate_ligands.py DiffSBDD/checkpoints/{model} --pdbfile {mutant_path} --outfile sdf-generated/{name}.sdf --resi_list {residues} --n_samples 100 --timesteps 100 --sanitize"

# Iterate through ligands and execute the command
for name, residues in residues_reference.items():
    cmd = base_cmd.format(model=model, mutant_path=mutant_path, name=name, residues=residues)
    print(f"Running: {cmd}")
    subprocess.run(cmd, shell=True)#, check=True)
