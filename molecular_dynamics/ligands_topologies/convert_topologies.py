#CODE REQUIRES CHARMM36 FF, NETWORKX 2.3, PYTHON 3.7 AND NUMPY
import subprocess
import os 

cwd = os.getcwd()
sep = os.sep

parent_dir = os.path.dirname(cwd)
charmm36_dir = parent_dir + sep + "charmm36-jul2022.ff" 
cmd = "python cgenff_charmm2gmx.py {} {}.mol2 {}.str {}"
mol_names = {"tafamidis": "3MI", "acoramidis": "16V"}
ligands = ["tafamidis", "acoramidis"]

if __name__ == "__main__":
    
    for ligand in ligands:

        print(f"Converting {ligand} topology files...")
        mol_name = mol_names[ligand]
        print(f"Using {mol_name} as the molecule name")
        
        cmd_formatted = cmd.format(
            mol_name,
            cwd + sep + ligand + sep + ligand,
            cwd + sep + ligand + sep + ligand,
            charmm36_dir
        )
        subprocess.run(cmd_formatted, shell=True, check=True)