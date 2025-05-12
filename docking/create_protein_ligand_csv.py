import os 
import pandas as pd 
import numpy as np 

cwd = os.getcwd()
sep = os.sep 


complex_types = ["monomer", "dimer", "tetramer"]
docking_path = cwd + sep + "docking" 
sdf_existing_path = docking_path + sep + "sdf-existing"
sdf_ligands_path = docking_path + sep + "sdf-ligands"
pdb_path = cwd + sep + "pdbs-alphafold"

pdb_files_dict = {complex_type: [file for file in os.listdir(pdb_path + sep + complex_type) if file.split(".")[-1] == "pdb"] for complex_type in complex_types}
sdf_files = [sdf_existing_path + sep + file for file in os.listdir(sdf_existing_path) if file.split(".")[-1] == "sdf"]
[sdf_files.append(sdf_ligands_path + sep + file) for file in os.listdir(sdf_ligands_path) if file.split(".")[-1] == "sdf"]

if __name__ == "__main__":
    
    for complex_type in complex_types:
        rows = []
        columns = ["complex_name", "protein_path", "ligand_description", "protein_sequence"]
        csv_path = docking_path + sep + "protein_ligand_{}.csv".format(complex_type)
        pdb_files = pdb_files_dict[complex_type]
        for pdb_file in pdb_files:
            pdb_filepath = pdb_path + sep + complex_type + sep + pdb_file
            for sdf_filepath in sdf_files:
                name = pdb_file.split(".")[0] + "_" + sdf_filepath.split(".")[0].split(sep)[-1]
                row = [name, pdb_filepath, sdf_filepath, None]
                rows.append(row)
        df = pd.DataFrame(rows, columns = columns)
        df.to_csv(csv_path)
    