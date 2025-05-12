import argparse
from pathlib import Path
import os 
import pandas as pd 
import numpy as np 
import subprocess
import time 
from rdkit import Chem
from rdkit.Chem import AllChem
import ast 
import re
from Bio import PDB
import pymol
from pymol import cmd

complex_types = ["tetramer"]#["monomer", "dimer", "tetramer"]



type_sdf = "optimized" #existing, generated or optimized

cwd = os.getcwd()
sep = os.sep 

sdf_dir = cwd + sep + f"sdf-{type_sdf}" #sdf-existing, sdf-optimized or sdf-generated

#if "_" in type_sdf:
#    type_sdf = type_sdf.split("_")[0]

docker = "autodock"

parent_dir = os.path.join(cwd, os.pardir)
parent_dir = os.path.abspath(parent_dir)
print(parent_dir)

pdb_dir = parent_dir + sep + "pdbs-alphafold" 
results_path = parent_dir + sep + "results" 
docked_structures = cwd + sep + "docked-outputs" + sep + type_sdf

diffdock_path = cwd + sep + "DiffDock" 
diffdock_inference = diffdock_path + sep +  "inference.py"
diffdock_conf = diffdock_path + sep + "default_inference_args.yaml"
diffdock_inference_cmd = "python {} --config {} --protein_ligand_csv {} --out_dir {}"

autodock_path = cwd + sep + "AutoDock" 
autodock_prepare_path = autodock_path + sep + "prepare_outputs" + sep + type_sdf
autodock_docked_outputs = autodock_path + sep + "docked-outputs" + sep + type_sdf
autodock_inference_cmd = "vina --receptor {} --ligand {} --out {} --exhaustiveness 32 --scoring {} --config {} --verbosity 2 > {}"
autodock_log_path = autodock_path + sep + "log" + sep + type_sdf
default_docking = "autodock"

os.makedirs(docked_structures, exist_ok=True)
os.makedirs(autodock_docked_outputs, exist_ok=True)
os.makedirs(autodock_path, exist_ok=True)
os.makedirs(autodock_log_path, exist_ok=True)
os.makedirs(autodock_prepare_path, exist_ok=True)
os.makedirs(results_path, exist_ok=True)

score_from_diffdock = False #if True, apply autodock vina on diffdock pose

def extract_vina_scores(log_file_path):
    scores = []
    pattern = re.compile(r"\s*\d+\s+([-\d\.]+)\s+")
    
    with open(log_file_path, 'r') as file:
        for line in file:
            match = pattern.match(line)
            if match:
                scores.append(float(match.group(1)))
    
    return scores

def calculate_center(pdb_file):
    # Parse the PDB file using Biopython's PDB parser
    parser = PDB.PPBuilder()
    structure = PDB.PDBParser(QUIET=True).get_structure('mol', pdb_file)
    
    # Initialize variables to accumulate coordinates
    total_x, total_y, total_z = 0, 0, 0
    atom_count = 0
    chains = ["A", "D"]
    # Iterate over all atoms in the structure
    for model in structure:
        for chain in model:
            resi_count = 0
            if chain.id in chains:
                for residue in chain:
                    resi_count += 1
                    if resi_count+1 == 109:
                        #print(residue)
                        for atom in residue:
                            # Accumulate the x, y, z coordinates of each atom
                            x, y, z = atom.coord
                            total_x += x
                            total_y += y
                            total_z += z
                            atom_count += 1
                    
    # Calculate the center by averaging the coordinates
    if atom_count > 0:
        center_x = total_x / atom_count
        center_y = total_y / atom_count
        center_z = total_z / atom_count
        return center_x, center_y, center_z
    else:
        raise ValueError("No atoms found in the PDB file")

def generate_vina_box(center, box_size):
    # Unpack center coordinates
    center_x, center_y, center_z = center
    
    # Unpack box size
    box_x, box_y, box_z = box_size
    
    # Format the Vina box-center and box-size parameters
    vina_center = f"{center_x:.2f} {center_y:.2f} {center_z:.2f}"
    vina_size = f"{box_x} {box_y} {box_z}"
    
    return vina_center, vina_size

def extract_ligand(json_file):
    with open(json_file, "r") as file:
        content = file.read()
    dictionary = ast.literal_eval(content)
    ligand_positions = dictionary.get("ligand_positions", [])[0]
    return ligand_positions

def prepare_sdf_ligand(ligand_position, ligandpath):
    """
    Prepares an SDF ligand file by ensuring 3D coordinates, adding hydrogens,
    and writing the processed molecules to a new file.

    Parameters:
        ligand_position (str): Path to the input SDF file or MOL block string.
        ligandpath (str): Path to save the processed SDF file.

    Returns:
        str: Path to the generated PDB file if successful, else None.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolops

    if isinstance(ligand_position, str):
        # Read molecules from the input SDF file
        supplier = Chem.SDMolSupplier(ligand_position, sanitize=False, strictParsing=True)
        valid_mols = []

        for idx, mol in enumerate(supplier):
            if mol is None:
                print(f"⚠️ Skipping invalid molecule at index {idx}")
                continue

            try:
                # Ensure 3D coordinates exist
                if mol.GetNumConformers() == 0:
                    print(f"Generating 3D coordinates for molecule at index {idx}")
                    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                else:
                    # Check if Z-coordinates are non-zero and mark as 3D
                    conf = mol.GetConformer()
                    if any(conf.GetAtomPosition(i).z != 0 for i in range(mol.GetNumAtoms())):
                        conf.Set3D(True)

                # Add hydrogens with proper handling of coordinates
                mol = rdmolops.AddHs(mol, addCoords=True)

                # Validate molecule state before saving
                if mol.GetNumAtoms() == 0:
                    raise ValueError("Molecule lost all atoms during processing")

                valid_mols.append(mol)

            except Exception as e:
                print(f"❌ Molecule {idx} failed processing: {str(e)}")
                continue

        # Write processed molecules to the output SDF file
        if valid_mols:
            with Chem.SDWriter(ligandpath) as writer:
                writer.SetKekulize(True)  # Preserve explicit double bonds
                for mol in valid_mols:
                    writer.write(mol)
            print(f"✅ Successfully processed {len(valid_mols)} molecules")
        else:
            print("⚠️ No valid molecules processed")
        return None

    else:
        # Handle MOL block input directly
        mol = Chem.MolFromMolBlock(ligand_position, sanitize=False)
        if mol is None:
            raise ValueError("Failed to parse molecule from input block")

        # Ensure 3D coordinates exist
        if mol.GetNumConformers() == 0:
            print("Generating 3D coordinates for input MOL block")
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        # Add hydrogens with proper handling of coordinates
        mol = rdmolops.AddHs(mol, addCoords=True)

        # Write processed molecule to the output SDF file
        with Chem.SDWriter(ligandpath) as writer:
            writer.SetKekulize(True)
            writer.write(mol)

        # Generate a PDB file for docking tools like Meeko or AutoDock-Vina
        filename_pdb = ligandpath.rsplit(".", 1)[0] + ".pdb"
        Chem.MolToPDBFile(mol, filename_pdb)
        return filename_pdb

def Receptor3DView(receptorPDB=None, boxPDB=None, ligPDB=None):
    pymol.finish_launching()  # Launch PyMOL if it's not already running

    if boxPDB:
        cmd.load(boxPDB, "box")
        cmd.show("sticks", "box")

    if receptorPDB:
        cmd.load(receptorPDB, "receptor")
        cmd.show("cartoon", "receptor")
        cmd.spectrum("b", "rainbow", "receptor")
        cmd.set("cartoon_transparency", 0.5, "receptor")

    if ligPDB:
        cmd.load(ligPDB, "ligand")
        cmd.show("sticks", "ligand")
    
    cmd.orient()  # Adjust the view for best fit
    cmd.zoom()    # Zoom to fit the loaded molecules

if __name__ == "__main__":

    if docker not in ["diffdock", "autodock"]:
        docker = "autodock"

    ligands = [file for file in os.listdir(sdf_dir) if file.split(".")[-1] == "sdf"]
    #ligands = ["tafamidis.sdf"]
    print(ligands)

    for complex_type in complex_types:
        print(complex_type)
        protein_ligand_csv = cwd + sep + "protein_ligand_{}.csv".format(complex_type)
        pdb_subdir = pdb_dir + sep + complex_type
        mutants = [file for file in os.listdir(pdb_subdir) if file.split(".")[-1] == "pdb"]
        #mutants = ["v142i-tetramer.pdb"]
        mutations = [mutant.split(".")[0] for mutant in mutants]

        for ligand in ligands:
            print("Ligand: ", ligand)
            df_rows = []
            columns = ["Mutation"]
            if docker == "autodock":
                n = 9
            else:
                n = 20
            [columns.append(f"Score_{i+1}") for i in range(n)]
            for mutant in mutants:  
                complex_name = mutant.split(".")[0]+"_"+ligand.split(".")[0]
                t_start = time.time()
                if docker == "diffdock":
                    continue# DO IT WITH NVIDIA APIs: use docking_diffdock.py script
                    output_name = "DiffDock"
                    docking_cmd = diffdock_inference_cmd.format(diffdock_inference, diffdock_conf, protein_ligand_csv, docked_structures + sep + output_name)
                    subprocess.run(docking_cmd, shell=True)
                    docking_results = os.listdir(docked_structures + sep + output_name + sep + complex_name)
                    diffdock_scores = [filename.replace(".sdf", "").split("-")[-1] for filename in docking_results]
                    print("Mutant: {}, Confidence: {}".format(mutant, diffdock_scores))
                else:
                    pdb_filepath = pdb_dir + sep + complex_type + sep + mutant
                    sdf_filepath = sdf_dir + sep + ligand
                    pdbqt_rec_filepath = autodock_prepare_path + sep + "{}".format(mutant.split(".")[0])
                    pdbqt_lig_filepath = autodock_prepare_path + sep + "{}.pdbqt".format(ligand.split(".")[0])

                    autodock_out = mutant.split(".")[0]+"_"+ligand.split(".")[0] 
                    logfile = autodock_log_path + sep + f"{autodock_out}.txt"
                    if f"{autodock_out}.txt" in os.listdir(autodock_log_path):
                        print("Log file already exists for {}. Skip docking, extracting vina scores from logs.".format(autodock_out)) 
                    else:
                        
                        if score_from_diffdock:
                            ligand_position = extract_ligand(f"{docked_structures}{sep}{autodock_out}.json")      
                            ligandsdf =  autodock_prepare_path + sep + "{}_{}.sdf".format(ligand.split(".")[0], mutant.split(".")[0])
                            ligandpdb = prepare_sdf_ligand(ligand_position, ligandsdf)
                        else:
                            ligandsdf = sdf_filepath + sep + ligand
                            print(ligandsdf, sdf_filepath)
                            ligandpdb = prepare_sdf_ligand(sdf_filepath, sdf_filepath)

                        ligandfile =  autodock_docked_outputs + sep + "{}_{}.pdbqt".format(mutant.split(".")[0], ligand.split(".")[0])
                        
                        if "tafamidis" in ligand:
                            box_size = 15, 20, 20
                            scoring = "vinardo"
                        else:
                            box_size = 15, 20, 20
                            scoring = "vinardo"
                   
                        center = calculate_center(pdb_filepath)
                        vina_center, vina_size = generate_vina_box(center, box_size)
                        #print(vina_size, vina_center)
                        prepare_rec_cmd = "mk_prepare_receptor.py -i {} -o {} -p --allow_bad_res -v \
                                        --box_size {} --box_center {}".format(pdb_filepath, pdbqt_rec_filepath, vina_size, vina_center)
                        prepare_lig_cmd = "mk_prepare_ligand.py -i {} -o {}".format(sdf_filepath, ligandfile)

                        print("Preparing Ligand")
                        subprocess.call(prepare_lig_cmd, shell=True)
                        print("Preparing Receptor")
                        subprocess.call(prepare_rec_cmd, shell=True)
            
                        outfile = autodock_docked_outputs + sep + autodock_out + ".pdbqt"
                        mutantfile = autodock_prepare_path + sep + "{}.pdbqt".format(mutant.split(".")[0])
                        print("Mutant file:", mutantfile)
                        print("Ligand file:", ligandfile)
                        print("Output file:", outfile)
                        print("Log file:", logfile)

                        print(mutant)
                        ligandpdb =  autodock_docked_outputs + sep + f"{mutant.split('.')[0]}_{ligand.split('.')[0]}.pdbqt"

                        #box_path = f"{autodock_prepare_path}{sep}{mutant.split('.')[0]}.box.pdb"
                        #Receptor3DView(receptorPDB=pdb_filepath, boxPDB=box_path, ligPDB=ligandpdb)
                        
                        configfile = f"{autodock_prepare_path}{sep}{mutant.split('.')[0]}.box.txt"
                        vina_cmd = autodock_inference_cmd.format(mutantfile, ligandfile, outfile, scoring, configfile, logfile)
                        print("Executing command:", vina_cmd)
                        subprocess.call(vina_cmd, shell=True)

                    scores = extract_vina_scores(logfile)
                    print("Mutant: {}, Vina score: {}".format(mutant, scores))
                    t_end = time.time() - t_start
                    print("Execution time for {}_{}: {} seconds".format(mutant.split(".")[0], ligand.split(".")[0], t_end))
                row = [mutant]
                [row.append(score) for score in scores]
                df_rows.append(row)
    
            df = pd.DataFrame(df_rows, columns=columns)
            df.set_index('Mutation', inplace=True)
            print(df)
            if docker == "autodock":
                df.to_csv(autodock_docked_outputs + sep + "vina_score_docking_{}_TTR{}.csv".format(ligand.split(".")[0], complex_type))
            else:
                continue
                #df.to_csv(docked_structures + sep + "diffdock_score_docking_{}_TTR{}.csv".format(ligand.split(".")[0], complex_type))              