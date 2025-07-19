import os 
from rdkit import Chem
from rdkit.Chem import QED
import pandas as pd


def hasSulfurAtoms(mol):
    """
    Check if a molecule contains sulfur atoms.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
    Returns:
        bool: True if sulfur atoms are present, False otherwise
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Atomic number for sulfur
            return True
    return False

def hasMonovalentCarbon(mol):
    """
    Check if a molecule contains monovalent carbon atoms.
    
    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
    Returns:
        bool: True if monovalent carbon atoms are present, False otherwise
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalValence() == 1:  # Atomic number for carbon
            return True
    return False

def process_ligands(input_path, output_path, remove_sulfur=False, remove_monovalent_C=True):
    """
    Processes ligands from SDF file, computes QED scores, 
    and saves the best ligand.
    
    Args:
        input_path (str): Path to input SDF file
        output_path (str): Path to save best ligand SDF
    Returns:
        float: QED score of best ligand
    """
    # Read SDF file
    suppl = Chem.SDMolSupplier(input_path)
    best_mol = None
    best_qed = -1.0
    
    # Process molecules
    for i, mol in enumerate(suppl):
        print(f"Processing molecule {i+1}/{len(suppl)}")
        if mol is not None:
            if remove_sulfur:
                # Check for sulfur atoms
                if hasSulfurAtoms(mol):
                    continue  # Skip molecules with sulfur atoms, we cannot use them in molecular dynamic simulations due to CGenff limitations
            if remove_monovalent_C:
                if hasMonovalentCarbon(mol):
                    print(f"Skipping molecule {i+1} due to monovalent carbon atoms.")
                    continue  # Skip molecules with monovalent carbon atoms

            # Compute QED score
            qed_score = QED.default(mol)
            if qed_score > best_qed:
                best_qed = qed_score
                best_mol = mol
    
    # Save best ligand
    if best_mol:
        for atom in best_mol.GetAtoms():
            print(f"Atom: {atom.GetSymbol()}, Valence: {atom.GetTotalValence()}")
        writer = Chem.SDWriter(output_path)
        writer.write(best_mol)
        writer.close()
    
    return best_qed


cwd = os.getcwd()
sep = os.sep 

sdf_types = ["existing", "generated", "optimized"]

if __name__ == "__main__":
    
    remove_sulfur = False  # Set to True to remove ligands with sulfur atoms
    remove_monovalent_C = False # Set to True to remove monovalent carbon atoms

    qed_dict = {}
    for sdf_type in sdf_types:
        ligands_path = cwd + sep + f"sdf-{sdf_type}"
        ligands_files = [f for f in os.listdir(ligands_path) if f.endswith(".sdf")]
        ligands_files = [f for f in ligands_files if not f.startswith("best_") and not f.startswith("ligand_")]
        ligands_files = [f for f in ligands_files if not f.endswith("_preprocessed.sdf")]
        
        for ligand_file in ligands_files:
            ligand_name = ligand_file.split(".")[0]
            ligand_path = ligands_path + sep + ligand_file
            ligand_path_best = ligands_path + sep + "best_{}.sdf".format(ligand_name)    
            # Process the ligand file as needed
            print(f"Processing {ligand_name} from {ligand_path}")
            
            if sdf_type in ["generated", "optimized"]:
                best_score = process_ligands(ligand_path, ligand_path_best, remove_sulfur=remove_sulfur, remove_monovalent_C=remove_monovalent_C)
                print(f"Best QED score for {ligand_name}: {best_score}")
            else:
                # For existing ligands, we assume they are already processed
                best_score = QED.default(Chem.MolFromMolFile(ligand_path))
                print(f"QED score for existing ligand {ligand_name}: {best_score}")

            qed_dict[ligand_name] = best_score
    
    # Save the QED scores to a CSV file
    qed_df = pd.DataFrame(list(qed_dict.items()), columns=["Ligand", "QED"])
    qed_df.to_csv(cwd + sep + "qed_scores.csv", index=False)
    print("QED scores saved to qed_scores.csv")
        
