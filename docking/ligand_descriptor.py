#!/usr/bin/env python3.9

"""
Molecular descriptor calculation library
Converted from Python 2 to Python 3.9
"""

import sys
import time
from typing import Tuple, List, Union, Optional
import os 
import pandas as pd
from rdkit.Chem import SDWriter

cwd = os.getcwd()
sep = os.sep
ligand_types = ["optimized", "existing", "generated"]
ligand_path = cwd + sep + "sdf-{}"

try:
    import openbabel
    from openbabel import pybel
except ImportError:
    print("OpenBabel and pybel are required. Install with: pip install openbabel-wheel")
    sys.exit(1)


def pains_finder(smile: str, smart_def_file: str) -> Tuple[str, List[str], List[List[int]]]:
    """
    Find PAINS (Pan Assay Interference Compounds) in a molecule.
    
    Args:
        smile: SMILES string of the molecule
        smart_def_file: Path to SMARTS definitions file
        
    Returns:
        Tuple of (pains_match, substr_name, highlight)
    """
    pymol = pybel.readstring("smi", smile)
    pymol.removeh()
    
    pains_match = 'No'
    substr_name = []
    highlight = []
    
    with open(smart_def_file, 'r') as input_file:
        for line in input_file.readlines():
            data = line.split()
            smarts = pybel.Smarts(data[0])
            matches = smarts.findall(pymol)
            
            if matches:
                pains_match = 'Yes'
                atoms = []
                
                for i in range(len(matches)):
                    for j in range(len(matches[i])):
                        atoms.append(matches[i][j])
                
                atoms.sort()
                highlight.append(atoms)
                
                info = ''
                for i in range(1, len(data)):
                    if i == 1:
                        info = data[i]
                    else:
                        info = info + ' ' + data[i]
                
                substr_name.append(info)
    
    return pains_match, substr_name, highlight


def get_inchi(smile: str) -> str:
    """Convert SMILES to InChI format."""
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("smi", "inchi")
    mol = openbabel.OBMol()
    conv.ReadString(mol, smile)
    inchi = conv.WriteString(mol)
    return inchi.rstrip('\n')


def get_inchi_key(smile: str) -> str:
    """Convert SMILES to InChI Key format."""
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("smi", "inchikey")
    mol = openbabel.OBMol()
    conv.ReadString(mol, smile)
    inchikey = conv.WriteString(mol)
    return inchikey.rstrip('\n')

def rings_count(obmol: openbabel.OBMol) -> Tuple[int, int, int]:
    """
    Count rings in molecule.
    
    Args:
        obmol: OpenBabel molecule object
        
    Returns:
        Tuple of (total_rings, aromatic_rings, aliphatic_rings)
    """
    n_rings = 0
    n_aromatic_rings = 0
    n_aliphatic_rings = 0
    
    for ring in openbabel.OBMolRingIter(obmol):
        n_rings += 1
        if ring.IsAromatic():
            n_aromatic_rings += 1
        else:
            n_aliphatic_rings += 1
    
    return n_rings, n_aromatic_rings, n_aliphatic_rings


def log_bb(log_p: float, tpsa: float) -> float:
    """
    Calculate blood-brain barrier permeability.
    
    Based on Vilar, S.; Chakrabarti, M.; Costanzi, S.
    Prediction of passive blood-brain partitioning.
    J Mol Graph Model, 2010, 8, 899-903.
    
    Args:
        log_p: Partition coefficient
        tpsa: Topological polar surface area
        
    Returns:
        Blood-brain barrier permeability score
    """
    model_1 = 0.5159 * log_p - 0.0277 * tpsa - 0.3462
    return model_1


def lipinski(hba: int, hbd: int, log_p: float, mw: float) -> Tuple[bool, bool]:
    """
    Calculate Lipinski's Rule of 3 and Rule of 5.
    
    Based on Lipinski, C.A.; Lombardo, F.; Dominy, B.W.; Feeney, P.J.
    Experimental and computational approaches to estimate solubility and
    permeability in drug discovery and development settings.
    Adv Drug Del Rev, 1997, 46, 3-26
    
    Args:
        hba: Hydrogen bond acceptors
        hbd: Hydrogen bond donors
        log_p: Partition coefficient
        mw: Molecular weight
        
    Returns:
        Tuple of (rule_of_3, rule_of_5)
    """
    n = 0
    
    if hba <= 5:
        n += 1
    if hbd <= 10:
        n += 1
    if log_p <= 5.0:
        n += 1
    if mw <= 500.0:
        n += 1
    
    rule_of_3 = n >= 3
    rule_of_5 = n == 4
    
    return rule_of_3, rule_of_5

def chiral_atoms(obmol: openbabel.OBMol) -> bool:
    """
    Check if molecule is chiral.
    
    Args:
        obmol: OpenBabel molecule object
        
    Returns:
        True if molecule is chiral
    """
    # Perceive stereochemistry from 3D coordinates
    openbabel.StereoFrom3D(obmol)

    # Now check for chiral atoms
    chiral_atoms = []
    for atom in openbabel.OBMolAtomIter(obmol):
        if atom.IsChiral():
            chiral_atoms.append(atom.GetId())
    return len(chiral_atoms) > 0, len(chiral_atoms)

def compute_hydrogen_bond_acceptors(smile: str) -> int:
    """
    Compute number of hydrogen bond acceptors in a SMILES string.
    
    Args:
        smile: SMILES string    
    
    Returns:
        Number of hydrogen bond acceptors
    """
    
    aa_obmol      = obmol
    #aa_obmol.AddHydrogens()
    n_HBA = 0

    for obatom in openbabel.OBMolAtomIter( aa_obmol ):
        if obatom.IsHbondAcceptor():
            n_HBA = n_HBA + 1
    
    return n_HBA

def is_chirality_specified(smile: str, n_chiral_atms: int) -> bool:
    """
    Check if chirality is specified in SMILES string.
    
    Args:
        smile: SMILES string
        n_chiral_atms: Number of chiral atoms
        
    Returns:
        True if chirality is specified
    """
    if n_chiral_atms < 1:
        return False
    
    defined_atms = 0
    def_pos = -2
    
    for n in range(len(smile)):
        cha = smile[n]
        if cha == '@':
            if def_pos != (n - 1):
                defined_atms += 1
            def_pos = n
    
    return defined_atms == n_chiral_atms


def double_bonds_ana(obmol: openbabel.OBMol, d_bonds: int, smile: str) -> Tuple[int, int, bool]:
    """
    Analyze double bonds in molecule.
    
    Args:
        obmol: OpenBabel molecule object
        d_bonds: Number of double bonds
        smile: SMILES string
        
    Returns:
        Tuple of (d_bonds_in_ring, d_bonds_out_ring, defined_EZ)
    """
    d_bonds_in_ring = 0
    d_bonds_out_ring = 0
    defined_ez = False
    n_def_signs = 0
    
    if d_bonds > 0:
        for bond in openbabel.OBMolBondIter(obmol):
            if bond.IsDouble():
                if bond.IsInRing():
                    d_bonds_in_ring += 1
                else:
                    sa = bond.GetBeginAtom().IsCarbon()
                    fa = bond.GetEndAtom().IsCarbon()
                    if sa and fa:
                        d_bonds_out_ring += 1
        
        if d_bonds_out_ring > 0:
            for n in range(len(smile)):
                cha = smile[n]
                if cha == '/' or cha == '\\':
                    n_def_signs += 1
            
            if n_def_signs == (d_bonds_out_ring * 2) or n_def_signs - 1 == (d_bonds_out_ring * 2):
                defined_ez = True
    
    return d_bonds_in_ring, d_bonds_out_ring, defined_ez


def is_ketone_group(bond: openbabel.OBBond) -> bool:
    """Check if bond is part of a ketone group."""
    atom_a = bond.GetBeginAtom()
    atom_type_a = atom_a.GetType()
    atom_b = bond.GetEndAtom()
    atom_type_b = atom_b.GetType()
    
    if atom_type_a == 'C2':
        carbon = atom_a
    else:
        carbon = atom_b
    
    if carbon.GetHvyValence() < 3:
        return False
    
    if carbon.GetHvyValence() == 3:
        for neighbour_atom in openbabel.OBAtomAtomIter(carbon):
            if neighbour_atom.GetType() == 'O3':
                return False
    
    return True


def is_aldehyde_group(bond: openbabel.OBBond) -> bool:
    """Check if bond is part of an aldehyde group."""
    atom_a = bond.GetBeginAtom()
    atom_type_a = atom_a.GetType()
    atom_b = bond.GetEndAtom()
    atom_type_b = atom_b.GetType()
    
    if atom_type_a == 'C2':
        carbon = atom_a
    else:
        carbon = atom_b
    
    return carbon.GetHvyValence() <= 2


def is_not_carboxyl_group(bond: openbabel.OBBond) -> bool:
    """Check if bond is not part of a carboxyl group."""
    atom_a = bond.GetBeginAtom()
    atom_type_a = atom_a.GetType()
    atom_b = bond.GetEndAtom()
    atom_type_b = atom_b.GetType()
    
    if atom_type_a == 'O2':
        oxygen = atom_a
    else:
        oxygen = atom_b
    
    return not oxygen.IsCarboxylOxygen()


def get_fingerprint(obmol: openbabel.OBMol) -> str:
    """Generate molecular fingerprint."""
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetOutFormat("fpt")
    ob_conversion.SetOptions("xs", ob_conversion.OUTOPTIONS)
    fpt = ob_conversion.WriteString(obmol)
    return fpt[1:]


def format_a_to_format_b(format_in: str, format_out: str, mol_in: str) -> str:
    """Convert molecule from one format to another."""
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats(format_in, format_out)
    obmol_out = openbabel.OBMol()
    ob_conversion.ReadString(obmol_out, mol_in)
    mol_out = ob_conversion.WriteString(obmol_out)
    return mol_out


def build_3d(format_type: str, smile: str, filename: str) -> str:
    """Build 3D structure from SMILES."""
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", format_type)
    obmol = openbabel.OBMol()
    ob_conversion.ReadString(obmol, smile)
    
    bd = openbabel.OBBuilder()
    bd.Build(obmol)
    obmol.SetDimension(3)
    obmol.AddHydrogens(False, True, 7)
    
    ff = openbabel.OBForceField.FindForceField("MMFF94")
    ff.Setup(obmol)
    ff.SteepestDescent(150, 1.0e-4)
    ff.WeightedRotorSearch(5, 25)
    ff.ConjugateGradients(250, 1.0e-6)
    ff.UpdateCoordinates(obmol)
    
    # Corrected output handling
    ob_conversion.WriteFile(obmol, filename)  # Directly write to file
    mol3d = ob_conversion.WriteString(obmol)  # Return string representation
    return mol3d


def build_2d(format_type: str, smile: str) -> str:
    """Build 2D structure from SMILES."""
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", format_type)
    obmol = openbabel.OBMol()
    ob_conversion.ReadString(obmol, smile)
    
    gen2d = openbabel.OBOp.FindType("Gen2D")
    gen2d.Do(obmol)
    obmol.DeleteHydrogens()
    obmol.SetTitle("CMLDID")
    
    mol2d = ob_conversion.WriteString(obmol)
    return mol2d


def is_complex(smile: str) -> bool:
    """
    Check if molecule is complex based on size and flexibility.
    
    Args:
        smile: SMILES string
        
    Returns:
        True if molecule is considered complex
    """
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInFormat("smi")
    obmol = openbabel.OBMol()
    ob_conversion.ReadString(obmol, smile)
    
    n_ha_atms = obmol.NumHvyAtoms()
    n_ha_rb = obmol.NumRotors()
    
    return n_ha_atms > 15 or n_ha_rb > 5


def check_web_data(web_user_rightconfig_spec: str, web_user_delivery: str, 
                  chiral: bool, specified_chirality: bool, 
                  d_bonds_out_ring: int, defined_ez: bool) -> bool:
    """
    Validate molecular data for web interface.
    
    Args:
        web_user_rightconfig_spec: User specification for configuration
        web_user_delivery: Delivery type
        chiral: Whether molecule is chiral
        specified_chirality: Whether chirality is specified
        d_bonds_out_ring: Number of double bonds outside rings
        defined_ez: Whether E/Z geometry is defined
        
    Returns:
        True if validation passes
        
    Raises:
        SystemExit: If validation fails
    """
    if web_user_rightconfig_spec == "NotChiral" and d_bonds_out_ring > 0.0:
        print("ERROR| Double bond(s) detected but \"No isomers\" has been declared")
        sys.exit()
    
    if web_user_rightconfig_spec == "NotChiral" and d_bonds_out_ring > 0.0 and not defined_ez:
        print("ERROR| Double bond(s) detected but geometry has been not specified in molecule structure")
        sys.exit()
    
    if web_user_rightconfig_spec == "NotChiral" and chiral:
        print("ERROR| Chiral atom(s) detected but \"No isomers\" has been declared")
        sys.exit()
    
    if web_user_delivery == "Theoretical":
        if chiral and not specified_chirality:
            print("ERROR| Chirality must be specified for theoretical compounds")
            sys.exit()
        
        if d_bonds_out_ring > 0.0 and not defined_ez:
            print("ERROR| Double bonds geometry must be specified for theoretical compounds")
            sys.exit()
        
        if web_user_rightconfig_spec == "No":
            print("ERROR| Theoretical entry must be declared (and designed) as right chirality and/or double bond geometry or \"No isomers\" ")
            sys.exit()
    
    if web_user_rightconfig_spec == "No" and web_user_delivery == "Pure":
        if chiral or d_bonds_out_ring > 0.0:
            print("ERROR| Your compound has unspecified chiral atom(s) and/or double bond(s), it can delivered as Mixture")
            sys.exit()
    
    if web_user_rightconfig_spec == "Yes":
        if not chiral and d_bonds_out_ring == 0.0:
            print("ERROR| No chiral atoms or double bonds detected, please specify \"No isomers\" to the question on the right design or modify your compound")
            sys.exit()
        
        if chiral and not specified_chirality:
            print("ERROR| Chiral atom(s) detected but the configuration has been not graphically specified")
            sys.exit()
        
        if d_bonds_out_ring > 0.0 and not defined_ez:
            print("ERROR| Double bond(s) detected but the geometry has been not graphically specified")
            sys.exit()
    
    return True


# Add this function before your main block
def calculate_druglikeness(smile):
    """
    Calculate drug-likeness score using RDKit's QED implementation
    
    Args:
        smile (str): SMILES string of the molecule
    
    Returns:
        dict: Dictionary containing QED score and molecular descriptors
    """
    from rdkit import Chem
    from rdkit.Chem.Descriptors import MolWt, MolLogP, NumHDonors, NumHAcceptors, TPSA, NumRotatableBonds
    from rdkit.Chem import QED
    
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        return None
    
    qed_score = QED.qed(mol)
    mw = MolWt(mol)
    logp = MolLogP(mol)
    hbd = NumHDonors(mol)
    hba = NumHAcceptors(mol)
    tpsa = TPSA(mol)
    rotb = NumRotatableBonds(mol)
    
    return {
        'qed_score': qed_score,
        'mw': mw,
        'logp': logp,
        'hbd': hbd,
        'hba': hba,
        'tpsa': tpsa,
        'rotb': rotb
    }

if __name__ == "__main__":

    dict = {}

    for type in ligand_types:
        path = ligand_path.format(type)
        #print(f"Processing {type} ligands in {path}")
        ligands = [file for file in os.listdir(path) if file.endswith('.sdf')]
        if type in ["optimized", "generated"]:
            ligands = [ligand for ligand in ligands if ligand.startswith("best_")]
            ligands = [f for f in ligands if not f.endswith("_preprocessed.sdf")]
        #print(f"Found {len(ligands)} ligands in {type} directory.")
        for ligand in ligands:

            if "_preprocessed" in ligand:
                continue
            ligand_path_full = os.path.join(path, ligand)
            try:
                mol = pybel.readfile("sdf", ligand_path_full).__next__()
                smile = mol.write("smi").strip().split("\t")[0]
                print(f"Processing {ligand} with SMILES: {smile}")
                
                obmol = pybel.readstring("smi", smile).OBMol
                chiral, n_chiral_atoms = chiral_atoms(obmol)
                is_specified = is_chirality_specified(smile, n_chiral_atoms)
                #print(f"Chirality specified: {is_specified}, Chiral atoms: {n_chiral_atoms}, Chiral: {chiral}")

                n_rings, n_aromatic_rings, n_aliphatic_rings = rings_count(obmol)
                #print(f"Total rings: {n_rings}, Aromatic rings: {n_aromatic_rings}, Aliphatic rings: {n_aliphatic_rings}")

                log_p = mol.calcdesc()['logP']
                tpsa = mol.calcdesc()['TPSA']
                bb_permeability = log_bb(log_p, tpsa)
                #print(f"LogP: {log_p}, TPSA: {tpsa}, Blood-brain barrier permeability: {bb_permeability}")

                hba = compute_hydrogen_bond_acceptors(smile)
                hbd = mol.calcdesc()['HBD']
                mw = mol.calcdesc()['MW']
                #print(f"HBA: {hba}, HBD: {hbd}, MW: {mw}")
                
                lipinski_rule_of_3, lipinski_rule_of_5 = lipinski(hba, hbd, log_p, mw)
                #print(f"Lipinski's Rule of 3: {lipinski_rule_of_3}, Rule of 5: {lipinski_rule_of_5}")

                ligand_preprocessed_path_full = os.path.join(path, f"{ligand.split('.')[0]}_preprocessed.sdf")
                mol_3d = build_3d("sdf", smile, ligand_preprocessed_path_full)

                # Calculate drug-likeness score
                druglikeness = calculate_druglikeness(smile)
                
                dict[ligand] = {
                    "smile": smile,
                    "chiral": chiral,
                    "n_chiral_atoms": n_chiral_atoms,
                    "is_specified": is_specified,
                    "n_rings": n_rings,
                    "n_aromatic_rings": n_aromatic_rings,
                    "n_aliphatic_rings": n_aliphatic_rings,
                    "log_p": log_p,
                    "tpsa": tpsa,
                    "bb_permeability": bb_permeability,
                    "hba": hba,
                    "hbd": hbd,
                    "mw": mw,
                    "lipinski_rule_of_3": lipinski_rule_of_3,
                    "lipinski_rule_of_5": lipinski_rule_of_5,

                    "qed_score": druglikeness['qed_score'] if druglikeness else None,
                    #"druglikeness_mw": druglikeness['mw'] if druglikeness else None,
                    #"druglikeness_logp": druglikeness['logp'] if druglikeness else None,
                    #"druglikeness_hbd": druglikeness['hbd'] if druglikeness else None,
                    #"druglikeness_hba": druglikeness['hba'] if druglikeness else None,
                    #"druglikeness_tpsa": druglikeness['tpsa'] if druglikeness else None,
                    #"druglikeness_rotb": druglikeness['rotb'] if druglikeness else None,
                }

            except Exception as e:
                print(f"Error processing {ligand}: {e}")

    df = pd.DataFrame.from_dict(dict, orient='index')
    output_file = cwd + sep + "ligand_descriptors.csv"
    df.to_csv(output_file, index_label='ligand')
    print(f"Descriptors saved to {output_file}")
