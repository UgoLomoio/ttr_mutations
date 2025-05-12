import pymol 
from pymol import cmd 
import os 

def Receptor3DView(receptorPDB=None, boxPDB=None, ligPDB=None):
    pymol.finish_launching()  # Launch PyMOL if it's not already running

    if boxPDB:
        for name, box in boxPDB.items():
            cmd.load(box, name)
            cmd.show("sticks", name)

    if receptorPDB:
        cmd.load(receptorPDB, "receptor")
        cmd.show("cartoon", "receptor")
        cmd.spectrum("chain", "rainbow", "receptor")
        cmd.set("cartoon_transparency", 0.5, "receptor")

    if ligPDB:
        for name, lig in ligPDB.items():
            cmd.load(lig, name)
            cmd.show("sticks", name)
    
    cmd.orient()  # Adjust the view for best fit
    cmd.zoom()    # Zoom to fit the loaded molecules

mutants = ["e74k"]
ligands = ["tafamidis", "acoramidis"]

cwd = os.getcwd()
sep = os.sep
complex_type = "tetramer"
parent_dir = os.path.join(cwd, os.pardir)
parent_dir = os.path.abspath(parent_dir)
pdb_dir = parent_dir + sep + "pdbs-alphafold" 
autodock_path = cwd + sep + "AutoDock" 
type_sdf = "existing"
autodock_prepare_path = autodock_path + sep + "prepare_outputs" + sep + type_sdf
autodock_ligand_path = autodock_path + sep + "docked-outputs" + sep + type_sdf

files = os.listdir()

if __name__ == "__main__":

    for mutant in mutants:  
        pdb_filepath = pdb_dir + sep + complex_type + sep + f"{mutant.split('.')[0]}-{complex_type}.pdb"
        print(pdb_filepath)
        ligandspdb = {}
        boxes = {}
        for ligand in ligands:
            ligandpdb =  autodock_ligand_path + sep + f"{mutant.split('.')[0]}-{complex_type}_{ligand.split('.')[0]}.pdbqt"
            print(ligandpdb)
            ligandspdb[ligand.split('.')[0]] = ligandpdb
            box_path = f"{autodock_prepare_path}{sep}{mutant.split('.')[0]}-{complex_type}.box.pdb"
            print(box_path)
            boxes[f"{ligand.split('.')[0]}-box"] = box_path

        Receptor3DView(receptorPDB=pdb_filepath, boxPDB=boxes, ligPDB=ligandspdb)
                        