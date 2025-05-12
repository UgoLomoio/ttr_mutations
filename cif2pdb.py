import os 
import pymol2
import shutil

cwd = os.getcwd()
sep = os.sep


pdb_dir = cwd + sep + "pdbs-alphafold"
cif_dir = cwd + sep + "cifs-alphafold"
cif_subdirs = [elem for elem in os.listdir(cif_dir) if os.path.isdir(cif_dir + sep + elem)]

def get_first_cif(path):

    cifs = [file for file in os.listdir(path) if file.split(".")[-1] == "cif"]
    cif = cifs[0]
    return cif

if __name__ == "__main__":


    for subdir in cif_subdirs:

        print(subdir)
        if subdir == "multimer":
            subdir == "tetramer"

        subpath = cif_dir + sep + subdir
        cif_file = get_first_cif(subpath)    
        cif_path = subpath + sep + cif_file
        print(cif_path)

        mutated = True if subdir.split("_")[0] == "mut" else False
        if mutated:
            _, mut, complex_type = subdir.split("_")
        else:
            if "multimer" in subdir:
                _, _, complex_type = subdir.split("_")
            else:
                _, _ = subdir.split("_")
                complex_type = "monomer"
                
        if complex_type == "multimer":
            complex_type = "tetramer"

        pdb_file = cif_file.replace(".cif", ".pdb").replace("multimer", "tetramer")
        os.makedirs(pdb_dir + sep + subdir, exist_ok=True)
        pdb_path = pdb_dir + sep + subdir + sep + pdb_file
 
        print(cif_path, pdb_path)
        
        with pymol2.PyMOL() as pymol:
            pymol.cmd.load(cif_path,'myprotein')
            pymol.cmd.save(pdb_path, selection='myprotein')
        
        complex_path = pdb_dir + sep + complex_type
        os.makedirs(complex_path, exist_ok=True)
        new_pdb_path = complex_path + sep + pdb_file
        print(complex_path, new_pdb_path)
        
        shutil.move(pdb_path, new_pdb_path)
        os.rmdir(pdb_dir + sep + subdir)