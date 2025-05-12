import os 
import shutil

cwd = os.getcwd()
sep = os.sep 

cifs_dir = cwd + sep + "cifs-alphafold"
af_output_dir = cwd + sep + "structures-alphafold"
af_subdirs = [elem for elem in os.listdir(af_output_dir) if os.path.isdir(af_output_dir + sep + elem)]
print(af_subdirs)

def get_first_cif(path):

    cifs = [file for file in os.listdir(path) if file.split(".")[-1] == "cif"]
    for cif in cifs:
        temp = cif.split(".")[0]
        temp = temp.split("_")
        if temp[-1] == "0":
            return cif
    return None 

if __name__ == "__main__":

    for subdir in af_subdirs:

        print(subdir)
        subpath = af_output_dir + sep + subdir
        cif_file = get_first_cif(subpath)    
        cif_path = subpath + sep + cif_file 
        print(cif_path)

        mutated = True if subdir.split("_")[0] == "mut" else False
        if mutated:
            if len(subdir.split("_")) == 2:
                _, mut = subdir.split("_")
                complex_type = "monomer"
            else:
                _, mut, complex_type = subdir.split("_")
        else:
            if "multimer" in subdir:
                _, _, complex_type = subdir.split("_")
            else:
                _, complex_type = subdir.split("_")
        if complex_type == "monomer":
            os.makedirs(cifs_dir + sep + subdir+"_monomer" , exist_ok=True)
            filename_new = cifs_dir + sep + subdir+"_monomer" + sep + cif_file
        else:
            os.makedirs(cifs_dir + sep + subdir, exist_ok=True)
            filename_new = cifs_dir + sep + subdir + sep + cif_file
        print(cif_path, filename_new)
        shutil.copyfile(cif_path, filename_new)