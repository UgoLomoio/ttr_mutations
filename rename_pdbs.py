#changing nomenclature from historical to standard: see https://www.ncbi.nlm.nih.gov/books/NBK1194/ 

"""
Historical protein numbering was based on the mature protein after cleavage of a 20-amino-acid signal sequence (e.g., p.Val50Met would be referred to as Val30Met)
We will use the standard mutation numbering by adding 20 to every aminoacid number formatted using historical numbering.
"""

import os 

cwd = os.getcwd()
sep = os.sep

pdb_path = cwd + sep + "pdbs-alphafold"

complex_types = ["tetramer", "dimer", "monomer"]

if __name__ == "__main__":

    for subdir in complex_types:
        
        print(subdir)
        subpath = pdb_path + sep + subdir
        if not os.path.exists(subpath):
            os.makedirs(subpath) 

        pdb_files = [filename for filename in os.listdir(subpath) if filename.split(".")[-1] == "pdb"]
        for filename in pdb_files:
            if "_" in filename:
                filename = filename.split(".")[0]
                if filename.split("_")[1].lower() == "wt": 
                    multimer = filename.split("_")[2] if filename.split("_")[2] in complex_types else None 
                    new_filename =  "{}-{}.pdb".format("wt", multimer) if multimer is not None else "{}.pdb".format("wt")
                else:
                    mut = filename.split("_")[2] 
                    or_aa = mut[0]
                    mut_number = int(mut[1:-1])
                    new_mut_number = str(mut_number + 20)
                    mut_aa = mut[-1]
                    new_mut = or_aa + str(new_mut_number) + mut_aa 
                    multimer = filename.split("_")[3] if filename.split("_")[3] in complex_types else None 
                    new_filename =  "{}-{}.pdb".format(new_mut, multimer) if multimer is not None else "{}.pdb".format(new_mut)
                    directory = multimer if multimer is not None else "monomer"
                    new_filepath = pdb_path + sep + directory + sep + new_filename
                    if not os.path.exists(pdb_path + sep + directory):
                        os.makedirs(pdb_path + sep + directory)
                print("{} -> {}".format(filename, new_filename))
  
                filepath = subpath + sep + filename + ".pdb"
                new_filepath = subpath + sep + new_filename
                os.rename(filepath, new_filepath)
       