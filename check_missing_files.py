import os

cwd = os.getcwd()
sep = os.sep 

pdb_path = cwd + sep + "pdbs-alphafold"
mutations_file = cwd + sep + "ttrmutations.txt"
complex_types = ["monomer", "tetramer"]

# Function to read a list from a file
def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        # Read lines and remove any leading/trailing whitespace characters (like newline)
        my_list = [line.strip().replace("A ", "") for line in file.readlines()]
    return my_list


if __name__ == "__main__":
    mutations = read_list_from_file(mutations_file)
    missing_mutations = {}
    for complex_type in complex_types:
        missing_mutations[complex_type] = []
        for mutation in mutations:
            subdir = pdb_path + sep + complex_type + sep
            filenames = os.listdir(subdir)
            found = False
            for filename in filenames:
                if "wt" in filename:
                    continue

                if complex_type == "tetramer":
                    mut_file = filename.split("-")[0].upper()
                else:
                    mut_file = filename.split(".")[0].upper()
                
                resi = int(mut_file[1:-1]) - 20 
                mut_file = mut_file[0] + str(resi) + mut_file[-1]
                #print(mut_file, mutation, len(mut_file), len(mutation))
                if mut_file == mutation:
                    found = True 
                    break
            if not found:
                missing_mutations[complex_type].append(mutation)
        
    for complex_type, missing_files in missing_mutations.items():
        missing_files_str = ""
        for file in missing_files:
            if file != missing_files[-1]:
                missing_files_str+= "{}, ".format(file)
            else:
                missing_files_str+=file
        print("Missing files for {}: {}".format(complex_type, missing_files_str))