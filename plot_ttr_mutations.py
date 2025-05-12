from pymol import cmd
import os 

cwd = os.getcwd()
sep = os.sep

input_file_path = cwd + sep + "ttrmutations.txt"
# Read the input file and process the mutations
with open(input_file_path, "r") as infile:
    mutations = infile.readlines()

mutations = [mutation.split(" ")[1][1:-2] for mutation in mutations]
chains = ["A", "B", "C", "D"]

complex_types = ["monomer", "dimer", "tetramer"]
pdbs_path = cwd + sep + "pdbs-alphafold"
pdb_files = {
    complex_type: os.listdir(pdbs_path + sep + complex_type) for complex_type in complex_types
}

pdb_files = {complex_type: [filename for filename in pdb_files[complex_type] if "wt" in filename] for complex_type in complex_types}

mutations_str = ""
for i, mutation in enumerate(mutations):
    if i == 0:
        mutations_str = mutation
    else:
        mutations_str += "+{}".format(mutation)
print(mutations_str)

for complex_type in complex_types:
    for file in pdb_files[complex_type]:
        cmd.load(pdbs_path + sep + complex_type + sep + file, file)
        cmd.do("set sphere_color, red")
        cmd.do("set sphere_scale, 0.5")
        if complex_type == "monomer":
            chains_to_do = "A"
        elif complex_type == "dimer":
            chains_to_do = "A+B"
        else: #TETRAMER
            chains_to_do = "A+B+C+D"
        line = "show spheres, (resi "+ mutations_str + " and chain "+ chains_to_do + " and name CA)"
        cmd.do(line)
        cmd.save(cwd + sep + "results" + sep + "{}.pse".format(file.split(".")[0]))
        cmd.delete(file)