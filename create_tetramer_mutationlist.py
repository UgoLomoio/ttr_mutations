import os 

cwd = os.getcwd()
sep = os.sep

if __name__ == "__main__":
    
    input_file_path = cwd + sep + "ttrmutations.txt"
    output_file_path = cwd + sep + "ttrmutations_tetramer.txt"
    output_file_path_dimer = cwd + sep + "ttrmutations_dimer.txt"

    # Read the input file and process the mutations
    with open(input_file_path, "r") as infile:
        mutations = infile.readlines()

    # Reformat the output so that each line contains the same mutation for chains A, B, C, and D
    reformatted_mutations = []

    for mutation in mutations:
        mutation = mutation.strip()  # Remove any leading/trailing whitespace
        chain, mut = mutation.split(" ", 1)  # Split into chain and mutation
        reformatted_mutations.append(f"A {mut};B {mut};C {mut};D {mut}")

    # Write the reformatted mutations to a new file
    with open(output_file_path, "w") as outfile:
        outfile.write("\n".join(reformatted_mutations) + "\n")
        
    # Reformat the output so that each line contains the same mutation for chains A, B
    reformatted_mutations = []

    for mutation in mutations:
        mutation = mutation.strip()  # Remove any leading/trailing whitespace
        chain, mut = mutation.split(" ", 1)  # Split into chain and mutation
        reformatted_mutations.append(f"A {mut};B {mut}")

    # Write the reformatted mutations to a new file
    with open(output_file_path_dimer, "w") as outfile:
        outfile.write("\n".join(reformatted_mutations) + "\n")



