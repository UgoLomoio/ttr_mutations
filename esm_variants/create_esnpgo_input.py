import pandas as pd

header = ">P02766"
fasta_sequence = "MASHRLLLLCLAGLVFVSEAGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE"

df = pd.read_csv("embeddings_umap.csv")
mutations = df["Variant"].values
mutations_str = ""
for i, mutation in enumerate(mutations):
    if mutation != "wt":
        aa_original = mutation[0]
        aa_mutated = mutation[-1]
        aa_index = mutation[1:-1]
        aa_index_hugo = int(aa_index) + 20
        mutation_hugo = f"{aa_original}{aa_index_hugo}{aa_mutated}"
    else:
        continue 
    if i == len(mutations) - 1:
        mutations_str += f"{mutation_hugo}"
    else:
        mutations_str += f"{mutation_hugo},"
header = f">{header} {mutations_str}"
print(header)

with open("esnpgo_input.fasta", "w") as fasta_file: 
    fasta_file.write(f"{header}\n")
    fasta_file.write(f"{fasta_sequence}\n")