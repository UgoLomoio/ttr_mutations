# Importing necessary library
from pathlib import Path
import os 

cwd = os.getcwd()
sep = os.sep 

def generate_fasta(wt_sequence, mutations, output_dir):
    """
    Generate FASTA files for a wild-type sequence and its single-mutated variants.

    Parameters:
    wt_sequence (str): The wild-type protein sequence.
    mutations (list): List of mutations in the format ["A50V", "D23N", ...].
    output_dir (str): Directory to save the FASTA files.
    n_copies (int): Number of protein copies.
    Returns:
    None
    """
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Write the wild-type sequence to a FASTA file
    wt_fasta_path = output_path / "WT.fasta"
    with open(wt_fasta_path, "w") as wt_file:
        wt_file.write(f">WT\n{wt_sequence}\n")

    print(f"Wild-type sequence saved to {wt_fasta_path}")

    # Generate single-mutated sequences and save them as FASTA files
    for mutation in mutations:
        original_residue = mutation[0]   # Original amino acid
        position = int(mutation[1:-1])  # Position of the mutation
        new_residue = mutation[-1]      # Mutated amino acid

        # Validate mutation format
        if wt_sequence[position - 1] != original_residue:
            print(f"Warning: Mutation {mutation} does not match the wild-type sequence.")
            continue

        # Create the mutated sequence
        mutated_sequence = (
            wt_sequence[:position - 1] + new_residue + wt_sequence[position:]
        )

        # Save the mutated sequence to a FASTA file
        mutation_fasta_path = output_path / f"Mut_{mutation}.fasta"
        with open(mutation_fasta_path, "w") as mut_file:
            mut_file.write(f">Mut_{mutation}\n{mutated_sequence}\n")

        print(f"Mutated sequence {mutation} saved to {mutation_fasta_path}")

# Function to read a list from a file
def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        # Read lines and remove any leading/trailing whitespace characters (like newline)
        my_list = [line.strip().replace("A ", "") for line in file.readlines()]
    return my_list

if __name__ == "__main__":
    # Wild-type sequence
    sequence_wt = "GPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE"

    mutations_list = read_list_from_file(cwd + sep + "ttrmutations.txt")


    # Output directory for FASTA files
    output_directory = "fasta_sequences"

    generate_fasta(sequence_wt, mutations_list, output_directory)







