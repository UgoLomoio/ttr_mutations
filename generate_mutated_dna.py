# Importing necessary library
from pathlib import Path
import os 
import pandas as pd 

cwd = os.getcwd()
sep = os.sep 


def process_delins(sequence_wt, mutation):
        id_del, temp = mutation.split("_")  # e.g., "220", "221delinsCT"
        id_del = int(id_del) - 1  # Convert to 0-based index
        id_ins, temp = temp.split("delins")  # e.g., "221", "CT"
        id_ins = int(id_ins) - 1  # Convert to 0-based index
        
        # Extract deletion and insertion residues
        del_res, ins_res = temp[0], temp[1]
        
        # Ensure indices are within the bounds of the sequence
        if id_ins < 0 or id_ins > len(sequence_wt):
            raise IndexError("Insertion index {} is out of bounds".format(id_ins))
        if id_del < 0 or id_del >= len(sequence_wt):
            raise IndexError("Deletion index {} is out of bounds".format(id_del))

        # Print initial sequence and mutation details for debugging
        print("Initial sequence:", sequence_wt)
        print("Attempting to insert '{}' at position {}".format(ins_res, id_ins))
        
        # Insert the nucleotide
        sequence = insert_mut(sequence_wt, ins_res, id_ins)
        print("Sequence after insertion:", sequence)

        # Check the nucleotide for deletion
        if sequence[id_del] == del_res:
            sequence = sequence[:id_del] + sequence[id_del + 1:]  # Delete the nucleotide
        else:
            raise Exception("Nucleotide in position {} is not {} (found {})".format(id_del + 1, del_res, sequence[id_del]))  # Adjust for human-readable index

        return sequence  # Return the modified sequence


def read_nucleotides(wt_file):

    nucleotides = ""
    with open(wt_file) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if i != 0:
            for char in line:
                nucleotides += char
    return nucleotides

def read_list_from_file(mutationlist_path):

    df = pd.read_excel(mutationlist_path)
    mutations = df["Nucleotide"]
    mutations_list = []
    for mutation in mutations:
        if "," in mutation:
            temp_muts = mutation.split(",")
            for temp_mut in temp_muts:
                mutations_list.append(temp_mut)
        else:
            mutations_list.append(mutation)
    return mutations_list

def insert_mut(string, mut, index):
    return string[:index] + mut + string[index:]

def generate_mutated(sequence_wt, mutations_list, output_directory):

    for mutation in mutations_list:

        print("Mutation: ", mutation)

        if "delins" in mutation: #220_221delinsCT

            sequence = process_delins(sequence_wt, mutation)

        elif "-" in mutation: # 265-266 TA>AT
            ids, muts = mutation.split(" ")
            id0, id1 = ids.split("-")
            id0 = int(id0) - 1
            id1 = int(id1) - 1
            org, mut = muts.split(">")
            org0, org1 = org
            mut0, mut1 = mut
            if sequence_wt[id0] == org0:
                sequence = (
                    sequence_wt[:id0] + mut0 + sequence_wt[id0:]
                )         
            else:
                raise Exception ("Nucleotide in position {} is not {} (found {})".format(id0, org0, sequence_wt[id0]))
            if sequence_wt[id1] == org1:
                sequence = (
                    sequence_wt[:id1] + mut1 + sequence_wt[id1:]
                )        
            else:
                raise Exception ("Nucleotide in position {} is not {} (found {})".format(id1, org1, sequence_wt[id1]))
            
        else: #idO>M
            temp, mut = mutation.split(">")
            org = temp[-1] 
            id = int(temp[:-1]) - 1
            if sequence_wt[id] == org:
                sequence = (
                    sequence_wt[:id] + mut + sequence_wt[id:]
                )
        
            else:
                raise Exception ("Nucleotide in position {} is not {} (found {})".format(id, org, sequence_wt[id]))
            
        with open(output_directory + sep + "{}.fna".format(mutation), "w") as f:
            f.writelines(sequence)

if __name__ == "__main__":

    # Wild-type sequence
    wt_file = cwd + sep + "gene.fna"
    sequence_wt = read_nucleotides(wt_file)
    mutations_list = read_list_from_file(cwd + sep + "aminoacid_nucleotide_mutations.xlsx")

    #Output directory for FASTA files
    output_directory = cwd + sep + "nucleotide_sequences"
    os.makedirs(output_directory, exist_ok=True)
    generate_mutated(sequence_wt, mutations_list, output_directory)