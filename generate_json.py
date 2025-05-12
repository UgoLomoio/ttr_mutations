from pathlib import Path
import os
import json

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"]="true"

def fasta2json(fasta_filepath, n_copies = 1):
    """
    Reads a FASTA file and writes its content to a JSON file inside the 'json_sequences' directory.

    Args:
        fasta_filepath (str): The path to the input FASTA file.
    """
    json_dir = "json_jobs"
    os.makedirs(json_dir, exist_ok=True)

    # Extract the filename without extension
    fasta_filename = os.path.basename(fasta_filepath)

    # Initialize the JSON structure
    json_data = []

    # Read the FASTA file and extract sequences
    try:
        with open(fasta_filepath, "r") as fasta_file:
            sequence = fasta_file.read().strip().split("\n", 1)[1]  # Get the sequence part only

            if n > 1:
                filename = fasta_filename.split(".")[0] + "_multimer"
            else:
                filename = fasta_filename.split(".")[0]

            json_filename = os.path.splitext(filename)[0] + ".json"

            # Populate JSON data with the extracted sequences
            json_data.append({
                "name": filename,
                #"modelSeeds": ["1589660288"],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence,
                            "count": n_copies
                        }
                    }
                ],
                "dialect": "alphafoldserver",
                "version": 1
            })

    except FileNotFoundError:
        print(f"Error: File '{fasta_filepath}' not found.")
        return
    except Exception as e:
        print(f"An error occurred: {e}")
        return

    # Write the JSON data to a file
    json_filepath = os.path.join(json_dir, json_filename)
    try:
        with open(json_filepath, "w") as json_file:
            json.dump(json_data, json_file, indent=4)
        print(f"JSON file created: {json_filepath}")
    except Exception as e:
        print(f"Failed to write JSON file: {e}")
        
    return json_filepath

if __name__ == "__main__":
    
    # Directory containing FASTA files
    fasta_directory = "fasta_sequences"  

    #Directory containing JSON files
    json_directory = "json_jobs"    
    
    n_copies = [1, 4]

    # Loop through all FASTA files in the input directory
    fasta_files = list(Path(fasta_directory).glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {fasta_directory}")

    for fasta_file in fasta_files:
    
        json_filepath = fasta2json(fasta_file)
    
        for n in n_copies:
            json_filepath = fasta2json(fasta_file, n_copies = n)
            

    
