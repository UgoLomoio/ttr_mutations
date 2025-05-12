import subprocess
from pathlib import Path

import os
import json

os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"]="true"

def fasta2json(fasta_filepath):
    """
    Reads a FASTA file and writes its content to a JSON file inside the 'json_sequences' directory.

    Args:
        fasta_filepath (str): The path to the input FASTA file.
    """
    # Ensure the 'json_sequences' directory exists
    json_dir = "json_sequences"
    os.makedirs(json_dir, exist_ok=True)

    # Extract the filename without extension
    fasta_filename = os.path.basename(fasta_filepath)
    json_filename = os.path.splitext(fasta_filename)[0] + ".json"

    # Initialize the JSON structure
    json_data = []

    # Read the FASTA file and extract sequences
    try:
        with open(fasta_filepath, "r") as fasta_file:
            sequences = []
            sequence = ""
            for line in fasta_file:
                line = line.strip()
                if line.startswith(">"):
                    if sequence:  # Save the previous sequence
                        sequences.append(sequence)
                        sequence = ""
                else:
                    sequence += line  # Continue building the sequence
            if sequence:  # Add the last sequence
                sequences.append(sequence)

            # Populate JSON data with the extracted sequences
            json_data.append({
                "name": fasta_filename,
                "modelSeeds": ["1589660288"],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": seq,
                            "count": len(sequences)
                        }
                    } for seq in sequences
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


def run_alphafold(fasta_dir, output_dir, alphafold_script, num_copies=1):
    """
    Run AlphaFold3 for each FASTA file.

    Parameters:
    fasta_dir (str): Directory containing FASTA files.
    output_dir (str): Directory to save AlphaFold outputs.
    alphafold_script (str): Path to the AlphaFold3 inference script.
    num_copies (int, optional): Number of protein copies for multimer models.

    Returns:
    None
    """
    # Ensure output directory exists
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    #Directory containing JSON files
    json_directory = "json_sequences"
    
    
    # Loop through all FASTA files in the input directory
    fasta_files = list(Path(fasta_dir).glob("*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        return

    for fasta_file in fasta_files:
    
        json_filepath = fasta2json(fasta_file)
        
        # Define the output directory for this FASTA file
        fasta_output_dir = output_path / fasta_file.stem 
        fasta_output_dir.mkdir(parents=True, exist_ok=True)

        print(f"Running AlphaFold for {fasta_file} ...")

        # create a temporary FASTA file if multimer
        temp_fasta_file = None
        if num_copies > 1:
            with open(fasta_file, "r") as original_fasta:
                sequence = original_fasta.read().strip().split("\n", 1)[1]  # Get the sequence part only

            # Create a temporary multimer FASTA file
            temp_fasta_file = fasta_output_dir / f"{fasta_file.stem}_multimer.fasta"
            with open(temp_fasta_file, "w") as temp_fasta:
                for i in range(num_copies):
                    temp_fasta.write(f">Chain{i+1}\n{sequence}\n")

            fasta_file = temp_fasta_file  # Use the temporary file
            json_filepath = fasta2json(fasta_file)
            
        # Run the AlphaFold inference script
        command = [
            "python", alphafold_script,
            "--json_path", str(json_filepath),
            "--output_dir", str(fasta_output_dir),
        ]

        try:
            subprocess.run(command, check=True)
            print(f"Completed: {fasta_file.stem}")
        except subprocess.CalledProcessError as e:
            print(f"Error running AlphaFold for {fasta_file.stem}): {e}")
        finally:
            # Clean up temporary multimer FASTA file
            if temp_fasta_file and temp_fasta_file.exists():
                temp_fasta_file.unlink()

if __name__ == "__main__":
    
    # Directory containing FASTA files
    fasta_directory = "fasta_sequences"  

    # Output directory for AlphaFold predictions
    output_directory = "alphafold_outputs"  # Alphafold3 output directory

    # Path to the AlphaFold3 inference script
    alphafold_inference_script = "alphafold3/run_alphafold.py"  # Replace with the actual script path

    # Run AlphaFold
    run_alphafold(fasta_directory, output_directory, alphafold_inference_script)
    run_alphafold(fasta_directory, output_directory, alphafold_inference_script, num_copies = 4)
