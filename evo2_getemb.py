import os 
import requests
import json
import base64 
import numpy as np 

cwd = os.getcwd()
sep = os.sep

embeddings_path = cwd + sep + 'evo2_embeddings'
fasta_sequences_path = cwd + sep + 'nucleotide_sequences'

if not os.path.exists(embeddings_path): 
    os.makedirs(embeddings_path)

# Get the embeddings

key = os.getenv("NGC_API_KEY") or input("Paste the Run Key: ")


def read_sequence(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
        sequence = lines[0]
    return sequence

if __name__ == "__main__":

    # Read the sequence from the fasta file
    for filename in os.listdir(fasta_sequences_path):

        print(filename)
        filepath = fasta_sequences_path + sep + filename
        filename = filename.split(".")[0]
        sequence = read_sequence(filepath)
        print(sequence)
        url="http://localhost:8000/biology/arc/evo2/forward"
        
        data = {
            "sequence": sequence,
            "output_layers": ["embedding_layer","unembed"]
        }
        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {key}",
            "Accept": "application/json"  # API returns Base64-encoded NPZ
        }

        # Send POST Request
        response = requests.post(url, headers=headers, json=data)

        # Check for errors
        if response.status_code == 200:
            response_json = response.json()
    
            if "data" in response_json:  # Ensure "data" key exists
                base64_str = response_json["data"]  # Base64 string
                npz_bytes = base64.b64decode(base64_str)  # Decode Base64
        
                # Save as NPZ file
                with open("output.npz", "wb") as f:
                    f.write(npz_bytes)
        
                print("NPZ file saved as output.npz")

                # Load NPZ file
                npz_data = np.load("output.npz")
                print("NPZ Contents:", npz_data.files)  # List of stored arrays
            else:
                print("Error: 'data' key not found in response.")
        else:
            print(f"Error {response.status_code}: {response.text}")

        break