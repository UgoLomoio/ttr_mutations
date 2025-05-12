import requests
from rdkit import Chem
import os
import json 
import pandas as pd 
import ast 

url = "https://health.api.nvidia.com/v1/biology/mit/diffdock"
header_auth = ""
if header_auth == "":
    raise Exception("Missing NVIDIA NIM APIKEY in 'header_auth' variable")

cwd = os.getcwd()
sep = os.sep 

docking_path = cwd + sep + "docking"
variants_path = cwd + sep + "pdbs-alphafold" + sep + "tetramer"
output_docking = docking_path + sep + "docked-outputs"
sdf_existing_path = docking_path + sep + "sdf-existing"

sdf_files = os.listdir(sdf_existing_path)
variant_files = os.listdir(variants_path)

def read_pdb(file_path):
    with open(file_path, "r") as file:
        file_data = file.readlines()

    # Join and format the data (you may want to keep the `\n` for line breaks)
    formatted_data = "".join(line.rstrip() + "\n" for line in file_data)

    #print(f"Preview of {file_path}:\n", formatted_data[-50:])
    return str(formatted_data)


def read_sdf(file_path):
    suppl = Chem.SDMolSupplier(file_path)
    molecules = [mol for mol in suppl if mol is not None]

    if not molecules:
        print("Error!")
        return ""


    sdf_text = Chem.MolToMolBlock(molecules[0])

    #print(f"Preview of {file_path}:\n", sdf_text[-50:])
    return str(sdf_text)


def _upload_asset(input):
    assets_url = "https://api.nvcf.nvidia.com/v2/nvcf/assets"

    headers = {
        "Authorization": header_auth,
        "Content-Type": "application/json",
        "accept": "application/json",
    }

    s3_headers = {
        "x-amz-meta-nvcf-asset-description": "diffdock-file",
        "content-type": "text/plain",
    }

    payload = {
        "contentType": "text/plain",
        "description": "diffdock-file"
        
    }

    response = requests.post(
        assets_url, headers=headers, json=payload, timeout=30
    )

    response.raise_for_status()

    asset_url = response.json()["uploadUrl"]
    asset_id = response.json()["assetId"]

    response = requests.put(
        asset_url,
        data=input,
        headers=s3_headers,
        timeout=300,
    )

    response.raise_for_status()
    return asset_id


def save_json_response(protein_id, ligand_id, output_path):
    headers = {
        "Content-Type": "application/json",
        "NVCF-INPUT-ASSET-REFERENCES": ",".join([protein_id, ligand_id]),
        "Authorization": header_auth
    }

    r = requests.post(url, headers=headers, json={
        "ligand": ligand_id,
        "ligand_file_type": "sdf",  # Adjust if needed
        "protein": protein_id,
        "num_poses": 20,
        "time_divisions": 20,
        "steps": 18,
        "save_trajectory": True,
        "is_staged": True
    })

    if r.status_code == 200:
        try:
            data = r.json()

            # Save only if the job was successful
            if data.get("status") == "failed":
                print("Docking job failed. JSON not saved.")
                return None

            with open(output_path, 'w') as json_file:
                json.dump(data, json_file, indent=4)

            print(f"JSON response saved to: {output_path}")
            return output_path  # Return file path for later use

        except json.JSONDecodeError:
            print("Error: Response is not valid JSON.")
    else:
        print(f"Request failed. Status code: {r.status_code}")

    return None

def read_diffdock_local(filepath):
    with open(filepath, "r") as file:
        content = file.read()  # Read entire file as a string

    # Convert JSON string to dictionary
    dictionary = ast.literal_eval(content)

    # Extract values
    confidence = dictionary.get("position_confidence", [])  # Use .get() to avoid errors
    trajectories = dictionary.get("trajectory", [])
    ligand_positions = dictionary.get("ligand_positions", [])

    return confidence, trajectories, ligand_positions

def save_confidence_scores(confidence_scores, variants, output_csv=None):
    """
    Saves a list of confidence scores (each list has 20 values) into a pandas DataFrame.
    """
    # Create column names: rank_1, rank_2, ..., rank_20
    num_poses = len(confidence_scores[0])  # Assuming all lists have 20 values
    columns =  ["Mutation"]
    [columns.append(f"rank_{i+1}") for i in range(num_poses)]

    rows = []
    for i, confidence_score in enumerate(confidence_scores):
        row = []
        variant = variants[i]
        row.append(variant)
        for c in confidence_score:
            row.append(c)
        rows.append(row)
    # Convert to DataFrame
    df = pd.DataFrame(rows, columns=columns)
    # Save to CSV if output path is provided
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"Data saved to {output_csv}")
    return df


if __name__ == "__main__":
            
    confidence_scores = []
    docked_results = os.listdir(output_docking)

    variants = [variant.split(".")[0].split("-")[0] for variant in variant_files]
    for sdf in sdf_files:
        if "acoramidis" in sdf:
            continue #skip, already done
        for variant in variant_files:
            filename = "{}_{}.json".format(variant.split(".")[0], sdf.split(".")[0])
            if filename in docked_results: #if docking results already retrived: read them locally, no api call needed
                print("Found local docking file for {}".format(filename))
                confidence, _, _ = read_diffdock_local(output_docking + sep + filename)
            else: #else: do api call
                print("API call for {}".format(filename))
                protein_id = read_pdb(variants_path + sep + variant)
                protein_id = _upload_asset(protein_id)
                ligand_id = read_sdf(sdf_existing_path + sep + sdf)
                ligand_id = _upload_asset(ligand_id)
                filepath = save_json_response(protein_id, ligand_id, output_docking + sep + filename)
                if filepath is not None:
                    confidence, _, _ = read_diffdock_local(filepath)
            confidence_scores.append(confidence)

        df = save_confidence_scores(confidence_scores, variants, output_csv=output_docking + sep + "diffdock_{}_existing.csv".format(sdf.split(".")[0]))
    #print(df)
