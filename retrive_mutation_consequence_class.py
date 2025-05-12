import os 
import requests
import pandas as pd 
import re

cwd = os.getcwd()
sep = os.sep 

input_file_path = cwd + sep + "ttrmutations.txt"
# Read the input file and process the mutations
with open(input_file_path, "r") as infile:
    mutations = infile.readlines()
mutations = [mutation.split(" ")[1].strip() for mutation in mutations]

def get_ttr_mutation_consequences(mutations: list):
    """
    Fetches the consequence class for multiple TTR mutations using the UniProt API.
    
    :param mutations: List of mutations in the forma "V30M" (Protein-level HGVS notation)
    :return: Dictionary mapping mutations to their consequence class
    """

    uniprot_id = "P02766"  # UniProt ID for human TTR
    url = "https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=100&location={}&accession={}"    
    mutation_consequences = {}
    
    for mutation in mutations:
        mutation_consequences[mutation] = {}
        wt = mutation[0]
        mut = mutation[-1]
        location = int(mutation[1:-1]) + 20
        print(wt, location, mut)
        response = requests.get(url.format(location, uniprot_id))
        if response.status_code != 200:
            print("Error fetching data from UniProt")
            return None
        data = response.json()[0]
        features = data.get("features")
        found = False
        for feature in features:
            # Extract "mutatedType""
            mutatedType = feature.get("mutatedType")
            print("{}->{}".format(mut, mutatedType))
            if mut == mutatedType:
                print("found")
                # Extract "clinicalSignificances" -> "type" if exists, otherwise "Unknown"
                clinical_type = feature.get("clinicalSignificances", {})
                if len(clinical_type)>0:
                    clinical_type = clinical_type[0].get("type", "Unknown")
                else:
                    clinical_type = "Unknown"
                mutation_consequences[mutation]["Clinical Significance"] = clinical_type
                found = True

        if not found:
            print("not found")
            mutation_consequences[mutation]["Clinical Significance"] = "Not found"

    return mutation_consequences

if __name__ == "__main__":

    consequences = get_ttr_mutation_consequences(mutations)
    for mutation, consequence in consequences.items():
        print(f"Consequence class for {mutation}: {consequence}")
    df = pd.DataFrame.from_dict(consequences, orient="index")
    df.to_csv(cwd + sep + "mutations_consequences.csv")