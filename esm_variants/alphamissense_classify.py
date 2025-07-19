#Starting with alphamissense outputs compute accuracy 

import os 
import pandas as pd
import numpy as np
from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score

def use_hgcn(df):
    df_copy = df.copy()
    for mutation in df['Mutation']:
        aa_original = mutation[0]
        aa_mutated = mutation[-1]
        aa_index = mutation[1:-1]
        aa_index_hugo = int(aa_index) + 20
        mutation_hugo = f"{aa_original}{aa_index_hugo}{aa_mutated}"
        #print(f"Processing mutation: {mutation}, HUGO mutation: {mutation_hugo}")
        df_copy.loc[df['Mutation'] == mutation, 'Mutation'] = mutation_hugo
    return df_copy

cwd = os.getcwd()
sep = os.sep

parent_dir = os.path.dirname(cwd)
results_path = parent_dir + sep + "results"

df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["Mutation", "Consequence"]
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely benign': 'Benign'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic'})
df_consequences = use_hgcn(df_consequences)

#alphamissense_output = cwd + sep + "AlphaMissense-Hotspot-P02766" + sep + "AlphaMissense-Hotspot-P02766.tsv"
alphamissense_output = cwd + sep + "alphamissense_P02766.tsv"
alphamissense_df = pd.read_csv(alphamissense_output, sep = "\t")

map_resiname = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

map_resiname_inv = {v: k for k, v in map_resiname.items()}


def three2one(mutation):
    """
    Convert a mutation from three-letter amino acid code to one-letter code.
    Example: 'p.Ala123Gly' -> 'A123G'
    """
    mutation = mutation.split(".")[1]  # Extract the mutation part after 'p.'
    if len(mutation) < 4:
        return mutation  # Return as is if the mutation is too short
    aa_original = mutation[:3]
    position = mutation[3:-3]
    aa_mutated = mutation[-3:]
    
    # Convert three-letter code to one-letter code
    aa_original_one = map_resiname_inv.get(aa_original, aa_original)
    aa_mutated_one = map_resiname_inv.get(aa_mutated, aa_mutated)
    
    return f"{aa_original_one}{position}{aa_mutated_one}"

#alphamissense_df = alphamissense_df.filter(['a.a.', 'position', 'benign variants', 'pathogenic variants', 'ambiguous variants'])
alphamissense_df = alphamissense_df.filter(['protein variant', 'pathogenicity score', 'pathogenicity class'])
alphamissense_df = alphamissense_df.replace({'pathogenicity class': {'likely_benign': 'Benign', 'likely_pathogenic': 'Pathogenic'}})
alphamissense_df = alphamissense_df.rename(columns={'protein variant': 'Mutation', 'pathogenicity score': 'Score', 'pathogenicity class': 'Prediction'})
alphamissense_df['Mutation'] = [three2one(mutation) for mutation in alphamissense_df['Mutation']]

def create_final_df(alphamissense_df):
    
    rows = [] 
    for index, row in alphamissense_df.iterrows():
        position = row['position']
        benign_variants = row['benign variants']
        if pd.isna(benign_variants):
            benign_variants = []
        else: 
            benign_variants = benign_variants.split(':')[1].split(',')
    
        pathogenic_variants = row['pathogenic variants']
        if pd.isna(pathogenic_variants):
            pathogenic_variants = []
        else:
            pathogenic_variants = pathogenic_variants.split(':')[1].split(',')

        aa = row['a.a.']
        for variants in benign_variants:
            for aa_mut in variants:
                mutation = f"{aa}{position}{aa_mut}"
                rows.append([mutation, 'Benign'])
        
        for variants in pathogenic_variants:
            for aa_mut in variants:
                mutation = f"{aa}{position}{aa_mut}"
                rows.append([mutation, 'Pathogenic'])
    
    df = pd.DataFrame(rows, columns=['Mutation', 'Prediction'])
    df = df.drop_duplicates(subset=['Mutation'])
    return df 

if __name__ == "__main__":

    df_consequences_unknown = df_consequences[df_consequences['Consequence'] == 'Unknown']
    df_consequences = df_consequences[df_consequences['Consequence'] != 'Unknown']
    
    #final_df = create_final_df(alphamissense_df)
    final_df = alphamissense_df
    final_df = final_df[final_df['Mutation'].isin(df_consequences['Mutation'])]

    y_pred = final_df['Prediction'].tolist()    
    y_true = df_consequences['Consequence'].tolist()

    y_pred = [1 if pred == 'Pathogenic' else 0 for pred in y_pred]
    y_true = [1 if true == 'Pathogenic' else 0 for true in y_true]
    map_class = {0: 'Benign', 1: 'Pathogenic'}

    print(f"Number of predictions: {len(y_pred)}")
    print(f"Number of true labels: {len(y_true)}")

    rocauc = roc_auc_score(y_true, y_pred)
    print(f"ROC AUC Score: {rocauc}")
    
    conf_matrix = confusion_matrix(y_true, y_pred)
    print("Confusion Matrix:")
    print(conf_matrix)
    
    class_report = classification_report(y_true, y_pred)
    print("Classification Report:")
    print(class_report)

    idx_missclassified = [i for i, (pred, true) in enumerate(zip(y_pred, y_true)) if pred != true]
    print(f"Number of misclassified samples: {len(idx_missclassified)}")
    print("Misclassified samples:") 
    for idx in idx_missclassified:
        print(f"Predicted: {map_class[y_pred[idx]]}, Target: {map_class[y_true[idx]]}, Mutation: {final_df['Mutation'].iloc[idx]}")
    