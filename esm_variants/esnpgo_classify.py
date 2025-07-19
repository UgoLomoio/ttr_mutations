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

esnpgo_df = pd.read_csv(cwd + sep + "esnpgo_predictions.tsv", sep="\t")

if __name__ == "__main__":


    df_consequences_unknown = df_consequences[df_consequences['Consequence'] == 'Unknown']
    df_consequences = df_consequences[df_consequences['Consequence'] != 'Unknown']
    esnpgo_df = esnpgo_df[esnpgo_df['Variation'].isin(df_consequences['Mutation'])]

    y_pred = esnpgo_df['Pred_class'].tolist()    
    y_true = df_consequences['Consequence'].tolist()

    y_pred = [1 if pred == 'Pathogenic' else 0 for pred in y_pred]
    y_true = [1 if true == 'Pathogenic' else 0 for true in y_true]

    rocauc = roc_auc_score(y_true, y_pred)
    print(f"ROC AUC Score: {rocauc}")
    
    conf_matrix = confusion_matrix(y_true, y_pred)
    print("Confusion Matrix:")
    print(conf_matrix)
    
    class_report = classification_report(y_true, y_pred)
    print("Classification Report:")
    print(class_report)
