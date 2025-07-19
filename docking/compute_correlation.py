import pandas as pd
import os
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import plotly.io as pio
import numpy as np 

pio.kaleido.scope.mathjax = None

plt.rcParams['font.size'] = 20  # General font size
plt.rcParams['font.weight'] = 'bold'  # All text bold
plt.rcParams['axes.titlesize'] = 22  # Title font size
plt.rcParams['axes.labelsize'] = 16  # X-Y label font size
plt.rcParams['legend.fontsize'] = 16  # Legend font size
plt.rcParams['xtick.labelsize'] = 12  # X-tick font size
plt.rcParams['ytick.labelsize'] = 12  # Y-tick font size

cwd = os.getcwd()
sep = os.sep 

docking_path = cwd + sep + "docking"
results_path = cwd + sep + "results"

sdf_type = "existing"
dockers = ["DiffDock", "AutoDock"] #DiffDock or AutoDock

sdf_path = f"sdf-{sdf_type}"
ligand_dir = docking_path + sep + sdf_path
ligands = ["tafamidis", "acoramidis"]

output_dir = os.path.join(cwd, "docking", "converted_pdbs")
docked_outputs = os.path.join(cwd, "docking", "AutoDock", "docked-outputs", "existing")  # Adjust path as needed
ligand_reference = os.path.join(output_dir, 'reference_ligand.pdb')


def agreement_score(dfs_ligands):

    mutations = []
    for ligand, dfs in dfs_ligands.items():
        for docker, df in dfs.items():
            for mutation in df["Mutation"].tolist():
                if mutation not in mutations:
                    mutations.append(mutation)

    ligands = list(dfs_ligands.keys())
    df_agreement = pd.DataFrame(columns=["Mutation", "Prediction_DiffDock", "Prediction_AutoDock", "Agreement"])
    for mutation in mutations:
        ligands_scores = {}
        for ligand in ligands:
            filter = dfs_ligands[ligand]["DiffDock"]["Mutation"] == mutation
            if filter.sum() == 0:
                print(f"Mutation {mutation} not found in {ligand} for DiffDock")
                continue
            # Get the corresponding values from both DataFrames
            vina_score = dfs_ligands[ligand]["AutoDock"][filter]["Vina Score"].values
            confidence = dfs_ligands[ligand]["DiffDock"][filter]["Confidence"].values
            ligands_scores[ligand] = {"Vina Score": vina_score[0], "Confidence": confidence[0]}
        
        # Check if both ligands are present
        if "tafamidis" not in ligands_scores or "acoramidis" not in ligands_scores:
            #print(f"Missing ligand scores for mutation {mutation}. Skipping...")
            continue
        # Find the best ligand
        vina_score_tafamidis = ligands_scores["tafamidis"]["Vina Score"]
        confidence_tafamidis = ligands_scores["tafamidis"]["Confidence"]
        vina_score_acoramidis = ligands_scores["acoramidis"]["Vina Score"]
        confidence_acoramidis = ligands_scores["acoramidis"]["Confidence"]
        if vina_score_tafamidis < vina_score_acoramidis:
            best_ligand_auto = "tafamidis"
        else:
            best_ligand_auto = "acoramidis"
        
        if confidence_tafamidis < confidence_acoramidis:
            best_ligand_diff = "tafamidis"
        else:
            best_ligand_diff = "acoramidis"
        
        df_agreement_temp = pd.DataFrame({
            "Mutation": [mutation],
            "Prediction_DiffDock": [best_ligand_diff],
            "Prediction_AutoDock": [best_ligand_auto],
            "Agreement": [best_ligand_diff == best_ligand_auto]
        })
        df_agreement = pd.concat([df_agreement, df_agreement_temp], ignore_index=True)

    agreement_score = df_agreement["Agreement"].sum() / len(df_agreement) * 100
    print(f"Agreement score: {agreement_score:.2f}%")

    return df_agreement

def create_dataframe(docking_csv, docker):
    
    df = pd.read_csv(docking_csv)
    df["Mutation"] = [mut.split("-")[0] for mut in df["Mutation"]]

    # Separate "wt" and other mutations
    df_wt = df[df["Mutation"] == "wt"]
    df_mut = df[df["Mutation"] != "wt"]

    # Extract numeric part and sort
    df_mut = df_mut.assign(
                num_part=df_mut["Mutation"].str.extract("(\d+)").astype(float)
    ).sort_values("num_part")

    # Drop the helper column and concatenate
    df = pd.concat([df_wt, df_mut.drop(columns=["num_part"])], ignore_index=True)

    # Melt the DataFrame for visualization
    df_melted = df.melt(id_vars=["Mutation"], var_name="Rank Type", value_name="Rank Value")

    if docker == "AutoDock":
        # Calculate the IQR
        Q1 = df_melted['Rank Value'].quantile(0.25)
        Q3 = df_melted['Rank Value'].quantile(0.75)
        IQR = Q3 - Q1

        # Define outlier condition (values outside 1.5 * IQR)
        outlier_condition = (df_melted['Rank Value'] < (Q1 - 1.5 * IQR)) | (df_melted['Rank Value'] > (Q3 + 1.5 * IQR))

        # Remove outliers
        df_melted = df_melted[~outlier_condition]

    return df_melted

def compute_correlation(diffdock_scores, autodock_scores):

    from scipy.stats import pearsonr
    
    print(len(diffdock_scores), len(autodock_scores))
    # Calculate Pearson correlation
    correlation, p_value = pearsonr(autodock_scores, diffdock_scores)
    print(f"Pearson correlation: {correlation:.4f}, p-value: {p_value:.4e}")
    return correlation, p_value

reference_dicts = {ligand: "V50M" for ligand in ligands}

if __name__ == "__main__":

    dfs_ligands = {}
    for ligand in ligands:
        ligand = ligand.split(".")[0]
        dfs = {}
        for docker in dockers:
                
            if docker == "DiffDock":
                docking_csv = docking_path + sep + "docked-outputs" + sep + sdf_type + sep + f"diffdock_{ligand}_existing.csv"
            else:
                docking_csv = docking_path + sep + "AutoDock" + sep + "docked-outputs" + sep + sdf_type + sep + f"vina_score_docking_{ligand}_TTRtetramer.csv"

            df_melted = create_dataframe(docking_csv, docker)
            df_melted["Mutation"] = [mut.upper() for mut in df_melted["Mutation"]]
                
            if docker == "AutoDock":

                to_df = []
                idx = 0
                score = df_melted[df_melted["Rank Type"] == f"Score_{str(idx+1)}"]
                df = df_melted.copy()
                df = df[df["Rank Type"] == f"Score_{str(idx+1)}"]
                df = df.drop("Rank Type", axis = 1)
                df.rename(columns={'Rank Value': 'Vina Score'}, inplace=True)
                df["Mutation"] = [mutant.upper() for mutant in df["Mutation"]]

            else:
                df = df_melted.copy()
                df = df[df["Rank Type"] == "rank_1"]
                df = df.drop("Rank Type", axis = 1)
                df.rename(columns={'Rank Value': 'Confidence'}, inplace=True)

            dfs[docker] = df

        mutants_autodock = dfs["AutoDock"]["Mutation"].tolist()
        mutants_diffdock = dfs["DiffDock"]["Mutation"].tolist()
        mutants_missing_autodock = set(mutants_diffdock) - set(mutants_autodock)
        mutants_missing_diffdock = set(mutants_autodock) - set(mutants_diffdock)
        print(f"Mutants missing in AutoDock for {ligand}: {mutants_missing_autodock}")
        print(f"Mutants missing in DiffDock for {ligand}: {mutants_missing_diffdock}")

        # Removing the missing mutants from the DataFrames
        dfs["AutoDock"] = dfs["AutoDock"][~dfs["AutoDock"]["Mutation"].isin(mutants_missing_diffdock)]
        dfs["DiffDock"] = dfs["DiffDock"][~dfs["DiffDock"]["Mutation"].isin(mutants_missing_autodock)]

        correlation, p_value = compute_correlation(dfs["AutoDock"]["Vina Score"], dfs["DiffDock"]["Confidence"])
        print(f"Correlation between AutoDock and DiffDock for {ligand}: {correlation:.4f}, p-value: {p_value:.4e}")

        dfs_ligands[ligand] = dfs

    agreement_score_df = agreement_score(dfs_ligands)
    print(agreement_score_df)
