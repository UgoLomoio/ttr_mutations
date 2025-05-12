import pandas as pd 
import os
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from sklearn import preprocessing
import re 
import numpy as np 


plt.rcParams['font.size'] = 20  # General font size
plt.rcParams['font.weight'] = 'bold'  # All text bold
plt.rcParams['axes.titlesize'] = 30  # Title font size
plt.rcParams['axes.labelsize'] = 25  # X-Y label font size
plt.rcParams['legend.fontsize'] = 22  # Legend font size
plt.rcParams['xtick.labelsize'] = 18  # X-tick font size
plt.rcParams['ytick.labelsize'] = 18  # Y-tick font size

cwd = os.getcwd()
sep = os.sep

basedir = cwd+sep+"mutant_stability"

filename = basedir + sep + "DynaMut2_Tetramer.xlsx"
results_path = cwd + sep + "results"

# Function to shift the numeric index by 20
def shift_mutation(mutation):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z])", mutation)
    if match:
        prefix, num, suffix = match.groups()
        new_num = str(int(num) + 20)  # Shift by 20
        return f"{prefix}{new_num}{suffix}"
    return mutation  # Return unchanged if pattern doesn't match

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


# Compute consensus classification
def classify_ddg(value):
    if value < -0.01:
        return -1.0 #"Destabilizing"
    elif value > 0.01:
        return +1.0 #"Stabilizing"
    else:
        return 0.0 #"No Effect"

if __name__ == "__main__":

    df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
    df_consequences.columns = ["Mutation", "Consequence"]
    unknown_values = ['Not found', 'Variant of uncertain significance', 'Unknown', 'Likely benign', 'Likely pathogenic']
    # Apply the function to the Mutation column
    df_consequences["Mutation"] = df_consequences["Mutation"].apply(shift_mutation)

    # Load the Excel file
    df = pd.read_excel(filename)
    # Apply the function to the Mutation column
    df["Mutation"] = [mutation.split(";")[0].split(" ")[1] for mutation in df["Mutation"]]
    df["Mutation"] = df["Mutation"].apply(shift_mutation)
    
    df_merged = pd.merge(df, df_consequences, on="Mutation")
    df_merged = df_merged[~df_merged['Consequence'].isin(unknown_values)]


    df_merged.set_index("Mutation", inplace=True)
    # Convert comma decimal format to dot
    for col in df_merged.columns:  
        if col == "Consequence":
            df_merged[col].replace({"Benign": 0.0, "Pathogenic": 1.0}, inplace=True)
        else:
            df_merged[col] = df_merged[col].astype(str).str.replace(',', '.').astype(float)
    
    for col in df_merged.columns:
        if col != "Consequence":
            df_merged[col] = df_merged[col].apply(classify_ddg)
    
    # Selecting only the numerical columns for correlation analysis
    correlation_matrix = df_merged[df_merged.columns].corr()

    # Save the correlation matrix as a CSV file
    correlation_matrix.to_csv(basedir + sep + 'correlation_effect_mutation_matrix.csv')

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5, annot_kws={'size': 14})
    plt.title('Correlation Heatmap between ΔΔG and mutation clinical significance', pad=20)
    silentremove(basedir + sep + 'correlation_ddg_significance_heatmap.png')
    plt.savefig(basedir + sep + 'correlation_ddg_significance_heatmap.pdf', dpi=300)  # Save figure
    plt.show()

    # Plotting the heatmap DDG
    # Create a custom colormap from red to white to green
    # Load the Excel file
    df1 = pd.read_excel(filename, index_col=None)
    df1["Mutation"] = [mutation.split(";")[0].split(" ")[1].strip() for mutation in df1["Mutation"]]
    df1["Mutation"] = df1["Mutation"].apply(shift_mutation)
    df1.drop("Sum ΔΔGStability", axis=1, inplace=True)
    df1.set_index("Mutation", inplace=True)
    df1.index = df1.index.astype(str)
    df1['sort_key'] = df1.index.str.extract(r'(\d+)')[0].astype(int)
    df1 = df1.sort_values(by='sort_key').drop(columns='sort_key')

    plt.figure(figsize=(10, 8))
    
    # Convert to DataFrame and sort by value for better visualization
    diff_series = pd.Series(df1['Prediction ΔΔGStability'])
    diff_series = diff_series.sort_values()

    # Create a multi-row heatmap for better label visibility
    # Determine optimal grid dimensions based on number of mutations
    n_mutations = len(diff_series)
    n_cols = min(20, n_mutations)  # Maximum 20 mutations per row
    n_rows = (n_mutations + n_cols - 1) // n_cols  # Ceiling division
            
    # Create a properly shaped DataFrame for the heatmap
    # Fill with NaN values where needed
    heatmap_data = np.full((n_rows, n_cols), np.nan)
    mutation_labels = np.empty((n_rows, n_cols), dtype=object)
            
    for i, (mutation, value) in enumerate(diff_series.items()):

        row_idx = i // n_cols
        col_idx = i % n_cols
        heatmap_data[row_idx, col_idx] = value
        mutation_labels[row_idx, col_idx] = mutation
            
    # Create a mask for NaN values
    mask = np.isnan(heatmap_data)

    colors = ["red", "white", "green"]
    cmap = LinearSegmentedColormap.from_list("custom_red_white_green", colors)     

    # Create figure with appropriate size
    fig, ax = plt.subplots(figsize=(n_cols * 0.8, n_rows * 0.6))
          
    # Plot heatmap
    sns.heatmap(heatmap_data, annot=False, cmap=cmap, center=0, 
                        mask=mask, cbar_kws={'label': f'ΔΔG differnce from WT)'}, 
                        linewidths=0.5, ax=ax)
            
    # Add mutation labels as x-tick labels and empty y-tick labels
    ax.set_xticks([])
    ax.set_yticks([])
            
    # Add mutation labels directly on the plot
    for i in range(n_rows):
        for j in range(n_cols):
            if not mask[i, j]:  # If not NaN
                ax.text(j + 0.5, i + 0.5, mutation_labels[i, j],
                        ha="center", va="center", color="black",
                        fontsize=8, fontweight='bold')

    plt.title('Variation of Free Gibbs Energy', pad=20)
    plt.ylabel('Mutations')
    plt.xticks(rotation=0)
    plt.yticks(rotation=0)
    plt.tight_layout()
    silentremove(basedir + sep + 'ddg_tetramer_dynamut2.png')
    plt.savefig(basedir + sep + 'ddg_tetramer_dynamut2.pdf', dpi=600)  # Save figure
    plt.show()