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

filename = basedir + sep + "Analysis-Mutations3A4D_cleaned_ChainA.xlsx"

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


# Function to shift the numeric index by 20
def shift_mutation(mutation):
    match = re.match(r"([A-Za-z])(\d+)([A-Za-z])", mutation)
    if match:
        prefix, num, suffix = match.groups()
        new_num = str(int(num) + 20)  # Shift by 20
        return f"{prefix}{new_num}{suffix}"
    return mutation  # Return unchanged if pattern doesn't match


# Compute consensus classification
def classify_ddg(value):
    if value < -0.01:
        return "Destabilizing"
    elif value > 0.01:
        return "Stabilizing"
    else:
        return "No Effect"

if __name__ == "__main__":
    # Load the Excel file
    df = pd.read_excel(filename)
    # Apply the function to the Mutation column
    df["Mutation"] = df["Mutation"].apply(shift_mutation)

    # Convert comma decimal format to dot
    for col in df.columns[1:]:  # Exclude 'Mutation' column
        df[col] = df[col].astype(str).str.replace(',', '.').astype(float)
    df.drop("MAESTRO", axis=1, inplace=True)
    
    # Selecting only the numerical columns for correlation analysis
    correlation_columns = df.columns[1:]  # Exclude 'Mutation' column
    correlation_matrix = df[correlation_columns].corr()

    # Save the correlation matrix as a CSV file
    correlation_matrix.to_csv(basedir + sep + 'correlation_matrix.csv')

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt='.2f', linewidths=0.5, annot_kws={'size': 14})
    plt.title('Correlation Heatmap of ΔΔG Predictions', pad=20)
    silentremove(basedir + sep + 'correlation_heatmap.png')
    plt.savefig(basedir + sep + 'correlation_heatmap.pdf', dpi=300)  # Save figure
    plt.show()

    # Apply classification
    classification_df = df.copy()
    for col in correlation_columns:
        classification_df[col] = df[col].apply(classify_ddg)

    # Save the new classification dataframe
    classification_df.to_csv(basedir + sep + 'classification_consensus.csv', index=False)
    
    del classification_df["Mutation"]
    # Create a DataFrame to store agreement counts
    agreement_percentages = pd.DataFrame(index=classification_df.columns, columns=classification_df.columns)

    # Compute agreement
    for col1 in classification_df.columns:
        for col2 in classification_df.columns:
            # Count agreements (same prediction)
            agreement_count = (classification_df[col1] == classification_df[col2]).sum()
            total_count = len(classification_df[col1])  # Total number of mutations
            # Calculate percentage of agreement
            agreement_percentages.loc[col1, col2] = (agreement_count / total_count) * 100
 
    # Convert the agreement percentage to numeric type
    agreement_percentages = agreement_percentages.fillna(0).astype(int)

    # Plotting the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(agreement_percentages, annot=True, fmt='d', cmap='YlGnBu', cbar_kws={'label': 'Agreement Percentage'},  annot_kws={'size': 14})
    plt.title('Agreement Between Prediction Methods', pad=20)
    plt.xlabel('Prediction Methods')
    plt.ylabel('Prediction Methods')
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.tight_layout()
    silentremove(basedir + sep + 'agreement_score_per_model.png')
    plt.savefig(basedir + sep + 'agreement_score_per_model.pdf', dpi=300)  # Save figure
    plt.show()

    # Plotting the heatmap DDG
    # Create a custom colormap from red to white to green
    colors = ["red", "white", "green"]
    cmap = LinearSegmentedColormap.from_list("custom_red_white_green", colors)
    plt.figure(figsize=(10, 8))
    df = df.set_index("Mutation")
    # Step 1: Extract numeric part from the index using regex
    df['sort_key'] = df.index.str.extract('(\d+)').astype(int)
    # Step 2: Sort the DataFrame by the sort_key
    df = df.sort_values(by='sort_key').drop(columns='sort_key')
    heatmap = sns.heatmap(df, annot=False, fmt='.2f', cmap=cmap, cbar_kws={'label': 'ΔΔG'})
    plt.title('Variation of Free Gibbs Energy', pad=20)
    plt.xlabel('Methods')
    plt.ylabel('Mutations')
    plt.rcParams['ytick.labelsize'] = 1 # Y-tick font size
    plt.rcParams['font.weight'] = 'normal'  # All text bold
    plt.rcParams['font.size'] = 1
    
    plt.xticks(rotation=45)
    n = len(list(df.index))
    plt.yticks(np.arange(0, n, 1), labels=list(df.index), rotation=0)
    plt.tight_layout()
    # Adjust colorbar ticks for spacing
    colorbar = heatmap.collections[0].colorbar
    colorbar.ax.tick_params(pad=10)  # Increase this value for more spacing
    silentremove(basedir + sep + 'ddg_models.png')
    plt.savefig(basedir + sep + 'ddg_models.pdf', dpi=300)  # Save figure
    plt.show()





