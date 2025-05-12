import os 
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt 
from matplotlib.colors import LinearSegmentedColormap


plt.rcParams['font.size'] = 20  # General font size
plt.rcParams['font.weight'] = 'bold'  # All text bold
plt.rcParams['axes.titlesize'] = 26  # Title font size
plt.rcParams['axes.labelsize'] = 18  # X-Y label font size
plt.rcParams['legend.fontsize'] = 18  # Legend font size
plt.rcParams['xtick.labelsize'] = 12  # X-tick font size
plt.rcParams['ytick.labelsize'] = 12  # Y-tick font size

cwd = os.getcwd()
sep = os.sep

results_path = cwd + sep + "results" 
tmscore_monomer_path = results_path + sep + "tm_scores_monomer_all.csv"
tmscore_tetramer_path = results_path + sep + "tm_scores_tetramer_all.csv" 

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def plot_heatmap(df, complex_type):

    # Plotting the heatmap 
    # Create a custom colormap from white to red
    colors = ["red", "white"]
    cmap = LinearSegmentedColormap.from_list("custom_white_red", colors)
    plt.figure(figsize=(10, 8))

    columns = df.columns
    if complex_type == "tetramer":
        columns = [column.split("-")[0] for column in columns]
    df["Mutation"] = columns
    df = df.set_index("Mutation")
    df.columns = columns
    heatmap = sns.heatmap(df, annot=False, fmt='.2f', cmap=cmap, cbar_kws={'label': 'TM-scores'})
    plt.title('TM-scores between TTR mutants in {} conformation'.format(complex_type), pad = 20)
    plt.xlabel('Mutations')
    plt.ylabel('Mutations')
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    plt.tight_layout()
    # Adjust colorbar ticks for spacing
    colorbar = heatmap.collections[0].colorbar
    colorbar.ax.tick_params(pad=10)  # Increase this value for more spacing
    silentremove(results_path + sep + 'tm_scores_all_{}.png'.format(complex_type))
    plt.savefig(results_path + sep + 'tm_scores_all_{}.pdf'.format(complex_type), bbox_inches="tight", dpi=600)  # Save figure
    plt.show()

if __name__ == "__main__":

    tmscore_monomer = pd.read_csv(tmscore_monomer_path)
    tmscore_tetramer = pd.read_csv(tmscore_tetramer_path)
    plot_heatmap(tmscore_monomer, "monomer")
    plot_heatmap(tmscore_tetramer, "tetramer")    
