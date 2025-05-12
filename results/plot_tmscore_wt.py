import os 
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt 


# Set global font properties
plt.rcParams['font.size'] = 20  # General font size
plt.rcParams['font.weight'] = 'bold'  # All text bold
plt.rcParams['axes.titlesize'] = 30  # Title font size
plt.rcParams['axes.labelsize'] = 25  # X-Y label font size
plt.rcParams['legend.fontsize'] = 22  # Legend font size
plt.rcParams['xtick.labelsize'] = 18  # X-tick font size
plt.rcParams['ytick.labelsize'] = 18  # Y-tick font size

cwd = os.getcwd()
sep = os.sep

results_path = cwd + sep + "results" 
tmscore_monomer_path = results_path + sep + "tm_scores_monomer.csv"
tmscore_tetramer_path = results_path + sep + "tm_scores_tetramer.csv" 

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

if __name__ == "__main__":

    tmscore_monomer = pd.read_csv(tmscore_monomer_path)
    tmscore_tetramer = pd.read_csv(tmscore_tetramer_path)
    tmscore_monomer['complex_type'] = 'monomer'
    tmscore_tetramer['complex_type'] = 'tetramer'
    dfs = [tmscore_monomer, tmscore_tetramer]
    df = pd.concat(dfs, ignore_index=True)

    # Plot the distribution
    plot = sns.displot(df, x="TM-score", hue="complex_type", kind="kde", fill=True, legend = False)

    # Customize the plot
    plt.title("TM-score Distribution for TTR monomer and tetramer")
    plt.xlabel("TM-score")
    plt.ylabel("Density")
    
    figure = plot.fig
    silentremove(results_path + sep + "tm_score_wt.png")
    plt.legend(df["complex_type"].unique(), loc='upper left')
    plt.show()
    figure.savefig(results_path + sep + "tm_score_wt.pdf", bbox_inches="tight", dpi = 600)