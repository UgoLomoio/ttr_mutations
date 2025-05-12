import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.cluster import DBSCAN
import pandas as pd 
import os 
import json 
import re 

plt.rcParams['font.size'] = 20  # General font size
plt.rcParams['font.weight'] = 'bold'  # All text bold
plt.rcParams['axes.titlesize'] = 22  # Title font size
plt.rcParams['axes.labelsize'] = 16  # X-Y label font size
plt.rcParams['legend.fontsize'] = 16  # Legend font size
plt.rcParams['xtick.labelsize'] = 12  # X-tick font size
plt.rcParams['ytick.labelsize'] = 12  # Y-tick font size

cwd = os.getcwd()
sep = os.sep

results_path = cwd + sep + "results" 
tmscore_monomer_path = results_path + sep + "tm_scores_monomer_all.csv"
tmscore_tetramer_path = results_path + sep + "tm_scores_tetramer_all.csv" 

results_path = cwd + sep + "results"
df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["Mutation", "Consequence"]
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
        
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

# Apply the function to the Mutation column
df_consequences["Mutation"] = df_consequences["Mutation"].apply(shift_mutation)

def extract_save_clusters(labels, mapper, complex_type, method):

    # Assuming hc_labels and db_labels are from previous clustering results
    clusters_dict = {}

    # Extract cluster indices from hierarchical clustering
    for idx, cluster_label in enumerate(labels):
        cluster_label = str(cluster_label)
        if cluster_label not in clusters_dict:
            clusters_dict[cluster_label] = []
        clusters_dict[cluster_label].append(mapper[idx])
    
    # Save dictionary to a text file
    with open(results_path + sep + "clusters_{}_{}.txt".format(method, complex_type), "w") as f:
        json.dump(clusters_dict, f, indent=4)

print("Cluster dictionary saved to clusters.txt")

if __name__ == "__main__":

    tmscore_monomer = pd.read_csv(tmscore_monomer_path)
    tmscore_tetramer = pd.read_csv(tmscore_tetramer_path)
    datas = {"monomer": tmscore_monomer, "tetramer": tmscore_tetramer}

    for complex_type, data in datas.items():

        mutations = data.columns
        mutations = [mutation.upper() for mutation in mutations]
        if complex_type == "tetramer":
            mutations = [column.split("-")[0] for column in mutations]
        # Hierarchical Clustering with Average Linkage
        linked = linkage(data, method='average')
        mapper = {i: mutation for i, mutation in enumerate(mutations)}
        # Plot the dendrogram
        plt.figure(figsize=(10, 5))
        dendrogram(linked, labels = mutations)
        plt.title("Hierarchical Clustering Dendrogram for TTR {} (Average Linkage)".format(complex_type))
        plt.xlabel("Mutations")
        plt.ylabel("Distance")


        color_map_consequence = {
            "Pathogenic": "red",
            "Likely pathogenic": "orange",
            "Unknown": "blue",
            "Likely benign": "teal",
            "Benign": "green"
        }
        # Define a mapping dictionary for colors
        color_map_x = {}
        for mutation in mutations:
            if mutation == "WT":
                consequence = "Benign"
            else:
                print(mutation, df_consequences["Mutation"])
                consequence = df_consequences[df_consequences["Mutation"] == mutation]["Consequence"].values[0]
            color_map_x[mutation] = color_map_consequence[consequence]
        # Get the x-tick labels and update their colors
        ax = plt.gca()  # Get the current axis
        for label in ax.get_xticklabels():
            text = label.get_text()
            if text in color_map_x:
                label.set_color(color_map_x[text])

        silentremove(results_path + sep + 'dendogram_{}.png'.format(complex_type))
        plt.savefig(results_path + sep + 'dendogram_{}.pdf'.format(complex_type), bbox_inches="tight", dpi=600)  # Save figure
        #plt.show()

        # Cut the dendrogram at a chosen threshold (e.g., 3 clusters)
        num_clusters = 3
        hc_labels = fcluster(linked, num_clusters, criterion='maxclust')
        extract_save_clusters(hc_labels, mapper, complex_type, "Hierarchical clustering Average Linkage")
        
        # Apply DBSCAN clustering
        dbscan = DBSCAN(eps=0.1, min_samples=50, metric = 'euclidean')
        db_labels = dbscan.fit_predict(data)
        extract_save_clusters(db_labels, mapper, complex_type, "DBSCAN clustering")

        # Visualize the clustering results
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Convert to NumPy array (if not already)
        data = np.array(data)

        # Hierarchical Clustering result
        sns.scatterplot(x=data[:, 0], y=data[:, 1], hue=hc_labels, palette="tab10", ax=axes[0])
        axes[0].set_title("Hierarchical Clustering Dendrogram (Average Linkage)".format(complex_type))

        # DBSCAN result
        sns.scatterplot(x=data[:, 0], y=data[:, 1], hue=db_labels, palette="tab10", ax=axes[1])
        axes[1].set_title("DBSCAN Clustering for TTR {} (Average Linkage)".format(complex_type))
        silentremove(results_path + sep + 'scatterplot_{}.png'.format(complex_type))
        plt.savefig(results_path + sep + 'scatterplot_{}.png'.format(complex_type), bbox_inches="tight", dpi=600)  # Save figure
        plt.show()

        # Output hierarchical and DBSCAN clustering results as a matrix
        results_matrix = np.column_stack((hc_labels, db_labels))
        #print("Hierarchical Clustering Labels and DBSCAN Labels:\n", results_matrix)
