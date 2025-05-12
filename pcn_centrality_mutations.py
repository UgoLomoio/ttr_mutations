import os 
import ast
import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd 

cwd = os.getcwd()
sep = os.sep 

centrality_measures = ["betweenness", "closeness", "degree_c", "eigenvector_c"]
pcn_outputs = cwd + sep + "ProteinContactNetworks" + sep + "Outputs" + sep + "Centralities" 

mutations_consequence_file = cwd + sep + "mutations_consequences.csv"
df = pd.read_csv(mutations_consequence_file, index_col = 0)

def plot_by(values, basedir, centrality_name, protein_name, mutations=False, zoom_intervals=None):
    global df 

    # Plotting the heatmap
    colors = ["purple", "sienna", "cornflowerblue", "darkolivegreen"]
    
    # If zoom intervals are provided, create multiple subplots
    if zoom_intervals:
        # Calculate how many plots we need (original + zoomed regions)
        num_plots = 1 + len(zoom_intervals)
        fig, axes = plt.subplots(num_plots, 1, figsize=(10, 5*num_plots))
        
        # Main plot (full view)
        ax_main = axes[0]
        for i, (label, value) in enumerate(values.items()):
            ax_main.plot(value, label=label, c=colors[i])
        
        # Add mutations to main plot if provided
        if mutations is not None:
            colors_mut = {
                "Pathogenic": "red",
                "Likely pathogenic": "orange",
                "Variant of uncertain significance": "gold",
                "Likely benign": "lime",
                "Benign": "green",
                "Not found": "black",
                "Unknown": "blue"
            }
            mutations_dict = {}
            y = {}
            c = {}
            for mutation, consequence in df.iterrows():
                consequence = consequence.values[0]
                if consequence not in list(c.keys()): 
                    c[consequence] = []
                    mutations_dict[consequence] = [] 
                    y[consequence] = []

                resi = int(mutation[1:-1])
                mutations_dict[consequence].append(resi)
                y[consequence].append(values["A"][resi])
                c[consequence].append(colors_mut[consequence])
            
            for consequence, mutations_c in mutations_dict.items():
                y_c = y[consequence]
                colors_c = c[consequence]
                ax_main.scatter(x=mutations_c, y=y_c, c=colors_c, label="{} Mutation".format(consequence))
        
        n = len(values["A"])
        ax_main.set_title('{} Centrality Values (Full View)'.format(centrality_name))
        ax_main.set_xlabel('Residue')
        ax_main.set_ylabel('Value')
        ax_main.set_xticks(np.arange(0, n, 10))
        ax_main.set_xticklabels(np.arange(21, 21+n, 10))
        ax_main.legend()
        
        # Create zoomed plots
        for idx, (start, end) in enumerate(zoom_intervals):
            ax_zoom = axes[idx+1]
            
            # Plot data in zoomed view
            for i, (label, value) in enumerate(values.items()):
                ax_zoom.plot(value, label=label, c=colors[i])
            
            # Add mutations to zoomed plot if provided
            if mutations is not None:
                for consequence, mutations_c in mutations_dict.items():
                    y_c = y[consequence]
                    colors_c = c[consequence]
                    ax_zoom.scatter(x=mutations_c, y=y_c, c=colors_c, label="{} Mutation".format(consequence))
            
            # Set the zoom limits
            ax_zoom.set_xlim(start, end)
            
            # Adjust y limits if needed (optional)
            # Calculate appropriate y limits based on the data in the zoomed region
            y_values = []
            for label, value in values.items():
                y_values.extend([v for i, v in enumerate(value) if start <= i <= end])
            
            if y_values:
                y_min = min(y_values) - 0.1 * (max(y_values) - min(y_values))
                y_max = max(y_values) + 0.1 * (max(y_values) - min(y_values))
                ax_zoom.set_ylim(y_min, y_max)
            
            ax_zoom.set_title('{} Centrality Values (Zoom: {}-{})'.format(centrality_name, start+21, end+21))
            ax_zoom.set_xlabel('Residue')
            ax_zoom.set_ylabel('Value')
            
            # Adjust x-ticks for better readability in zoomed region
            tick_spacing = max(1, (end - start) // 10)
            ax_zoom.set_xticks(np.arange(start, end+1, tick_spacing))
            ax_zoom.set_xticklabels(np.arange(start+21, end+21+1, tick_spacing))
            
            # Only show legend if it's not too crowded
            if end - start > 20:
                ax_zoom.legend()
        
        plt.tight_layout()
        plt.savefig('{}_{}_with_zooms.png'.format(centrality_name, protein_name), dpi=300)
        
    else:
        # Original single plot behavior
        plt.figure(figsize=(10, 8))
        for i, (label, value) in enumerate(values.items()):
            plt.plot(value, label=label, c=colors[i])
            
        if mutations is not None:
            colors_mut = {
                "Pathogenic": "red",
                "Likely pathogenic": "orange",
                "Variant of uncertain significance": "gold",
                "Likely benign": "lime",
                "Benign": "green",
                "Not found": "black",
                "Unknown": "blue"
            }
            mutations_dict = {}
            y = {}
            c = {}
            for mutation, consequence in df.iterrows():
                consequence = consequence.values[0]
                if consequence not in list(c.keys()): 
                    c[consequence] = []
                    mutations_dict[consequence] = [] 
                    y[consequence] = []

                resi = int(mutation[1:-1])
                mutations_dict[consequence].append(resi)
                y[consequence].append(values["A"][resi])
                c[consequence].append(colors_mut[consequence])
            
            for consequence, mutations_c in mutations_dict.items():
                y_c = y[consequence]
                colors_c = c[consequence]
                plt.scatter(x=mutations_c, y=y_c, c=colors_c, label="{} Mutation".format(consequence))
        
        n = len(values["A"])
        plt.title('{} Centrality Values'.format(centrality_name))
        plt.xlabel('Residue')
        plt.ylabel('Value')
        plt.xticks(np.arange(0, n, 10), labels=np.arange(21, 21+n, 10))
        plt.tight_layout()
        plt.legend()
        plt.savefig('{}_{}.png'.format(centrality_name, protein_name), dpi=300)


def split_by_chain(centrality_dict):

    new_dict = {}
    for residue, value in centrality_dict.items():
        residue, chain = residue.split(" ")
        resi = residue[3:]
        if chain not in list(new_dict.keys()):
            new_dict[chain] = []
        else: 
            new_dict[chain].append(value)
    return new_dict


def top_k_from_dict(input_dict, k):
    """
    Returns the top k items from a dictionary based on their values.
    
    Parameters:
    input_dict (dict): The input dictionary from which to extract top items.
    k (int): The number of top items to return.
    
    Returns:
    dict: A dict of top K (key, value) sorted by value.
    """
    # Sort the dictionary by value in descending order and get the top k items
    top_k = dict(sorted(input_dict.items(), key=lambda item: item[1], reverse=True)[:k])
    return top_k



if __name__ == "__main__":

    top_k = {}
    for centrality_measure in centrality_measures:
        centrality_output = pcn_outputs + sep + centrality_measure + sep + "Txt"
        files = os.listdir(centrality_output)

        for file in files:
            if "wt" not in file:
                continue
            print(file)
            if file[-3:] == "txt":
                protein = file.split("_")[0]
                if protein not in list(top_k.keys()):
                    top_k[protein] = {}
                content = open(centrality_output + sep + file, 'r').read()
                centrality_dict = ast.literal_eval(content)
                values = split_by_chain(centrality_dict)
                centrality_name = centrality_measure.split("_")[0]
                centrality_name = centrality_name[0].upper() + centrality_name[1:]
                #plot_by(values, centrality_output, centrality_name, protein, mutations = True)
                plot_by(values, centrality_output, centrality_name, protein, mutations = True, zoom_intervals=[(0, 20)])
                if centrality_name not in list(top_k[protein].keys()):
                    temp = list(top_k_from_dict(centrality_dict, 10).keys())
                    top_k[protein][centrality_name] = temp
                
    
    df = pd.DataFrame.from_dict(top_k)
    df.to_csv("top_10_centralities.csv")
    print(df)