import os 
import ast
import matplotlib.pyplot as plt 
import numpy as np 
import seaborn as sns 
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib import gridspec
import re 
import pandas as pd 

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

centrality_measures = ["betweenness", "closeness", "degree_c", "eigenvector_c"]
pcn_outputs = cwd + sep + "ProteinContactNetworks" + sep + "Outputs" + sep + "Centralities" 
complex_types = ["monomer", "tetramer"]

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


def extract_numeric_key(key):
    """
    Extract the numeric part from a key like 'a101v' for sorting.
    Ensure key is a string before processing.
    """
    key = str(key)
    match = re.search(r"\d+", key)
    return int(match.group()) if match else float('inf')

def difference_between_dicts(dict1, dict2):
    """
    Compute the difference between two dictionaries with the same keys.
    Returns a dictionary with key-wise differences.
    """
    if dict1.keys() != dict2.keys():
        raise ValueError("Dictionaries must have the same keys to compute differences")
    
    return {key: dict1[key] - dict2[key] for key in dict1}

def filter_by_chain(data, chain="A"):
    """
    Filters the input dictionary to only include entries for the specified chain.
    
    Parameters:
    - data (dict): Dictionary with keys containing chain identifiers.
    - chain (str): Chain identifier to filter by (default is "A").

    Returns:
    - dict: Filtered dictionary containing only entries for the specified chain.
    """
    return {key: value for key, value in data.items() if key.split()[1] == chain}

def plot_heatmap(dict_of_dicts, title, row_threshold=5, col_threshold=5, max_labels=50, 
                 zoom_intervals=None, zero_threshold=0.01):
    """
    Plot a heatmap filtering both rows and columns based on non-zero counts.
    
    Parameters:
    -----------
    dict_of_dicts : dict
        Dictionary of dictionaries containing the data to plot
    title : str
        Title of the plot
    row_threshold : int
        Minimum number of non-zero values a row (mutant) must have to be included
    col_threshold : int
        Minimum number of non-zero values a column (residue) must have to be included
    max_labels : int
        Maximum number of labels to display on axes
    zoom_intervals : list of tuples or None
        List of (start, end) tuples defining x-axis regions to zoom in on
    zero_threshold : float
        Values with absolute magnitude below this threshold will be considered as zeros
    """
    keys = sorted(next(iter(dict_of_dicts.values())).keys(), key=extract_numeric_key)
    files = sorted(dict_of_dicts.keys(), key=extract_numeric_key)
    
    # Create the full data matrix
    full_data = np.array([[dict_of_dicts[file].get(key, 0) for key in keys] for file in files])
    
    # Count non-zero values for each column (residue) and row (mutant)
    non_zero_cols = np.sum(np.abs(full_data) > zero_threshold, axis=0)
    non_zero_rows = np.sum(np.abs(full_data) > zero_threshold, axis=1)
    
    # Filter columns (residues) that have at least col_threshold non-zero values
    significant_col_indices = np.where(non_zero_cols >= col_threshold)[0]
    
    # Filter rows (mutants) that have at least row_threshold non-zero values
    significant_row_indices = np.where(non_zero_rows >= row_threshold)[0]
    
    # Handle empty results by adjusting thresholds if needed
    if len(significant_col_indices) == 0:
        print(f"Warning: No residues have {col_threshold} or more non-zero values. Adjusting threshold.")
        max_nonzero_cols = np.max(non_zero_cols)
        if max_nonzero_cols > 0:
            significant_col_indices = np.where(non_zero_cols >= max_nonzero_cols)[0]
        else:
            significant_col_indices = np.arange(min(10, len(keys)))
    
    if len(significant_row_indices) == 0:
        print(f"Warning: No mutants have {row_threshold} or more non-zero values. Adjusting threshold.")
        max_nonzero_rows = np.max(non_zero_rows)
        if max_nonzero_rows > 0:
            significant_row_indices = np.where(non_zero_rows >= max_nonzero_rows)[0]
        else:
            significant_row_indices = np.arange(min(10, len(files)))
    
    # Filter the data, keys, and files
    filtered_data = full_data[np.ix_(significant_row_indices, significant_col_indices)]
    filtered_keys = [keys[i] for i in significant_col_indices]
    filtered_files = [files[i] for i in significant_row_indices]
    
    # Create a custom colormap from red to white to green
    colors = ["red", "white", "green"]
    cmap = LinearSegmentedColormap.from_list("custom_red_white_green", colors)
    
    # Find the absolute maximum value to create symmetric bounds
    vmax = max(abs(filtered_data.min()), abs(filtered_data.max()))
    vmin = -vmax
    
    # Use TwoSlopeNorm to properly center the colormap at 0.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    if zoom_intervals is None:
        # Single plot with filtered data
        fig = plt.figure(figsize=(max(8, len(filtered_keys)/10), max(8, len(filtered_files)/5)))
        ax = sns.heatmap(filtered_data, xticklabels=filtered_keys, yticklabels=filtered_files, square=False, 
                         cmap=cmap, annot=False, fmt='.2f', norm=norm,
                         cbar_kws={'label': 'Centrality Difference'})
        plt.xlabel(f"TTR Residues (filtered: ≥{col_threshold} non-zero values)")
        plt.ylabel(f"Mutants (filtered: ≥{row_threshold} non-zero values)")
        
        # Limit the number of labels on axes
        if len(filtered_keys) > max_labels:
            ax.set_xticks(ax.get_xticks()[::max(1, len(filtered_keys)//max_labels)])
            ax.set_xticklabels([filtered_keys[i] for i in range(0, len(filtered_keys), max(1, len(filtered_keys)//max_labels))])
        plt.xticks(rotation=90, fontsize=10)
        plt.yticks(rotation=0, fontsize=8)
        plt.title(f"{title} (Filtered heatmap)")
        colorbar = ax.collections[0].colorbar
        colorbar.ax.tick_params(pad=10)
    else:
        # Create a figure with multiple subplots: main plot + zoomed regions
        n_plots = 1 + len(zoom_intervals)
        fig = plt.figure(figsize=(max(8, len(filtered_keys)/10), 8 * n_plots))
        
        # Create a grid for the subplots
        gs = gridspec.GridSpec(n_plots, 1, height_ratios=[1] + [0.7] * len(zoom_intervals))
        
        # Main plot (filtered view)
        ax_main = plt.subplot(gs[0])
        main_heatmap = sns.heatmap(filtered_data, xticklabels=filtered_keys, yticklabels=files, 
                                  square=False, cmap=cmap, annot=False, fmt='.2f', norm=norm,
                                  cbar_kws={'label': 'Centrality Difference'}, ax=ax_main)
        ax_main.set_xlabel("TTR Residues (filtered)")
        ax_main.set_ylabel("Mutants")
        
        # Limit the number of labels on axes for main plot
        if len(filtered_keys) > max_labels:
            ax_main.set_xticks(ax_main.get_xticks()[::max(1, len(filtered_keys)//max_labels)])
            ax_main.set_xticklabels([filtered_keys[i] for i in range(0, len(filtered_keys), max(1, len(filtered_keys)//max_labels))])
        ax_main.set_xticklabels(ax_main.get_xticklabels(), rotation=90, fontsize=10)
        ax_main.set_yticklabels(ax_main.get_yticklabels(), rotation=0, fontsize=5)
        ax_main.set_title(f"{title} (Showing residues with ≥{row_threshold} non-zero values)")
        
        # Create zoomed plots - adapt zoom intervals to filtered indices
        for i, (start, end) in enumerate(zoom_intervals):
            # Map original residue numbers to filtered indices
            if isinstance(start, str) or isinstance(end, str):
                # If residues are provided as strings (e.g., "A50")
                if start in filtered_keys:
                    start_idx = filtered_keys.index(start)
                else:
                    # Find closest available residue
                    start_idx = 0
                    
                if end in filtered_keys:
                    end_idx = filtered_keys.index(end)
                else:
                    # Find closest available residue
                    end_idx = len(filtered_keys) - 1
            else:
                # If indices are provided directly, map them to filtered indices
                # This is approximate since we've filtered columns
                original_start_key = keys[min(start, len(keys)-1)]
                original_end_key = keys[min(end, len(keys)-1)]
                
                if original_start_key in filtered_keys:
                    start_idx = filtered_keys.index(original_start_key)
                else:
                    start_idx = 0
                    
                if original_end_key in filtered_keys:
                    end_idx = filtered_keys.index(original_end_key)
                else:
                    end_idx = len(filtered_keys) - 1
            
            # Ensure indices are within bounds
            start_idx = max(0, min(start_idx, len(filtered_keys)-1))
            end_idx = max(start_idx, min(end_idx, len(filtered_keys)-1))
            
            # Create zoomed subplot
            ax_zoom = plt.subplot(gs[i+1])
            
            # Extract the data for the zoomed region
            zoom_data = filtered_data[:, start_idx:end_idx+1]
            zoom_keys = filtered_keys[start_idx:end_idx+1]
            
            # Create heatmap for zoomed region
            zoom_heatmap = sns.heatmap(zoom_data, xticklabels=zoom_keys, yticklabels=files,
                                      square=False, cmap=cmap, annot=False, fmt='.2f', norm=norm,
                                      cbar_kws={'label': 'Centrality Difference'}, ax=ax_zoom)
            
            # Adjust labels for zoomed region
            ax_zoom.set_xlabel("TTR Residues (Zoomed)")
            ax_zoom.set_ylabel("Mutants")
            
            # Set appropriate number of ticks for zoomed region
            zoom_max_labels = min(max_labels, len(zoom_keys))
            if zoom_max_labels > 0:
                ax_zoom.set_xticks(np.linspace(0, len(zoom_keys)-1, zoom_max_labels).astype(int))
                ax_zoom.set_xticklabels([zoom_keys[i] for i in np.linspace(0, len(zoom_keys)-1, zoom_max_labels).astype(int)])
            
            ax_zoom.set_xticklabels(ax_zoom.get_xticklabels(), rotation=90, fontsize=10)
            ax_zoom.set_yticklabels(ax_zoom.get_yticklabels(), rotation=0, fontsize=5)
            
            # Set title for zoomed region
            zoom_title = f"{title} (Zoom: {zoom_keys[0]} to {zoom_keys[-1]})"
            ax_zoom.set_title(zoom_title)
        
        plt.tight_layout()
    
    return fig


if __name__ == "__main__":


    means = []
    conf_df = []
    chains_df = []
    centr_df = []
    stds = []
    maxs = []
    mins = []

    for complex_type in complex_types:    
        print(complex_type)
        for centrality_measure in centrality_measures:
            print(centrality_measure)
            centrality_output = pcn_outputs + sep + centrality_measure + sep + "Txt"
            files = os.listdir(centrality_output)
            if complex_type == "monomer":
                files = [file for file in files if "tetramer" not in file]
                chains = ["A"]
            else:
                files = [file for file in files if "tetramer" in file]
                chains = ["A", "B", "C", "D"]


            #print(files)
            centrality_dicts_chains = {}
            for file in files:
                protein = file.split("_")[0]    
                centrality_dicts_chains[protein] = {}
                for chain in chains: 
                    if file[-3:] == "txt":
                        content = open(centrality_output + sep + file, 'r').read()
                        centrality_dict = ast.literal_eval(content)
                        centrality_dict = filter_by_chain(centrality_dict, chain=chain)
                        #print(chain, list(centrality_dict.keys())[0])
                        centrality_dict = {(int(residue.split(" ")[0][3:]) + 20): value for residue, value in centrality_dict.items()}
                        centrality_dicts_chains[protein][chain] = centrality_dict
                 
            if complex_type == "monomer":
                centrality_dict_wt = centrality_dicts_chains["wt"]
            else:
                centrality_dict_wt = centrality_dicts_chains["wt-tetramer"]

            diff_dicts_chain = {}
            for variant, centrality_dict_mutant in centrality_dicts_chains.items():
                for chain, centrality_dict_mutant in centrality_dict_mutant.items():
                    if chain not in diff_dicts_chain:
                        diff_dicts_chain[chain] = {}
                    if "wt" not in variant:
                        if complex_type == "tetramer":
                            variant = variant.split("-")[0]
                        diff_dict = difference_between_dicts(centrality_dict_wt[chain], centrality_dict_mutant)
                        diff_dicts_chain[chain][variant] = diff_dict

            variants = list(centrality_dicts_chains.keys())
            for chain in chains:
                title = "{} centrality".format(centrality_measure)
                fig = plot_heatmap(diff_dicts_chain[chain], title, col_threshold=5, row_threshold=5)
                values_ = [value for value in diff_dicts_chain[chain][variant].values() for variant in variants]
                
                name = "{}_{}_chain{}.pdf".format(centrality_measure, complex_type, chain)
                std = np.std(values_)
                mean = np.mean(values_)
                
                means.append(round(mean, 4))
                stds.append(round(std, 4))
                centr_df.append(centrality_measure.split("_")[0])
                chains_df.append(chain)
                conf_df.append(complex_type)
                maxs.append(round(max(values_), 4))
                mins.append(round(min(values_), 4))
                silentremove(results_path + sep + name.split(".")[0]+".png")
                fig.savefig(results_path + sep + name, bbox_inches="tight", dpi = 600)

    df = pd.DataFrame({
        "Conformation": conf_df,
        "Chain": chains_df,
        "Centrality Measure": centr_df,
        "Mean": means,
        "Std": stds,
        "Max": maxs,
        "Min": mins
    })
    df.to_csv(results_path + sep + "centrality_globalstats.csv")