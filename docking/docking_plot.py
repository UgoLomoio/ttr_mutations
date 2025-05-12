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

sdf_types = ["generated"] #existing, optimized or generated

if "existing" in sdf_types:
    dockers = ["DiffDock", "AutoDock"] #DiffDock or AutoDock
else:
    if "optimized" in sdf_types:
        if "existing" not in sdf_types:
            sdf_types.append("existing")
    if "generated" in sdf_types:
        if "existing" not in sdf_types:
            sdf_types.append("existing")
    dockers = ["AutoDock"]

ligands_dict = {}
for sdf_type in sdf_types:
    print(sdf_type)
    sdf_path = f"sdf-{sdf_type}"
    ligand_dir = docking_path + sep + sdf_path
    if "optimized" in sdf_types:
        if sdf_type == "existing":
            ligands = ["tafamidis"]
            ligands_dict[sdf_type] = ligands
        else:
            ligands = [ligand.split(".")[0] for ligand in os.listdir(ligand_dir) if ".sdf" in ligand]
            ligands_dict[sdf_type] = ligands 
    else:
        ligands = [ligand.split(".")[0] for ligand in os.listdir(ligand_dir) if ".sdf" in ligand]
        ligands_dict[sdf_type] = ligands
    print(f"Ligands {sdf_type}: ", ligands)


output_dir = os.path.join(cwd, "docking", "converted_pdbs")
docked_outputs = os.path.join(cwd, "docking", "AutoDock", "docked-outputs", "existing")  # Adjust path as needed
ligand_reference = os.path.join(output_dir, 'reference_ligand.pdb')

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

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


def plot_normal(values):
    
    # Convert to DataFrame and sort by value for better visualization
    diff_series = pd.Series(values)
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
         
    # Confidence: Red for negative, Green for positive
    cmap = LinearSegmentedColormap.from_list('cmap_confidence', ['red', 'lawngreen', 'white'], N=256)
           
    # Create figure with appropriate size
    fig, ax = plt.subplots(figsize=(n_cols * 1.2, n_rows * 1))
          
    # Plot heatmap
    sns.heatmap(heatmap_data, annot=False, cmap=cmap, center=0, 
                        mask=mask, cbar_kws={'label': 'DiffDock Best Confidence'}, 
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
                        fontsize=14, fontweight='bold')
         
    plt.title(f'DiffDock best docking confidence score between mutants and {ligand}', pad = 0.5)
        
    plt.savefig(results_path + sep + f"conf_{ligand}.pdf", dpi=600, bbox_inches='tight')
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()

def plot_difference(diff, mode = "reference", title = "Title", filename = "diffwt_vina.pdf"):

    modes = ["reference", "optimized"]
    if mode not in modes:
        mode = "reference"

    # Convert to DataFrame and sort by value for better visualization
    diff_series = pd.Series(diff)
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
         
    if docker == "AutoDock":
        # Define colormap: Red for positive, Green for negative
        cmap = LinearSegmentedColormap.from_list('cmap_vina', ['green', 'white', 'red'], N=256)
    else:
        # Confidence: Red for negative, Green for positive
        cmap = LinearSegmentedColormap.from_list('cmap_confidence', ['red', 'white', 'green'], N=256)
           
    # Create figure with appropriate size
    fig, ax = plt.subplots(figsize=(n_cols * 0.8, n_rows * 0.8))
          
    # Plot heatmap
    sns.heatmap(heatmap_data, annot=False, cmap=cmap, center=0, 
                        mask=mask, cbar_kws={'label': title}, 
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
                        fontsize=10, fontweight='bold')
         
    plt.title(title, pad = 0.5)
        
    # Save the figure to a file
    #if docker == "AutoDock":
        #ilentremove(results_path + sep + f"diffwt_vina_{ligand}_{mode}.png")
        #plt.savefig(results_path + sep + f"diffwt_vina_{ligand}_{mode}.pdf", dpi=600, bbox_inches='tight')
    #else:
        #silentremove(results_path + sep + f"diffwt_conf_{ligand}_{mode}.png")
        #plt.savefig(results_path + sep + f"diffwt_conf_{ligand}_{mode}.pdf", dpi=600, bbox_inches='tight')
    # Adjust layout and show the plot
    plt.savefig(results_path + sep + filename, dpi=600, bbox_inches='tight')
    plt.tight_layout()
    plt.show()

def get_best_ligand(ligand_path, ref_ligand_path, output_path=None):
    from pymol import cmd
    import numpy as np
    import os
    
    # Clear PyMOL session
    cmd.delete('all')
    cmd.reinitialize()
    
    # Load reference ligand
    cmd.load(ref_ligand_path, 'ref_ligand')
    cmd.remove('hydrogens')  # Remove hydrogens if present


    print(f"Loaded reference ligand from {os.path.basename(ref_ligand_path)}")
    
    # Load PDBQT file with multiple poses
    cmd.load(ligand_path, 'docked_poses', format='pdbqt')
    n_states = cmd.count_states('docked_poses')
    print(f"Found {n_states} poses in {os.path.basename(ligand_path)}")
    
    # Calculate RMSD for each pose
    rmsd_values = []
    
    for state in range(1, n_states + 1):
        # Create separate object for this state
        obj_name = f'pose_{state}'
        cmd.create(obj_name, 'docked_poses', state, 1)
        
        # Get heavy atom counts for reference and verification
        ref_heavy = cmd.count_atoms('ref_ligand and not elem H')
        pose_heavy = cmd.count_atoms(f'{obj_name} and not elem H')
        print(f"Pose {state-1}: Reference has {ref_heavy} heavy atoms, pose has {pose_heavy} heavy atoms")
        
        try:
            # Use different settings to improve atom matching
            cmd.set('retain_order', 0)  # Don't require same atom order
            
            # Try different alignment approaches
            # 1. Match by atom coordinates rather than names
            cmd.select('ref_heavy', 'ref_ligand and not elem H')
            cmd.select('pose_heavy', f'{obj_name} and not elem H')
            
            # Try alignment with flexible atom matching
            alignment = cmd.align('pose_heavy', 'ref_heavy', cycles=5, transform=1, object=f'aln_{state}')
            
            if alignment and len(alignment) >= 1 and 0 < alignment[0] < 100:
                rmsd = alignment[0]  # RMSD after alignment
                atoms_aligned = alignment[1]  # Number of atom pairs aligned
                print(f"Pose {state-1}: Successfully aligned {atoms_aligned} atoms with RMSD = {rmsd:.4f} Å")
                rmsd_values.append((state-1, rmsd))
            else:
                # Try element-based matching if normal alignment fails
                print(f"Standard alignment failed, trying element-based matching...")
                total_rmsd = 0
                total_atoms = 0
                
                for elem in ['C', 'N', 'O', 'S', 'P', 'F']:
                    cmd.select(f'ref_{elem}', f'ref_ligand and elem {elem}')
                    cmd.select(f'pose_{elem}', f'{obj_name} and elem {elem}')
                    
                    ref_count = cmd.count_atoms(f'ref_{elem}')
                    pose_count = cmd.count_atoms(f'pose_{elem}')
                    
                    if ref_count > 0 and pose_count > 0:
                        try:
                            # Align by element type
                            elem_align = cmd.align(f'pose_{elem}', f'ref_{elem}')
                            if elem_align and elem_align[0] < 100:
                                total_rmsd += elem_align[0] * min(ref_count, pose_count)
                                total_atoms += min(ref_count, pose_count)
                                print(f"  {elem} atoms: RMSD = {elem_align[0]:.4f} Å")
                        except:
                            pass
                
                if total_atoms > 0:
                    avg_rmsd = total_rmsd / total_atoms
                    print(f"Pose {state-1}: Element-based RMSD = {avg_rmsd:.4f} Å")
                    rmsd_values.append((state-1, avg_rmsd))
                else:
                    # If all else fails, at least use molecular centers
                    try:
                        ref_center = cmd.centerofmass('ref_ligand')
                        pose_center = cmd.centerofmass(obj_name)
                        
                        # Calculate distance between centers
                        import math
                        center_dist = math.sqrt(sum([(a-b)**2 for a, b in zip(ref_center, pose_center)]))
                        print(f"Pose {state-1}: Using center distance = {center_dist:.4f} Å")
                        rmsd_values.append((state-1, center_dist))
                    except Exception as e:
                        print(f"Could not calculate any RMSD for pose {state-1}: {e}")
        
        except Exception as e:
            print(f"Error in RMSD calculation for pose {state-1}: {e}")
        
        # Clean up
        cmd.delete(obj_name)
    
    # Find the best pose (lowest RMSD)
    if rmsd_values:
        best_pose = min(rmsd_values, key=lambda x: x[1])
        print(f"Best pose is {best_pose[0]} with RMSD = {best_pose[1]:.4f} Å")
        
        # Save the best pose if requested
        if output_path:
            cmd.create('best_pose', 'docked_poses', best_pose[0]+1, 1)
            cmd.save(output_path, 'best_pose')
            print(f"Best pose saved to {output_path}")
        
        return best_pose[0]
    else:
        print("No valid RMSD values calculated")
        return None

reference_dicts = {ligand: "V50M" for _, ligands in ligands_dict.items() for ligand in ligands}

if __name__ == "__main__":

    dfs = {}

    for sdf_type, ligands in ligands_dict.items():

        dfs[sdf_type] = {}
        for docker in dockers:
            dfs[sdf_type][docker] = {}

            for ligand in ligands:

                ligand = ligand.split(".")[0]
                
                if docker == "DiffDock":
                    docking_csv = docking_path + sep + "docked-outputs" + sep + sdf_type + sep + f"diffdock_{ligand}_existing.csv"
                else:
                    docking_csv = docking_path + sep + "AutoDock" + sep + "docked-outputs" + sep + sdf_type + sep + f"vina_score_docking_{ligand}_TTRtetramer.csv"

                df_melted = create_dataframe(docking_csv, docker)
                df_melted["Mutation"] = [mut.upper() for mut in df_melted["Mutation"]]
                

                if docker == "AutoDock":

                    to_df = []
                    if sdf_type == "existing" or ligand in ["tafamidis", "acoramidis"]:
                        align_ligand = True
                    else:
                        align_ligand = False
                    
                    if not align_ligand:
                        idx = 0
                        score = df_melted[df_melted["Rank Type"] == f"Score_{str(idx+1)}"]
                        df = df_melted.copy()
                        df = df[df["Rank Type"] == f"Score_{str(idx+1)}"]
                        df = df.drop("Rank Type", axis = 1)
                        df.rename(columns={'Rank Value': 'Vina Score'}, inplace=True)
                        df["Mutation"] = [mutant.upper() for mutant in df["Mutation"]]
                        df.set_index('Mutation', inplace=True)
                    else:
                        mutants = df_melted["Mutation"].unique()
                        for mutant in mutants:
                            filename = f"{mutant.lower()}-tetramer_{ligand}.pdbqt"
                            filapath = os.path.join(docked_outputs, filename)
                            idx = get_best_ligand(filapath, ligand_reference, output_path=None)
                            score = df_melted[df_melted["Rank Type"] == f"Score_{str(idx+1)}"]
                            data = [mutant, score["Rank Value"].values[0]]
                            to_df.append(data)
                            #break
                        df = pd.DataFrame(to_df, columns=["Mutation", "Vina Score"])
                    dfs[sdf_type][docker][ligand] = df
                else:
                    df = df_melted.copy()
                    df = df[df["Rank Type"] == "rank_1"]
                    df = df.drop("Rank Type", axis = 1)
                    df.rename(columns={'Rank Value': 'Confidence'}, inplace=True)
                    dfs[sdf_type][docker][ligand] = df

                # Box Plot
                fig_box = px.box(df_melted, x="Mutation", y="Rank Value", color="Mutation", title="Box Plot of Ranks")
                    # Customize the hoverlabel
                fig_box.update_layout(hoverlabel=dict(
                                                font=dict(size=20, family="Arial"),
                                                align="left"))
                
                # Customize the legend
                fig_box.update_layout(legend = dict(font = dict(family = "Courier", size = 30, color = "black")),
                            legend_title = dict(font = dict(family = "Courier", size = 50, color = "blue")))

                fig_box.update_layout(title=dict(text='Box plot {} score for {} and TTR mutants'.format(docker, ligand)),
                                yaxis_zeroline=False, xaxis_zeroline=False, 
                                xaxis_showgrid=False, yaxis_showgrid=False,
                                template='plotly_white',
                                #plot_bgcolor='rgba(0, 0, 0, 0)',
                                #paper_bgcolor='rgba(0, 0, 0, 0)',
                                title_font_size=50,
                                xaxis_title_font_size=40,
                                yaxis_title_font_size=40,
                                xaxis_tickfont_size=30,
                                yaxis_tickfont_size=30,
                                height = 1080,
                                width = 1920,
                )
                file_path = os.path.join(cwd, f"boxplot_{docker}_{ligand}.html")
                fig_box.write_html(file_path)
                # Save the figure as pdf with 600 DPI
                silentremove(cwd + sep + f"boxplot_{docker}_{ligand}.png")
                pio.write_image(fig_box, cwd + sep + f"boxplot_{docker}_{ligand}.pdf", format="pdf", scale=1)  # scale increases resolution

                fig_box.show()          
  

    for docker in dockers:
      for sdf_type, ligands in ligands_dict.items():
        for ligand in ligands:
            ligand = ligand.split(".")[0]
            
            df = dfs[sdf_type][docker][ligand]
            print(docker, sdf_type, ligand)
            print(df)
            if sdf_type == "optimized":
                original_ligand = ligand.split("_")[0]
                if original_ligand in dfs["existing"][docker]:
                    df_no_opt = dfs["existing"][docker][original_ligand]
                    df_no_opt.set_index('Mutation', inplace=True)
                else:
                    raise Exception(f"Original ligand {original_ligand} not in sdf-existing directory.")

            if sdf_type == "generated": 
                reference_ligand = "acoramidis"
                if reference_ligand in dfs["existing"][docker]:
                    df_acora = dfs["existing"][docker][reference_ligand]
                    df_acora.set_index('Mutation', inplace=True)
                else:
                    raise Exception(f"Reference ligand {reference_ligand} not in sdf-existing directory.")
                
            # Step 1: Set index to 'Mutation' and select the reference row
            if "Mutation" in df.columns:
                df.set_index('Mutation', inplace=True)
            reference = reference_dicts[ligand]
            reference_row = df.loc[reference]

            # Drop the reference row after selecting it
            df.drop(reference, inplace=True)

            # Step 2: Compute the difference between reference and all other rows
            if docker == "AutoDock":
                score_label = "Vina Score"
                d = "vina"
            else:
                score_label = "Confidence"
                d = "confidence"
            diff = df[score_label] - reference_row[score_label]
            if docker == "DiffDock":
                values = df[score_label]
                plot_normal(values)
            else:
                plot_difference(diff, mode = "reference", title= f"{docker}: Difference against V50M", filename = f"diffwt_{d}_{ligand}.pdf")

            if sdf_type == "optimized":
                diff_opt = df[score_label] - df_no_opt[score_label]
                #print(diff_opt)
                plot_difference(diff_opt, mode = "optimized", title= f"{docker}: Optimized - Original Tafamidis", filename = f"diffwt_{d}_optimized.pdf")
            elif sdf_type == "generated":
                diff_gen = df[score_label] - df_acora[score_label]
                plot_difference(diff_gen, mode = "optimized", title = f"{docker}: Generated - Acoramidis", filename = f"diffwt_{d}_generated.pdf")
            elif sdf_type == "existing":
                if ligand == "acoramidis":
                    df_tafa = dfs[sdf_type][docker]['tafamidis']
                    df_acora = dfs[sdf_type][docker]['acoramidis']
                    diff_gen = df_acora[score_label] - df_tafa[score_label]
                    plot_difference(diff_gen, mode = "optimized", title= f"{docker}: Acoramidis - Tafamidis", filename = f"diffwt_{d}_acoramidis_tafa.pdf")
            
            values = df.index
            lines = []
            for value in values: 
                lines.append(value + "\n")

            with open(cwd + sep + "sorted_mutations.txt", "w") as file:
                file.writelines(lines)
        

    if sdf_type == "existing" and "AutoDock" in dockers:
        # Find the best Vina score for each ligand
        best_tafamidis_score = dfs[sdf_type]["AutoDock"]['tafamidis']['Vina Score'].mean()
        best_acoramidis_score = dfs[sdf_type]["AutoDock"]['acoramidis']['Vina Score'].mean()

        # Compare to determine which ligand has the better score
        best_ligand = 'tafamidis' if best_tafamidis_score < best_acoramidis_score else 'acoramidis'
        print(f"The best ligand is {best_ligand} with a Vina score. Tafamidis: {best_tafamidis_score} vs Acoramidis: {best_acoramidis_score})")
