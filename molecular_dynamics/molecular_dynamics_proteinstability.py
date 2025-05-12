import os
import subprocess
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import glob 
from gromacs.fileformats.xpm import XPM

dirname = os.path.dirname(os.path.abspath(__file__))
if dirname != "molecular_dynamics":
    os.chdir(dirname)

cwd = os.getcwd()
sep = os.sep 
print(f"Current working directory: {cwd}")

mol_names = {"tafamidis": "3mi", "acoramidis": "16v"}

parent_dir = os.path.join(cwd, os.pardir)
parent_dir = os.path.abspath(parent_dir)

tetramer = True

# Define paths and parameters
type_sdf = "existing"
OUT_DIR = cwd + sep + "outputs"
if tetramer:
    PROTEIN_DIR = parent_dir + sep + "pdbs-alphafold" + sep + "tetramer"
else:
    PROTEIN_DIR = parent_dir + sep + "pdbs-alphafold" + sep + "monomer"
LIGAND_DIR = parent_dir + sep + "docking" + sep + f"sdf-{type_sdf}"
GROMACS_DIR = "/usr/local/gromacs/bin/gmx"
TRAJ_DIR = cwd + sep + "trajectories"
ENER_DIR = cwd + sep + "energies"
RMSD_DIR = cwd + sep + "rmsd"
MDP_DIR = cwd + sep + "mdp_files"
TEMP_DIR = cwd + sep + "temp_files"
TRAJ_PLOT = cwd + sep + "traj_plots"
LIGAND_TOP = cwd + sep + "ligands_topologies"
ENER_PLOT = cwd + sep + "energies_plots"
MUTANT_TOP = cwd + sep + "topologies"
CLUST_PLOT = cwd + sep + "cluster_plots"
PCA_PLOT = cwd + sep + "pca_plots"

ligands = [file.split(".")[0] for file in os.listdir(LIGAND_DIR) if file.split(".")[-1] == "sdf"]
protein_todo = ["wt", "a45t", "v50m", "v142i", "t80a"]
if protein_todo is not None:
    mutants = []
    for mutant in protein_todo:
        mutants.append(mutant+"-tetramer")
else:
    mutants = [file.split(".")[0] for file in os.listdir(PROTEIN_DIR) if file.split(".")[-1] == "pdb"]


# Ensure log directory exists
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(TRAJ_DIR, exist_ok=True)
os.makedirs(ENER_DIR, exist_ok=True)
os.makedirs(RMSD_DIR, exist_ok=True)
os.makedirs(TRAJ_PLOT, exist_ok=True)
os.makedirs(MDP_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)
os.makedirs(LIGAND_TOP, exist_ok=True)
os.makedirs(ENER_PLOT, exist_ok=True)
os.makedirs(MUTANT_TOP, exist_ok=True)
os.makedirs(CLUST_PLOT, exist_ok=True)
os.makedirs(PCA_PLOT, exist_ok=True)

def plot_rmsd(t, rmsd, out_file):

    # Importing the necessary libraries for plotting
    import matplotlib.pyplot as plt
    import numpy as np

    # Converting time from nanoseconds to picoseconds   
    t = t * 1000

    # Creating a figure and axis objects for the plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    # Filling the area between the x-values (res) and y-values (rmsf) with a semi-transparent black color
    ax.fill_between(t,rmsd, color="black", linestyle="-", alpha=0.3)

    # Plotting the line representing the RMSF values
    ax.plot(t,rmsd, color="black", linestyle="-")

    # Setting labels for the x-axis (time in ps) and y-axis (RMSD value)
    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"C$_\alpha$ RMSD (nm)")

    # Saving the plot as a PDF image with higher resolution (600 dpi)
    plt.savefig(f"{out_file.split('.')[0]}.pdf", format="pdf", dpi=600)

    # Displaying the plot
    #plt.show()

def plot_rmsf(t, rmsf, out_file):

    # Importing the necessary libraries for plotting
    import matplotlib.pyplot as plt
    import numpy as np

    #t is residue number 

    # Creating a figure and axis objects for the plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    # Filling the area between the x-values (res) and y-values (rmsf) with a semi-transparent black color
    ax.fill_between(t,rmsf, color="black", linestyle="-", alpha=0.3)

    # Plotting the line representing the RMSF values
    ax.plot(t,rmsf, color="black", linestyle="-")

    # Setting labels for the x-axis and y-axis (RMSF value)
    ax.set_xlabel('Residue ID (atom)')
    ax.set_ylabel(r"C$_\alpha$ RMSF (nm)")

    # Saving the plot as a PDF image with higher resolution (600 dpi)
    plt.savefig(f"{out_file.split('.')[0]}.pdf", format="pdf", dpi=600)

    # Displaying the plot
    #plt.show()

def plot_energy(t, em, out_file):

    # Importing the necessary libraries for plotting
    import matplotlib.pyplot as plt
    import numpy as np

    # Converting time from nanoseconds to picoseconds   
    t = t * 1000

    # Creating a figure and axis objects for the plot
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)

    # Filling the area between the x-values (res) and y-values (em) with a semi-transparent black color
    #ax.fill_between(t, em, color="black", linestyle="-", alpha=0.3)

    # Plotting the line representing the Potential Energy values
    ax.plot(t, em, color="black", linestyle="-")

    # Setting labels for the x-axis (time in ps) and y-axis (EM value)
    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"C$_\alpha$ Potential Energy (kJ mol−1)")

    # Saving the plot as a PDF image with higher resolution (600 dpi)
    plt.savefig(f"{out_file.split('.')[0]}.pdf", format="pdf", dpi=600)

    # Displaying the plot
     #plt.show()


def plot_cluster(filename, mutant, ligand):

    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    # Read the XPM file
    matrix = read_xpm(filename)
    matrix_pca = pca.fit_transform(matrix)
    
    # Perform KMeans clustering
    n_clusters = 3
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    kmeans.fit(matrix)
    labels = kmeans.labels_
    silhouette_avg = silhouette_score(matrix, labels)
    print(f"Silhouette Score: {silhouette_avg}")
    
    # Plot the clusters
    plt.figure(figsize=(8, 6))
    plt.scatter(matrix_pca[:, 0], matrix_pca[:, 1], c=labels, cmap='viridis', s=50)
    if ligand is None:
        plt.title(f"KMeans Clustering for {mutant}")
    else:
        plt.title(f"KMeans Clustering for {mutant} and {ligand}")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.colorbar(label='Cluster Label')
    if ligand is None:
        plt.savefig(f"{CLUST_PLOT}/{mutant}_clusters.pdf", format="pdf", dpi=600)
    else:
        plt.savefig(f"{CLUST_PLOT}/{mutant}_{ligand}_clusters.pdf", format="pdf", dpi=600)
     #plt.show()

def prepare_mol(mol_file):
    filename = mol_file.split(".")[0]
    mol_file_fix = os.path.join(TEMP_DIR, f"{filename}_fix.mol2")
    cmd = f"perl sort_mol2_bonds.pl {mol_file} {mol_file_fix}"
    subprocess.run(cmd, shell=True, check=True)
    return mol_file_fix

def read_xpm(filename):
    
    if not os.path.exists(filename):
        raise Exception(f"File {filename} does not exist.")
    else:
        print(f"Reading file {filename}")
        matrix = XPM(filename, reverse=True)
        return matrix.array
    

def plot_pca(PSF, DCD, mutant, ligand, n_components=3):

    import MDAnalysis as mda
    from MDAnalysis.analysis import pca, align
    import warnings
    import pandas as pd 

    # suppress some MDAnalysis warnings about writing PDB files
    warnings.filterwarnings('ignore')

    u = mda.Universe(PSF, DCD)
    #print("Aligning trajectory")
    aligner = align.AlignTraj(u, u, select='backbone',
                            in_memory=True).run()
    #print("Calculating PCA")
    pc = pca.PCA(u, select='backbone',
             align=True, mean=None,
             n_components=None).run()
    #print("PCA done")
    backbone = u.select_atoms('backbone')
    n_bb = len(backbone)
    #print('There are {} backbone atoms in the analysis'.format(n_bb))
    #print(pc.p_components.shape)

    fig = plt.figure(figsize=(10, 5))
    plt.plot(pc.cumulated_variance[:10])
    plt.xlabel('Principal component')
    plt.ylabel('Cumulative variance')
    plt.savefig(f"{PCA_PLOT}/{mutant}_{ligand}_pca_cumvar.pdf", dpi=600)

    transformed = pc.transform(backbone, n_components=n_components)
    df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(n_components)])
    df['Time (ps)'] = df.index * u.trajectory.dt
    
    g = sns.PairGrid(df, hue='Time (ps)',
                 palette=sns.color_palette('Oranges_d',
                                           n_colors=len(df)))
    g.map(plt.scatter, marker='.')
    plt.savefig(f"{PCA_PLOT}/{mutant}_{ligand}_pca.pdf", dpi=600)

    pc1 = pc.p_components[:, 0]
    trans1 = transformed[:, 0]
    projected = np.outer(trans1, pc1) + pc.mean.flatten()
    coordinates = projected.reshape(len(trans1), -1, 3)
    proj1 = mda.Merge(backbone)
    proj1.load_new(coordinates, order="fac")
    return proj1
    
def run_gromacs_experiment(mutant, ligand):
    """Run a GROMACS experiment for a given mutant and ligand"""


    if ligand is not None:
            mutant_solv = f"{TEMP_DIR}/{mutant}_{ligand}_solv.gro"
            mutant_ions = f"{TEMP_DIR}/{mutant}_{ligand}_ions.gro"
            em_gro = f"{TEMP_DIR}/{mutant}_{ligand}_em.gro"
            em_trr = f"{TEMP_DIR}/{mutant}_{ligand}_em.trr"
            em_edr = f"{ENER_DIR}/{mutant}_{ligand}_em.edr"
            em_log = f"{TEMP_DIR}/{mutant}_{ligand}_em.log"
            nvt_gro = f"{TEMP_DIR}/{mutant}_{ligand}_nvt.gro"
            nvt_trr = f"{TEMP_DIR}/{mutant}_{ligand}_nvt.trr"
            nvt_edr = f"{TEMP_DIR}/{mutant}_{ligand}_nvt.edr"
            nvt_log = f"{TEMP_DIR}/{mutant}_{ligand}_nvt.log"
            npt_gro = f"{TEMP_DIR}/{mutant}_{ligand}_npt.gro"
            npt_trr = f"{TEMP_DIR}/{mutant}_{ligand}_npt.trr"
            npt_edr = f"{TEMP_DIR}/{mutant}_{ligand}_npt.edr"
            npt_log = f"{TEMP_DIR}/{mutant}_{ligand}_npt.log"
            md_trr = f"{TRAJ_DIR}/{mutant}_{ligand}.trr"
            md_xtc = f"{TRAJ_DIR}/{mutant}_{ligand}.xtc"
            md_edr = f"{ENER_DIR}/{mutant}_{ligand}_md.edr"
            md_log = f"{MDP_DIR}/{mutant}_{ligand}_md.log"
            md_xpm = f"{TRAJ_DIR}/{mutant}_{ligand}.xpm"
            md_tpr = f"{TRAJ_DIR}/{mutant}_{ligand}.tpr"
            md_gro = f"{MUTANT_TOP}/{mutant}_{ligand}_md.gro"
            out_file = f"{RMSD_DIR}/{mutant}_{ligand}.xvg"

    else: 
            mutant_solv = f"{TEMP_DIR}/{mutant}_solv.gro"
            mutant_ions = f"{TEMP_DIR}/{mutant}_ions.gro"
            em_gro = f"{TEMP_DIR}/{mutant}_em.gro"
            em_trr = f"{TEMP_DIR}/{mutant}_em.trr"
            em_edr = f"{ENER_DIR}/{mutant}_em.edr"
            em_log = f"{TEMP_DIR}/{mutant}_em.log"
            nvt_gro = f"{TEMP_DIR}/{mutant}_nvt.gro"
            nvt_trr = f"{TEMP_DIR}/{mutant}_nvt.trr"
            nvt_edr = f"{TEMP_DIR}/{mutant}_nvt.edr"
            nvt_log = f"{TEMP_DIR}/{mutant}_nvt.log"
            npt_gro = f"{TEMP_DIR}/{mutant}_npt.gro"
            npt_trr = f"{TEMP_DIR}/{mutant}_npt.trr"
            npt_edr = f"{TEMP_DIR}/{mutant}_npt.edr"
            npt_log = f"{TEMP_DIR}/{mutant}_npt.log"
            md_trr = f"{TRAJ_DIR}/{mutant}.trr"
            md_xtc = f"{TRAJ_DIR}/{mutant}.xtc"
            md_edr = f"{ENER_DIR}/{mutant}_md.edr"
            md_log = f"{MDP_DIR}/{mutant}_md.log"
            md_xpm = f"{TRAJ_DIR}/{mutant}.xpm"
            md_tpr = f"{TRAJ_DIR}/{mutant}.tpr"
            md_gro = f"{MUTANT_TOP}/{mutant}_md.gro"
            out_file = f"{RMSD_DIR}/{mutant}.xvg"

    # Check if log file already exists
    count = 0
    if os.path.exists(out_file):
        print(f"Log file for {mutant} and {ligand} RMSDs already exists. Skipping experiment.")
        t, rmsd = read_xvg(out_file)
        if ligand is None:
            plot_rmsd(t, rmsd, f"{TRAJ_PLOT}/{mutant}_rmsd.pdf")
        else:
            plot_rmsd(t, rmsd, f"{TRAJ_PLOT}/{mutant}_{ligand}_rmsd.pdf")
        count += 1
    if os.path.exists(md_xpm):
        print(f"Log file for {mutant} and {ligand} RMSD matrix already exists. Skipping experiment.")
        #matrix = read_xpm(md_xpm)
        #plot_heatmap(matrix, ligand, mode = "steps", mutant=mutant)
        count += 1
    if os.path.exists(f"{em_edr.split('.')[0]}_potential_energy.xvg"):
        print(f"Log file for {mutant} and {ligand} energy already exists. Skipping experiment.")
        filename = f"{em_edr.split('.')[0]}_potential_energy.xvg"
        t_e, em = read_xvg(filename)
        if ligand is None:
            plot_energy(t_e, em, f"{ENER_PLOT}/{mutant}_energy.pdf")
        else:
            plot_energy(t_e, em, f"{ENER_PLOT}/{mutant}_{ligand}_energy.pdf")
        count += 1
    if count == 3:
        print(f"All log files for {mutant} and {ligand} already exist. Skipping experiment.")
        t, rmsd = read_xvg(out_file)
        #skipping this part
        """
        if ligand is None:
            filename = f"{TEMP_DIR}/{mutant}_cluster.xpm"
        else:
            filename = f"{TEMP_DIR}/{mutant}_{ligand}_cluster.xpm"
        if os.path.exists(filename):
            print(f"Cluster file for {mutant} and {ligand} already exists. Skipping experiment.")
        else:
            cmd_cluster = f"gmx cluster -s {npt_gro} -dm {md_xpm} -o {filename} -tu ns"
            subprocess.run(cmd_cluster, shell=True, check=True, stdout=subprocess.DEVNULL)
        plot_cluster(filename, mutant, ligand)

        # Visualize frames
        if ligand is None:
            filename = f"{TRAJ_PLOT}/{mutant}_frames.pdf"
        else:
            filename = f"{TRAJ_PLOT}/{mutant}_{ligand}_frames.pdf"
        print(f"Creating PCA projections for {mutant} and {ligand}")
        plot_pca(md_gro, md_xtc, mutant, ligand)
        return t, rmsd
        """
        return t, rmsd
    # Run GROMACS experiment
    try:

        if ligand is not None:
            npt_mdp = f"npt_c.mdp"  
            nvt_mdp = f"nvt_c.mdp"
            md_mdp = f"md_c.mdp" 
            em_mdp = f"em_c.mdp"
            ions_mdp = f"ions_c.mdp"

        else:
            npt_mdp = f"npt.mdp"  
            nvt_mdp = f"nvt.mdp"
            md_mdp = f"md.mdp" 
            em_mdp = f"em.mdp"
            ions_mdp = f"ions.mdp"

        print(f"Running GROMACS experiment on transthyretin: {mutant}")
        mutant_path = f"{PROTEIN_DIR}/{mutant}.pdb"
        mutant_processed_path = f"{TEMP_DIR}/{mutant}.gro"
        topol_file = f"{mutant}.top"  

        # Use TIP3P water model and specify force field (Amber99SB-ILDN or CHARMM36)
        cmd_gmx = f"gmx pdb2gmx -f {mutant_path} -o {mutant_processed_path} -water tip3p -p {topol_file} -ff charmm36-jul2022"
        subprocess.run(cmd_gmx, shell=True, check=True, stdout=subprocess.DEVNULL)

        if ligand is not None:
            
            ligand_name = ligand.split(".")[0]
            mol_name = mol_names[ligand_name]
            # Fix the ligand PDB and GRO file names
            lig_pdb = f"{LIGAND_TOP}{sep}{mol_name}_ini.pdb" 
            lig_gro = f"{mol_name}.gro"
            lig_top = f"{mol_name}.top"
            lig_ipt = f"{LIGAND_TOP}{sep}{mol_name}.itp"
            lig_prm = f"{LIGAND_TOP}{sep}{mol_name}.prm"

            # Convert ligand PDB to GRO format
            if f"{ligand_name}_gmx.gro" not in os.listdir(LIGAND_TOP):
                cmd_mol_edit = f"gmx editconf -f {lig_pdb} -o {lig_gro}"
                subprocess.run(cmd_mol_edit, shell=True, check=True, stdout=subprocess.DEVNULL)
            else:
                print("Skipping GRO creation, already exists")

            complex_gro = f"{TEMP_DIR}/{mutant}_{ligand}.gro"
            complex_top_file = f"{mutant}_{ligand}.top"

            # Copy mutant GRO coordinates to complex
            with open(mutant_processed_path, "r") as mutant_gro:
                lines_mut = mutant_gro.readlines()

            # Copy ligand GRO coordinates to complex
            coordinates_lig = []
            with open(lig_gro, "r") as ligand_gro:
                lines = ligand_gro.readlines()
                for line in lines:
                    if mol_name.upper() in line:
                        coordinates_lig.append(line)
                n_atoms_lig = len(coordinates_lig)
                
            # Add ligand to complex GRO file
            with open(complex_gro, "a") as complex_gro_file:
                
                # Write mutant coordinates
                complex_gro_file.write(lines_mut[0])
                n_atoms = int(lines_mut[1]) + n_atoms_lig
                complex_gro_file.write(" "+str(n_atoms)+"\n")
                complex_gro_file.writelines(lines_mut[2:-1])
                #complex_gro_file.write("\n")
                # Write ligand coordinates
                for line in coordinates_lig:
                    complex_gro_file.write(line)
                # write final line
                #complex_gro_file.write("\n")
                complex_gro_file.write(lines_mut[-1])
            
            # Copy mutant topology to complex
            with open(topol_file, "r") as mutant_top:
                lines = mutant_top.readlines()

            #Add ligand to topology and parameters to complex topology
            lines_to_add1 = ['\n', '; Include ligand topology\n', f'#include "{lig_ipt}"\n']
            lines_to_add2 = ['\n', '; Include ligand parameters\n', f'#include "{lig_prm}"\n']

            with open(complex_top_file, "w") as complex_top:    
                for line in lines: 
                    complex_top.write(line)
                    if "./charmm36-jul2022.ff/forcefield.itp" in line:
                        for line_to_add in lines_to_add2:
                            complex_top.write(line_to_add)
                    if "#endif" in line:
                        for line_to_add in lines_to_add1:
                            complex_top.write(line_to_add)
                complex_top.write(f"{mol_name.upper()} \t\t 1 \n")
                
            mutant_processed_path = complex_gro
            topol_file = complex_top_file
            print(f"Successfully prepared ligand {ligand_name} and updated topology")


        # Create a box with nm distance
        protein_gro = f"{TEMP_DIR}/{mutant}_newbox.gro"
        cmd_editconf = f"gmx editconf -f {mutant_processed_path} -o {protein_gro} -c -d 2 -bt cubic"#octahedron
        subprocess.run(cmd_editconf, shell=True, check=True, stdout=subprocess.DEVNULL)
        
        if md_xtc in os.listdir(TRAJ_DIR):
            print(f"Log file for {mutant} and {ligand} already exists. Skipping experiment.")
        else:
            # Solvate with TIP3P water model
            cmd_solvate = f"gmx solvate -cp {protein_gro} -cs tip3p.gro -o {mutant_solv} -p {topol_file}"
            subprocess.run(cmd_solvate, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Add ions at 0.15M NaCl concentration
            ions_tpr = f"{TEMP_DIR}/ions.tpr"
            cmd_grompp = f"gmx grompp -f {ions_mdp} -c {mutant_solv} -p {topol_file} -o {ions_tpr}"
            subprocess.run(cmd_grompp, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Add ions with genion (neutral system + 0.15M NaCl)
            cmd_genion = f"echo 'SOL' | gmx genion -s {ions_tpr} -o {mutant_ions} -p {topol_file} -pname NA -nname CL -neutral -conc 0.15"
            subprocess.run(cmd_genion, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Energy minimization
            em_tpr = f"{TEMP_DIR}/em.tpr"
            cmd_grompp_em = f"gmx grompp -f {em_mdp} -c {mutant_ions} -p {topol_file} -o {em_tpr} -maxwarn 1"
            subprocess.run(cmd_grompp_em, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Run energy minimization
            cmd_mdrun_em = f"gmx mdrun -s {em_tpr} -c {em_gro} -o {em_trr} -e {em_edr} -g {em_log} -v"
            subprocess.run(cmd_mdrun_em, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Calculate potential energy
            xvg_file = f"{em_edr.split('.')[0]}_potential_energy.xvg"
            cmd_energy = f"echo 'Potential' | gmx energy -f {em_edr} -o {xvg_file}"
            subprocess.run(cmd_energy, shell=True)
            t_e, em = read_xvg(xvg_file)
            if ligand is None:
                filename = f"{ENER_PLOT}/{mutant}_energy.pdf"
            else:
                filename = f"{ENER_PLOT}/{mutant}_{ligand}_energy.pdf"
            plot_energy(t_e, em, filename)

            if ligand is not None:
                # Create index file for ligand
                lig_index = f"{TEMP_DIR}/index_{mol_name}.ndx"
                index_lig_cmd = f"echo 'q' | gmx make_ndx -f {lig_gro} -o {lig_index}"
                subprocess.run(index_lig_cmd, shell=True, check=True)  
                
                # Generate position restraints
                genrestrs_cmd = f"echo '0' | gmx genrestr -f {lig_gro} -n {lig_index} -o posre_{mol_name}.itp -fc 1000 1000 1000"
                subprocess.run(genrestrs_cmd, shell=True, check=True)
                
                # Include the ligand position restraint in the topology file
                lines_to_add_posre = ['\n', '; Ligand position restraint\n', '#ifdef POSRES\n', f'#include "posre_{mol_name}.itp"\n', '#endif\n']
                
                # Correct way to modify the topology file
                with open(topol_file, "r") as top_read:
                    lines = top_read.readlines()
                
                with open(topol_file, "w") as top_write:
                    for line in lines:
                        top_write.write(line)
                        if line.strip() == f'#include "{mol_name}.itp"':
                            for line_add in lines_to_add_posre:
                                top_write.write(line_add)
                
                # Create a temporary file with make_ndx commands
                with open(f"{TEMP_DIR}/make_ndx_input.txt", "w") as f:
                    f.write("1 | 13\n")  # Combine protein and ligand
                    f.write("name 21 Protein_LIG\n")  # Name the new group
                    f.write("11 & ! 13\n")  # Create Water_and_ions (non-Protein without ligand)
                    f.write("name 22 Water_and_ions_only\n")
                    f.write("q\n")  # Quit

                # Run make_ndx with the input file
                index_file = f"{TEMP_DIR}/index.ndx"
                cmd_merge = f"gmx make_ndx -f {em_gro} -o {index_file} < {TEMP_DIR}/make_ndx_input.txt"
                subprocess.run(cmd_merge, shell=True, check=True)


            # NVT equilibration
            nvt_tpr = f"{TEMP_DIR}/nvt.tpr"
            if ligand is not None:
                cmd_grompp_nvt = f"gmx grompp -f {nvt_mdp} -c {em_gro} -p {topol_file} -n {index_file} -r {em_gro} -o {nvt_tpr} -maxwarn 1" 
            else:
                cmd_grompp_nvt = f"gmx grompp -f {nvt_mdp} -c {em_gro} -p {topol_file} -r {em_gro} -o {nvt_tpr} -maxwarn 1" 
            subprocess.run(cmd_grompp_nvt, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Run NVT equilibration
            cmd_mdrun_nvt = f"gmx mdrun -s {nvt_tpr} -c {nvt_gro} -o {nvt_trr} -e {nvt_edr} -g {nvt_log} -v"
            subprocess.run(cmd_mdrun_nvt, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Sometimes, NVT equilibration can cause dissociation of the tetramer into two dimers. One of them moving out of the box.
            # To fix this issue we do as follows:
            # 1. Make the whole protein again
            cmd_fixnvt = f"echo '0\n0' | gmx trjconv -s {nvt_tpr} -f {nvt_gro} -o {nvt_gro} -pbc whole"
            subprocess.run(cmd_fixnvt, shell=True, check=True, stdout=subprocess.DEVNULL)
            # 2. Center the protein in the box
            cmd_centernvt = f"echo '0\n0' | gmx trjconv -s {nvt_tpr} -f {nvt_gro} -o {nvt_gro} -pbc mol -center"
            subprocess.run(cmd_centernvt, shell=True, check=True, stdout=subprocess.DEVNULL)


            # NPT equilibration
            npt_tpr = f"{TEMP_DIR}/npt.tpr"
            if ligand is not None:
                cmd_grompp_npt = f"gmx grompp -f {npt_mdp} -c {nvt_gro} -r {nvt_gro} -p {topol_file} -n {index_file} -o {npt_tpr} -maxwarn 1"
            else:
                cmd_grompp_npt = f"gmx grompp -f {npt_mdp} -c {nvt_gro} -r {nvt_gro} -p {topol_file} -o {npt_tpr} -maxwarn 1"
            subprocess.run(cmd_grompp_npt, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Run NPT equilibration
            cmd_mdrun_npt = f"gmx mdrun -s {npt_tpr} -c {npt_gro} -o {npt_trr} -e {npt_edr} -g {npt_log} -v"
            subprocess.run(cmd_mdrun_npt, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Sometimes, NPT equilibration can cause dissociation of the tetramer into two dimers. One of them moving out of the box.
            # To fix this issue we do as follows:
            # 1. Make the whole protein again
            cmd_fixnpt = f"echo '0\n0' | gmx trjconv -s {npt_tpr} -f {npt_gro} -o {npt_gro} -pbc whole"
            subprocess.run(cmd_fixnpt, shell=True, check=True, stdout=subprocess.DEVNULL)
            # 2. Center the protein in the box
            cmd_centernpt = f"echo '0\n0' | gmx trjconv -s {npt_tpr} -f {npt_gro} -o {npt_gro} -pbc mol -center"
            subprocess.run(cmd_centernpt, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Production MD
            if ligand is not None:
                cmd_grompp_md = f"gmx grompp -f {md_mdp} -c {npt_gro} -p {topol_file} -n {index_file} -o {md_tpr} -maxwarn 1"
            else:
                cmd_grompp_md = f"gmx grompp -f {md_mdp} -c {npt_gro} -p {topol_file} -o {md_tpr} -maxwarn 1"
            subprocess.run(cmd_grompp_md, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Run production MD with GPU acceleration
            cmd_mdrun_md = f"gmx mdrun -s {md_tpr} -c {md_gro} -o {md_trr} -x {md_xtc} -e {md_edr} -g {md_log} -v -ntmpi 1 -nb gpu"
            subprocess.run(cmd_mdrun_md, shell=True, check=True, stdout=subprocess.DEVNULL)

            #Sometimes, PBC artifacts can occur in the trajectory. To fix this issue we do as follows:
            # 1. Make the whole protein again 
            cmd_trjconv1 = f"echo '0\n0' | gmx trjconv -s {md_tpr} -f {md_xtc} -o {md_xtc} -pbc whole -tu ns"
            subprocess.run(cmd_trjconv1, shell=True, check=True, stdout=subprocess.DEVNULL)
            # 2. Remove jumping artifacts accross the periodic boundaries
            cmd_trjconv2 = f"echo '0\n0' | gmx trjconv -s {md_tpr} -f {md_xtc} -o {md_xtc} -pbc nojump -tu ns"
            subprocess.run(cmd_trjconv2, shell=True, check=True, stdout=subprocess.DEVNULL)
            # 3. Center the protein in the box
            cmd_trjconv3 = f"echo '0\n0' | gmx trjconv -s {md_tpr} -f {md_xtc} -o {md_xtc} -center -pbc mol -ur compact"
            subprocess.run(cmd_trjconv3, shell=True, check=True, stdout=subprocess.DEVNULL)

            # Quite same treatment for the gro topology file:
            # 1. Make the whole protein topology 
            cmd_fixgro = f"echo '0\n0' | gmx trjconv -s {md_tpr} -f {md_gro} -o {md_gro} -pbc whole"
            subprocess.run(cmd_fixgro, shell=True, check=True, stdout=subprocess.DEVNULL)
            # 2. Center the protein in the box
            cmd_centergro = f"echo '0\n0' | gmx trjconv -s {md_tpr} -f {md_gro} -o {md_gro} -pbc mol -center"
            subprocess.run(cmd_centergro, shell=True, check=True, stdout=subprocess.DEVNULL)


        # Calculate RMSD
        cmd_rms = f"echo '1\n1' | gmx rms -s {md_tpr} -o {out_file} -m {md_xpm} -f {md_xtc} -nopbc -tu ns"
        subprocess.run(cmd_rms, shell=True)

        # Read RMSD data
        t, rmsd = read_xvg(out_file)
        matrix = read_xpm(md_xpm)
        if ligand is None:
            filename = f"{TRAJ_PLOT}/{mutant}_rmsd.pdf"
        else:
            filename = f"{TRAJ_PLOT}/{mutant}_{ligand}_rmsd.pdf"

        # Plot RMSD
        plot_rmsd(t, rmsd, filename)
        # Plot RMSD matrix
        plot_heatmap(matrix, ligand, mode = "steps", mutant=mutant)

        if ligand is None:
            filename = f"{RMSD_DIR}/{mutant}_rmsf.xvg"
            plotname = f"{TRAJ_PLOT}/{mutant}_rmsf.pdf"
        else:
            filename = f"{RMSD_DIR}/{mutant}_{ligand}_rmsf.xvg"
            plotname = f"{TRAJ_PLOT}/{mutant}_{ligand}_rmsf.pdf"

        # Compute RMSF 
        cmd_rmsf = f"echo '1\n1' | gmx rmsf -s {md_tpr} -f {md_xtc} -o {filename}"
        subprocess.run(cmd_rmsf, shell=True, check=True, stdout=subprocess.DEVNULL)
        # Read RMSF data
        t_f, rmsf = read_xvg(filename=filename)
        # Plot RMSF
        plot_rmsf(t_f, rmsf, plotname)

        if ligand is None:
            filename = f"{TEMP_DIR}/{mutant}_cluster.xpm"
        else:
            filename = f"{TEMP_DIR}/{mutant}_{ligand}_cluster.xpm"

        cmd_cluster = f"gmx cluster -s {md_gro} -dm {md_xpm} -o {filename} -tu ns"
        subprocess.run(cmd_cluster, shell=True, check=True, stdout=subprocess.DEVNULL)
        plot_cluster(filename, mutant, ligand)

        """
        # Visualize frames
        if ligand is None:
            filename = f"{TRAJ_PLOT}/{mutant}_frames.pdf"
        else:
            filename = f"{TRAJ_PLOT}/{mutant}_{ligand}_frames.pdf"
        
        visualize_gromacs_frames_old(md_gro, md_xtc, n_frames=10, selection='polymer', output_file=filename)
        """
        plot_pca(md_gro, md_xtc, mutant, ligand, n_components=3)
        #create_movie(md_gro, md_xtc, output_file=filename.replace(".pdf", ".gif"))
        return t, rmsd
    
    except subprocess.CalledProcessError as e:
        print(f"Error running GROMACS for {mutant} and {ligand}: {e}")
        return None

def read_xvg(filename):
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return None
    else:
        print(f"Reading file {filename}")
        x,y = np.loadtxt(filename,comments=["#", "@"],unpack=True)
        return x, y

def compare_rmsd_mutants(rmsds, ligand = None):

    fig = plt.figure(figsize=(10, 8))

    for name, results in rmsds.items():
        t = results["t"]
        rmsd = results["rmsd"]
        plt.plot(t, rmsd, label=name)
        
    plt.xlabel("Time (ps)")
    plt.ylabel("RMSD (nm)")
    plt.title("RMSD of mutants")
    plt.legend()
    if ligand is None:
        plt.savefig(f"{OUT_DIR}/rmsd_mutants.pdf", format="pdf", dpi=600)
    else:
        plt.savefig(f"{OUT_DIR}/{ligand}_rmsd_mutants.pdf", format="pdf", dpi=600)


def compute_rmsd_matrix(ligand):    
    """Compute RMSD matrix for a given ligand across all mutants"""
    
    import pickle 

    matrix = np.zeros((len(mutants), len(mutants)))
    mins = []
    maxs = [] 
    matrixs = {}
    rmsds = {}
    for i, mutant in enumerate(mutants):
        if mutant not in rmsds:

            if ligand is None:
                md_xpm = f"{TRAJ_DIR}/{mutant}.xpm"
            else:
                md_xpm = f"{TRAJ_DIR}/{mutant}_{ligand}.xpm"

            name = f"{mutant}_{ligand}" if ligand is not None else mutant
            rmsds[name] = {"t": None, "rmsd": None}

            t, rmsd = run_gromacs_experiment(mutant, ligand)     
            matrix = read_xpm(md_xpm)
            rmsds[name]["t"] = t
            rmsds[name]["rmsd"] = rmsd

            min_ = np.min(matrix)
            max_ = np.max(matrix)
            mins.append(min_)
            maxs.append(max_)
            matrixs[name] = matrix
        
    min_ = np.min(mins)
    max_ = np.max(maxs)
    for mutant, matrix in matrixs.items():
        plot_heatmap(matrix, ligand, mode = "steps", mutant=mutant, min = min_, max = max_)
    
    compare_rmsd_mutants(rmsds, ligand)
    
    print("Saving RMSD matrix")
    # Calculate RMSD matrix
    saved_mutants = rmsds.keys()
    matrix = np.zeros((len(saved_mutants), len(saved_mutants)))
    for i, mutant1 in enumerate(saved_mutants):
        for j, mutant2 in enumerate(saved_mutants):
            val1 = np.max(rmsds[mutant1]["rmsd"])
            val2 = np.max(rmsds[mutant2]["rmsd"])
            #print(val1, val2)
            matrix[i, j] = matrix[j, i] = abs(val1 - val2)

    with open(f'{RMSD_DIR}/rmsds_dictionary.pkl', 'wb') as f:
        pickle.dump(rmsds, f)
    return matrix

def plot_heatmap(matrix, ligand, mode = "mutants", mutant = None, mutants = None, min = None, max = None):
    """Plot heatmap for RMSD matrix"""

    if mode in ["mutants", "steps"]:
        
        plt.figure(figsize=(10, 8))
        if mode == "mutants":
            mutants = [mutant.split("-")[0] for mutant in mutants]
            ticklabels = [mutant.split("-")[0] for mutant in mutants]
            xtickslabels = ticklabels
            yticklabels = ticklabels
            xlabel = "Mutants"
            ylabel = "Mutants"
            if ligand is not None:
                title = f"Pairwise RMSD Distances for {ligand}"
                filename = f"{OUT_DIR}/{ligand}_heatmap.pdf"
            else:
                title = "Pairwise RMSD Distances without ligand"
                filename = f"{OUT_DIR}/heatmap_noligand.pdf"
        else:
            steps = [str(i) for i in range(1, len(matrix) + 1)]
            title = "Pairwise RMSD Distances for mutant {} trajectories".format(mutant)
            xlabel = "Times (ps)"
            ylabel = "Times (ps)"
            if ligand is not None:
                filename = f"{OUT_DIR}/{mutant}_{ligand}_heatmap_steps.pdf"
            else:
                filename = f"{OUT_DIR}/{mutant}_heatmap_steps.pdf"
            xtickslabels = steps
            yticklabels = steps
        
        matrix = matrix.astype(np.float32)
        annot = False
        # Create a heatmap using seaborn
        if mode == "steps" and min is not None and max is not None:
            sns.heatmap(matrix, annot=annot, cmap="RdYlGn_r", xticklabels=xtickslabels, yticklabels=yticklabels, cbar_kws={'label': 'RMSD (nm)'}, vmin=min, vmax=max)            
        else:
            sns.heatmap(matrix, annot=annot, cmap="RdYlGn_r", xticklabels=xtickslabels, yticklabels=yticklabels, cbar_kws={'label': 'RMSD (nm)'})
        
        plt.tight_layout()
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(filename, format="pdf", dpi=600)

    else:
        raise ValueError("Invalid mode. Choose either 'mutants' or 'steps'.")

def remove_temp_files(all = True, deep = False):
    
    if deep:
        dir_to_clean = [TEMP_DIR, MDP_DIR, TRAJ_DIR, RMSD_DIR, ENER_DIR, ENER_PLOT, TRAJ_PLOT, PCA_PLOT, CLUST_PLOT, MUTANT_TOP, OUT_DIR]
        for dirpath in dir_to_clean:
            for filename in os.listdir(dirpath):
                print(f"Removing file: {filename} from {dirpath}")
                os.remove(os.path.join(dirpath, filename))  
    if all:
        print("Cleaning up temporary files...")
        files = glob.glob(TEMP_DIR + sep + "*")
        for f in files:
            os.remove(f)
    
    files = glob.glob(cwd + sep + "*")
    for f in files: 
        filename = os.path.basename(f)
        if filename not in ["md.mdp", "nvt.mdp", "npt.mdp", "3mi.top", "16v.top", "tip3p.gro", "em.mdp", "ions.mdp", "md_c.mdp", "nvt_c.mdp", "npt_c.mdp", "em_c.mdp", "ions_c.mdp", "molecular_dynamics_proteinstability.py", "cgenff_charmm2gmx.py", "sort_mol2_bonds.pl", "visualize_trajectories.ipynb"]:
            if os.path.isfile(f):
                print(f"Removing file: {f}")
                os.remove(f)
    print("Temporary files cleaned up.")

# Main execution
if __name__ == "__main__":

    # Clean up temporary files
    deep_clean = False #keep it False if you want to keep the files between one mutant and another
    remove_temp_files(all=True, deep = deep_clean)
    
    start_time = datetime.now()
    print(f"Starting experiments at {start_time}")

    """
    # Run GROMACS experiments for each mutant without ligand 
    rmsd_matrix = compute_rmsd_matrix(ligand = None)
    plot_heatmap(rmsd_matrix, None, mode="mutants", mutants=mutants)
    """
    # Run GROMACS experiments for each mutant with ligand
    
    for ligand in ligands:
        print(f"\nProcessing ligand: {ligand}")
        rmsd_matrix = compute_rmsd_matrix(ligand)
        plot_heatmap(rmsd_matrix, ligand, mutant=None, mode="mutants", mutants=mutants)

    
    end_time = datetime.now()
    print(f"\nExperiments completed at {end_time}")
    print(f"Total runtime: {end_time - start_time}")

    # Clean up temporary files again
    #remove_temp_files(all = False, deep = False)