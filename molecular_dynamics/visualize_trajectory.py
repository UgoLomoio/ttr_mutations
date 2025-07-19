import os 
from pymol import cmd
import pymol

pymol.finish_launching(['pymol'])  # Launch PyMOL

cwd = os.getcwd()
sep = os.sep

topology_dir = cwd + sep + "topologies"
trajectory_dir = cwd + sep + "trajectories"
visual_trajectories_dir = cwd + sep + "visual_trajectories"
if not os.path.exists(visual_trajectories_dir):
    os.makedirs(visual_trajectories_dir)


def visualize_trajectory(top_file, traj_file, output_file):
    """
    Visualizes a molecular dynamics trajectory using PyMOL.

    Parameters:
    top_file (str): Path to the topology file.
    traj_file (str): Path to the trajectory file.
    output_file (str): Path to save the visualization output.
    """

    # Load the topology and trajectory files
    cmd.load(top_file, "topology")
    cmd.load_traj(traj_file)

    # Set visualization options
    cmd.show("cartoon", "topology")
    cmd.set("cartoon_transparency", 0.2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.hide("spheres")
    cmd.hide("lines")

    # Save the session as a PSE file
    cmd.save(output_file)
    print(f"Visualization saved to {output_file}")

if __name__ == "__main__":
    
    mutant = "v50m"
    ligand = "tafamidis"
    #ligand = "diflunisal"

    if ligand:
        top_file = f"{topology_dir}{sep}{mutant}-tetramer_{ligand}_md.gro"
        traj_file = f"{trajectory_dir}{sep}{mutant}-tetramer_{ligand}_fix.xtc"
        output_file = f"{visual_trajectories_dir}{sep}{mutant}-tetramer_{ligand}_visualization.pse"
    else:
        top_file = f"{topology_dir}{sep}{mutant}-tetramer_md.gro"
        traj_file = f"{trajectory_dir}{sep}{mutant}-tetramer_fix.xtc"
        output_file = f"{visual_trajectories_dir}{sep}{mutant}-tetramer_visualization.pse"
        
    visualize_trajectory(top_file, traj_file, output_file)