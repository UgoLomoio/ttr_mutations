from molecular_dynamics_proteinstability import * 
import os
import subprocess
from pathlib import Path
from pymol import cmd
import pymol 

pymol.finish_launching(['pymol'])  # Launch PyMOL

def read_rmsf_from_pdb_bfactor(pdb_file):
    """
    Alternative: Read RMSF values directly from PDB B-factor column
    """
    rmsf_dict = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                resid = int(line[22:26].strip())
                bfactor = float(line[60:66].strip())
                rmsf_dict[resid] = bfactor
    return rmsf_dict

def set_bfactor_by_rmsf(object_name, rmsf_dict):
    """
    Set B-factor values to RMSF in PyMOL object
    """
    for chain, residues in rmsf_dict.items():
        for resid, data in residues.items():
            rmsf = data['rmsf']
            # Set b-factor for all atoms in the residue
            cmd.alter(f"{object_name} and chain {chain} and resi {resid}", f"b={rmsf}")
    cmd.rebuild()

def color_by_rmsf(object_name, palette='hot', selection='name CA'):
    """
    Color protein by RMSF using spectrum command
    """
    cmd.spectrum('b', palette, f'{object_name} and {selection}')


def correct_pbc_artifacts(xtc_file, tpr_file, output_file, index_file=None, 
                         center_group="Protein", output_group="System"):
    """
    Correct PBC artifacts from XTC trajectory files using GROMACS trjconv.
    
    Parameters:
    -----------
    xtc_file : str
        Path to input XTC trajectory file
    tpr_file : str
        Path to TPR structure file (required for PBC correction)
    output_file : str
        Path to output corrected XTC file
    index_file : str, optional
        Path to index file for custom group selection
    center_group : str
        Group to center on (default: "Protein")
    output_group : str
        Group to output (default: "System")
    
    Returns:
    --------
    bool
        True if correction successful, False otherwise
    """
    
    # Check if required files exist
    if not os.path.exists(xtc_file):
        print(f"Error: XTC file {xtc_file} not found")
        return False
    
    if not os.path.exists(tpr_file):
        print(f"Error: TPR file {tpr_file} not found")
        return False
    
    # Create temporary files for multi-step correction
    temp_dir = Path(output_file).parent
    temp_whole = temp_dir / "temp_whole.xtc"
    temp_nojump = temp_dir / "temp_nojump.xtc"
    
    try:
        # Step 1: Make molecules whole
        print("Step 1: Making molecules whole...")
        cmd1 = [
            "gmx", "trjconv",
            "-s", tpr_file,
            "-f", xtc_file,
            "-o", str(temp_whole),
            "-pbc", "whole"
        ]
        
        if index_file:
            cmd1.extend(["-n", index_file])
        
        result1 = subprocess.run(cmd1, input="0\n", text=True, 
                               capture_output=True, check=True)
        
        # Step 2: Remove jumps
        print("Step 2: Removing jumps...")
        cmd2 = [
            "gmx", "trjconv",
            "-s", tpr_file,
            "-f", str(temp_whole),
            "-o", str(temp_nojump),
            "-pbc", "nojump"
        ]
        
        if index_file:
            cmd2.extend(["-n", index_file])
        
        result2 = subprocess.run(cmd2, input="0\n", text=True, 
                               capture_output=True, check=True)
        
        # Step 3: Center and apply molecular PBC
        print("Step 3: Centering and applying molecular PBC...")
        cmd3 = [
            "gmx", "trjconv",
            "-s", tpr_file,
            "-f", str(temp_nojump),
            "-o", output_file,
            "-pbc", "mol",
            "-center",
            "-ur", "compact"
        ]
        
        if index_file:
            cmd3.extend(["-n", index_file])
        
        # Prepare input for group selection
        group_input = f"1\n0\n"  # Typically: 1 for Protein (center), 0 for System (output)
        
        result3 = subprocess.run(cmd3, input=group_input, text=True, 
                               capture_output=True, check=True)
        
        print(f"PBC correction completed successfully. Output: {output_file}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error during PBC correction: {e}")
        print(f"Command output: {e.stdout}")
        print(f"Command error: {e.stderr}")
        return False
    
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False
    
    finally:
        # Clean up temporary files
        for temp_file in [temp_whole, temp_nojump]:
            if temp_file.exists():
                temp_file.unlink()

def get_rmsf_dict_from_xvg(xvg_file, pdb_file):

    _, rmsf_values = read_xvg(xvg_file)
    rmsf_values = {key+1:value  for key, value in enumerate(rmsf_values)}
    rmsf_dict = {}
    with open(pdb_file, 'r') as f:  
        for line in f:
            if line.startswith('ATOM'):
                resid = int(line[22:26].strip())
                resname = line[17:20].strip()
                chain = line[21].strip()
                if chain not in rmsf_dict:
                    rmsf_dict[chain] = {}
                rmsf_dict[chain][resid] = {'resname': resname, 'rmsf': rmsf_values[resid]}
    return rmsf_dict

def plot_by_rmsf(protein_file, rmsf_file, trj_file=None, palette='red_white_blue'):

    cmd.load(protein_file, "ttr")
    if trj_file is not None:
        cmd.load_traj(trj_file)
        cmd.hide("lines")
        cmd.hide("spheres")
    
    cmd.show("cartoon", "ttr")  
    rmsf_values = get_rmsf_dict_from_xvg(rmsf_file, protein_file)
    set_bfactor_by_rmsf("ttr", rmsf_values)
    if "_" in palette:
        palette_list = palette.split("_")
    else:
        palette_list = None
    color_by_rmsf("ttr", palette=palette)
    min_ = min(rmsf_values[chain][resid]['rmsf'] for chain in rmsf_values for resid in rmsf_values[chain])
    max_ = max(rmsf_values[chain][resid]['rmsf'] for chain in rmsf_values for resid in rmsf_values[chain])
    
    if palette_list is not None:
        cmd.do("ramp_new colorbar, none, [{}, {}], {}".format(min_, max_, palette_list))
    else:
        cmd.do("ramp_new colorbar, none, [{}, {}], {}".format(min_, max_, palette))

if __name__ == "__main__":

    if None not in ligands:
        ligands.append(None)
    
    for mutant in mutants:
        print(f"Mutant: {mutant}")
        for ligand in ligands:
            print(f"Ligand: {ligand}")
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

            #Remove PBC artifacts 
            if ligand is not None:
                md_xtc_fix = f"{mutant}_{ligand}_fix.xtc"
            else:
                md_xtc_fix = f"{mutant}_fix.xtc"
                
            md_xtc_fix = os.path.join(TRAJ_DIR, md_xtc_fix)
            if not os.path.exists(md_xtc_fix):
                print(f"Correcting PBC artifacts for {md_xtc}...")
                _ = correct_pbc_artifacts(
                        xtc_file=md_xtc,
                        tpr_file=md_tpr,
                        output_file=md_xtc_fix
                )

            md_xtc = md_xtc_fix
            
            if not os.path.exists(out_file):
                print(f"Calculating RMSD for {mutant} with ligand {ligand}...")
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
            

            plot_rmsd(t, rmsd, filename)
            #plot_heatmap(matrix, ligand, mode = "steps", mutant=mutant)

            if ligand is None:
                rmsf_file = f"{RMSD_DIR}/{mutant}_rmsf.xvg"
                plotname = f"{TRAJ_PLOT}/{mutant}_rmsf.pdf"
            else:
                rmsf_file = f"{RMSD_DIR}/{mutant}_{ligand}_rmsf.xvg"
                plotname = f"{TRAJ_PLOT}/{mutant}_{ligand}_rmsf.pdf"

            # Check if RMSF file already exists
            if not os.path.exists(rmsf_file):
                print(f"Calculating RMSF for {mutant} with ligand {ligand}...")
                # Compute RMSF 
                cmd_rmsf = f"echo '1\n1' | gmx rmsf -s {md_tpr} -f {md_xtc} -o {rmsf_file}"
                subprocess.run(cmd_rmsf, shell=True, check=True, stdout=subprocess.DEVNULL)

            t_f, rmsf = read_xvg(filename=rmsf_file)
            plot_rmsf(t_f, rmsf, plotname)

            if ligand is None:
                protein_file = f"{PROTEIN_DIR}/{mutant}.pdb"
                trj_file = None
            else:
                protein_file = f"{MUTANT_TOP}/{mutant}_{ligand}_md.gro"
                trj_file = md_xtc

            print(f"Coloring structure by RMSF for {mutant}")
            plot_by_rmsf(protein_file, rmsf_file, trj_file=trj_file)

            if ligand is None:
                filename = f"{TEMP_DIR}/{mutant}_cluster.xpm"
            else:
                filename = f"{TEMP_DIR}/{mutant}_{ligand}_cluster.xpm"

            # Check if cluster file already exists
            if not os.path.exists(filename):
                print(f"Calculating clusters for {mutant} with ligand {ligand}...")
                # Perform clustering
                cmd_cluster = f"gmx cluster -s {md_gro} -dm {md_xpm} -o {filename} -tu ns"
                subprocess.run(cmd_cluster, shell=True, check=True, stdout=subprocess.DEVNULL)
            print(f"Plotting clusters for {mutant} with ligand {ligand} to {filename}...")
            plot_cluster(filename, mutant, ligand)

