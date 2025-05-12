import os
import sys
import pymol
from pymol import cmd
import pandas as pd

# Initialize PyMOL in quiet mode
pymol.finish_launching(['pymol', '-qc'])

# Define paths
cwd = os.getcwd()
docked_outputs = os.path.join(cwd, "AutoDock", "docked-outputs", "existing")  # Adjust path as needed
#ref_path = os.path.join(cwd, "AutoDock", "3tct_D_3MI.sdf")  # Your reference structure

# Function to convert PDBQT to PDB using PyMOL
def convert_pdbqt_to_pdb(pdbqt_file, output_dir):
    base_name = os.path.basename(pdbqt_file).replace('.pdbqt', '')
    pdb_file = os.path.join(output_dir, f"{base_name}.pdb")
    
    # Clear PyMOL
    cmd.delete('all')
    
    # Load PDBQT file
    cmd.load(pdbqt_file, 'docked')
    
    # Save as PDB
    cmd.save(pdb_file, 'docked')
    
    print(f"Converted {pdbqt_file} to {pdb_file}")
    return pdb_file

# Create output directory for converted PDBs
output_dir = os.path.join(cwd, "converted_pdbs")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load reference structure
cmd.delete('all')
ref_ligand_path = os.path.join(output_dir, 'reference_ligand.pdb')
cmd.load(ref_ligand_path, 'reference')

# Extract ligand from reference structure
cmd.select('ref_ligand', 'reference and organic')  # Adjust selection as needed
cmd.create('ref_ligand', 'ref_ligand')
cmd.save(os.path.join(output_dir, 'reference_ligand.pdb'), 'ref_ligand')

# Process all PDBQT files
results = []
files = [f for f in os.listdir(docked_outputs) if f.endswith('.pdbqt')]
print(f"Found {len(files)} PDBQT files")

for i, file in enumerate(files):
    pdbqt_path = os.path.join(docked_outputs, file)
    pdb_path = pdbqt_path.split(".")[0] + ".pdb"
    if not os.path.exists(pdb_path):  
        # Convert PDBQT to PDB using PyMOL
        pdb_path = convert_pdbqt_to_pdb(pdbqt_path, output_dir)
    
    # Calculate RMSD
    cmd.delete('all')
    cmd.load(ref_ligand_path, 'ref_ligand')
    cmd.load(pdb_path, 'docked_pose')
    
    # Align docked pose to reference ligand and get RMSD
    try:
        # First alignment without fitting (direct RMSD)
        rmsd_direct = cmd.rms_cur('docked_pose', 'ref_ligand')
        
        # Second alignment with fitting (minimized RMSD)
        alignment = cmd.align('docked_pose', 'ref_ligand')
        rmsd_fitted = alignment[0]  # First element is the RMSD after fitting
        
        results.append({
            'pose': i+1,
            'file': file,
            'rmsd_direct': rmsd_direct,
            'rmsd_fitted': rmsd_fitted
        })
        
        print(f"Pose {i+1} ({file}): Direct RMSD = {rmsd_direct:.4f} Å, Fitted RMSD = {rmsd_fitted:.4f} Å")
    except Exception as e:
        print(f"Error calculating RMSD for {file}: {e}")

# Create DataFrame and save results
df = pd.DataFrame(results)
results_path = os.path.join(cwd, "rmsd_results.csv")
df.to_csv(results_path, index=False)
print(f"Results saved to {results_path}")

# Identify successful poses (typically RMSD < 2.0 Å is considered successful)
successful_poses = [r for r in results if r["rmsd_fitted"] < 2.0]
print(f"\nNumber of successful poses (RMSD < 2.0 Å): {len(successful_poses)}")

# Sort results by RMSD
sorted_results = sorted(results, key=lambda x: x["rmsd_fitted"])

# Print top 5 poses
print("\nTop 5 poses by RMSD:")
for i, result in enumerate(sorted_results[:5]):
    print(f"Rank {i+1}: {result['file']} - RMSD = {result['rmsd_fitted']:.4f} Å")

"""
# Create a PyMOL session with the best pose aligned to the reference
if sorted_results:
    best_pose = sorted_results[0]
    best_pose_file = best_pose['file']
    best_pose_path = os.path.join(output_dir, best_pose_file.replace('.pdbqt', '.pdb'))
    
    cmd.delete('all')
    cmd.load(ref_ligand_path, 'ref_ligand')
    cmd.load(best_pose_path, 'best_pose')
    cmd.align('best_pose', 'ref_ligand')
    
    # Color the molecules
    cmd.color('green', 'ref_ligand')
    cmd.color('cyan', 'best_pose')
    
    # Save session
    session_path = os.path.join(cwd, "best_pose_alignment.pse")
    cmd.save(session_path)
    print(f"PyMOL session with best pose alignment saved to {session_path}")
"""

# Clean up PyMOL
cmd.quit()
"""
    # Optional: Analyze protein-ligand interactions for best pose
    # Requires additional libraries like PLIP
    try:
        from plip.structure.preparation import PDBComplex
        from plip.exchange.report import BindingSiteReport
        
        # For the best pose (lowest RMSD)
        best_pose_idx = min(range(len(results)), key=lambda i: results[i]["rmsd_minimized"])
        best_pose_file = f"best_pose_{best_pose_idx+1}.pdb"  # You need to save this file first
        
        # Analyze interactions
        complex_struct = PDBComplex()
        complex_struct.load_pdb(best_pose_file)
        for ligand in complex_struct.ligands:
            complex_struct.characterize_complex(ligand)
        
        # Print interaction summary
        for key, site in sorted(complex_struct.interaction_sets.items()):
            binding_site = BindingSiteReport(site)
            print(f"\nInteractions for best pose (RMSD = {results[best_pose_idx]['rmsd_minimized']:.4f} Å):")
            print(f"Hydrophobic interactions: {len(binding_site.hydrophobic_info)}")
            print(f"Hydrogen bonds: {len(binding_site.hbond_info)}")
            print(f"Salt bridges: {len(binding_site.saltbridge_info)}")
            print(f"Pi-stacking: {len(binding_site.pistacking_info)}")
            print(f"Pi-cation interactions: {len(binding_site.pication_info)}")
            print(f"Halogen bonds: {len(binding_site.halogen_info)}")
    except ImportError:
        print("\nPLIP not installed. Skipping interaction analysis.")


    # Load structures
    ref = mda.Universe("experimental_tafamidis.pdb")
    docked = mda.Universe("docked_pose.pdb")

    # Select ligand atoms
    ref_ligand = ref.select_atoms("resname TAF")  # Replace TAF with your ligand residue name
    docked_ligand = docked.select_atoms("resname TAF")

    # Align the docked pose to the reference
    alignment = mda.analysis.align.alignto(docked_ligand, ref_ligand)
    print(f"RMSD after alignment: {alignment[1]:.4f} Å")

    # Visualize
    view = nv.show_mdanalysis(ref)
    view.add_component(docked_ligand.positions, default=False)
    view.add_representation('licorice', selection='all', component=0, color='green')
    view.add_representation('licorice', selection='all', component=1, color='red')
    view.center()
    view
    """