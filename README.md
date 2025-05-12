Requirements:

Operative System: Linux, MACOS, maybe Windows (?)
Python 3.9
Anaconda3


Installation:

1. Install Anaconda3 following this guide: 
    https://docs.conda.io/projects/conda/en/stable/user-guide/install/linux.html


2. Create a new conda enviroment called "ttr" with Python 3.10:

    conda create -n ttr python=3.9

    Before Installing AutoDock Vina on windows:

        Download boost from: https://www.boost.org/doc/libs/1_78_0/more/getting_started/windows.html
        Download swing from: https://www.swig.org/Doc1.3/Windows.html

        Unpack zip files inside the "ttr" conda enviroment directory (/home/$username$/anaconda3/envs/ttr)
        
    On Linux\MACOS:

        conda activate ttr
        conda install -c conda-forge numpy swig boost-cpp libboost sphinx sphinx_rtd_theme

3. Install PyMOL in the ttr enviroment using conda:

    conda activate ttr
    conda install -c conda-forge -c schrodinger pymol-bundle

4. Use pip to install other python reqiurements in the ttr enviroment:

    conda activate ttr
    pip install -r requirements.txt


5. If you want to generate mutated proteins using AlphaFold3 in a local machine, install AlphaFold3 inference code from: 

    git clone https://github.com/google-deepmind/alphafold3.git
    cd alphafold3
    pip install -r requirements.txt

(REMEMBER TO REQUEST ACCESS TO MODEL PARAMETERS: af.bin, then put this file inside the "/alphafold3/models/" directory)

6. To install DiffSBDD follow the guide of the official repository:

    cd docking
    git clone https://github.com/arneschneuing/DiffSBDD.git
    
    conda create -n DiffSBDD
    conda activate DiffSBDD
    conda install cudatoolkit=10.2 -c pytorch
    pip install torch
    conda install -c conda-forge pytorch-lightning
    conda install -c conda-forge wandb
    conda install -c conda-forge rdkit
    conda install -c conda-forge biopython
    conda install -c conda-forge imageio
    conda install -c anaconda scipy
    pip install torch-geometric
    pip install torch-scatter
    conda install -c conda-forge openbabel
    conda install seaborn

7. To install DiffDock:

    cd docking
    git clone https://github.com/gcorso/DiffDock.git

    Then create the enviroment using the scripts available in the README.md file of DiffDock:

        conda env create --file environment.yml
        conda activate diffdock

    Install Requirements:
    
        cd DiffDock
        conda install cudatoolkit -c nvidia
        conda install -c nvidia cuda-nvcc
        sudo ln -s /usr/local/cuda-8.0 /usr/local/cuda 
        pip install torch==1.13.1+cu117
        pip install -r requirements.txt

Order of execution to reproduce: 

Create Mutated structures: 

1. python generate_fasta.py
2. python generate_json.py 
3. Manual step: Use json_jobs to obtain Structures using AlphaFold Server (Beta). or automatic step using "python generate_structure.py" (needs alphafold3 inference code available on github and good hardware)
4. python popolate_cifs.py
5. python cif2pdb.py
6. python rename_pdbs.py

Align Single-point variants against wt using TM-Align:

1. All steps above
2. python tmalign_wt.py

Use DynaMut2 to compute impact of single mutations on protein stability: both monomer, dimer and tetramer

1. python create_tetramer_mutationlist.py
2. Manual step using dynamut2 server (https://biosig.lab.uq.edu.au/dynamut2/) or use dynamut2_inference.py

Docking with available ligands using DiffDock + AutoDock Vina: 

Install DiffDock-L following the steps from https://github.com/gcorso/DiffDock

1. cd docking
2. python dock_ligands.py or python docking_diffdock.py (uses NVIDIA NIM, requires a valid APIKEY)
3. python docking_plot.py

Create new ligands using DiffSBDD: (some changes to the original code are needed, follow the exceptions raised to resolve all the issues.)

1. cd docking
2. git clone https://github.com/arneschneuing/DiffSBDD.git
3. python generate_ligands.py


ESM2 Embedding:

1. python esm2_getembeddings_classify.py 

    This script uses NVIDIA NIM to compute embeddings of mutants. 

Molecular Dynamics:

Needs GROMACS installed, GPU compatibility is recommended.

1. cd molecular_dynamics
2. python molecular_dynamics_proteinstability.py (Needs Gromacs, performs automatically all the steps of topology preparation, computes trajectories and create plots)
3. visualize_trajectories.ipynb, a notebook that can be opened with jupyterlab that visualizes trajectories and creates final plots.


TMAlign:

Download TMalign_cpp.exe file from: https://zhanggroup.org/TM-align/TMalign_cpp.gz
Copy and paste the TMalign_cpp.exe file inside "tmalign_exe" directory
Now you can use our scripts: 

1. python compute_tmscores_all.py (or compute_tmscores_wt.py)
2. cd results
3. python plot_tmscore_all.py (or plot_tmscore_wt.py)


Protein Contact Networks:

0. Create a conda enviroment following the instructions in the README.md file contained inside "ProteinContactNetworks" directory. 
1. Use PCN-Miner command-line interface to extract leiden communities and compute node centrality such as betweenness, closeness, etc.
2. python plot_centrality_all.py
