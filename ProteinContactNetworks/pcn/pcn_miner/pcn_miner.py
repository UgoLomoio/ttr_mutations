#IMPORT LIBRARIES
import numpy as np
from networkx import from_numpy_array, to_numpy_array
from urllib.request import urlretrieve
from scipy.spatial.distance import euclidean
from operator import itemgetter
from sklearn import metrics
from scipy.linalg import eigh, eig
import os
from sys import platform
from multiprocessing import Pool
import datetime
import itertools
import matplotlib.pyplot as plt

#tk GUI progress bar
import tkinter as tk
from tkinter import ttk

#centrality measures
from networkx.algorithms.centrality import degree_centrality, eigenvector_centrality, closeness_centrality, betweenness_centrality, betweenness_centrality_subset

#embedding algorithms
from gem.embedding.hope     import HOPE
from gem.embedding.lap      import LaplacianEigenmaps
from node2vec import Node2Vec

#algorithms for clustering
from sklearn.cluster import SpectralClustering, KMeans
from fcmeans import FCM

#community detection algorithms

from cdlib.algorithms import infomap as infomap_cdlib
from cdlib.algorithms import spinglass as spinglass_cdlib
from cdlib.algorithms import leiden as leiden_cdlib
from cdlib.algorithms import louvain as louvain_cdlib
from cdlib.algorithms import walktrap as walktrap_cdlib

#from networkx.algorithms.community.centrality import girvan_newman as girvan_newman_
from networkx.algorithms.community.asyn_fluid import asyn_fluidc as asyn_fluidc_
from networkx.algorithms.community import greedy_modularity_communities

#handle different path separators
if platform == "linux" or platform == "linux2":
    # linux
    add_slash_to_path = '/'
elif platform == "darwin":
    # OS X
    add_slash_to_path = '/'
elif platform == "win32":
    # Windows...
    add_slash_to_path = '\\'

def checkIfFilesExists(files, initial_choice, proteins_path, adj_path = None, comp_adj_fr = None, window = None):
    """
        Given a list of 'pdb' files or 'adj' files the algorithm:
            in case of initial_choice='adj' will check if the adjacency cfx file of the given protein really exists in the "adj_path" and it will check if the pdb file exists in the protein_path
            in case of initial_choice='pdb' will check if the pdb file exists in the protein_path
        Parameters:
            files: list of string, list of files
            initial_choice: string, describe the type of the files: 'pdb' if they are pdb files or 'adj' if they are adjacency matrix files
            protein_path:
            adj_path: default None, file path for adjacency matrix files, it has to be not None only if initial_choice='adj'.
            comp_adj_fr: tk.Frame, is the frame of the GUI that contains the progress bar.
            window: tk.Tk, is the window of the GUI.
        Returns: None
    """

    not_existing_pdb_files = []
    all_pdb_files_exists = True

    for file in files:

        if initial_choice=="pdb":
            p_name = file[:(len(file)-4)]
        elif initial_choice=="adj":
            filename_splitted = (file.split(".txt"))[0].split("_") #adj = 6vxx_adj_mat_min_max.txt
            p_name = filename_splitted[0]

        pdb_path = "{}{}.pdb".format(proteins_path, p_name)
        if((not os.path.isfile(pdb_path))):
            all_pdb_files_exists = False
            not_existing_pdb_files.append(file)

    if(not all_pdb_files_exists):
        for file in not_existing_pdb_files:
            if initial_choice=="pdb":
                p_name = file[:(len(file)-4)]
            elif initial_choice=="adj":
                filename_splitted = (file.split(".txt"))[0].split("_") #adj = 6vxx_adj_mat_min_max.txt
                p_name = filename_splitted[0]
            print(("protein {} pdb file missing, fetching on PDB database...").format(p_name))
            urlretrieve("http://files.rcsb.org/download/{}.pdb".format(p_name), "{}{}.pdb".format(proteins_path, p_name))
            print(("protein {} fetched successfully").format(p_name))
        all_pdb_file_exists = True

    if(initial_choice == "adj"):
        not_existing_adj_files = []
        all_adj_files_exists = True
        for file in files:

            file_path = "{}{}".format(adj_path, file)
            if((not os.path.isfile(file_path))):
                all_adj_files_exists = False
                not_existing_adj_files.append(file)

        if (not all_adj_files_exists):
            for filename in not_existing_adj_files:
                filename_splitted = (filename.split(".txt"))[0].split("_") #adj = 6vxx_adj_mat_min_max.txt
                p_name = filename_splitted[0]
                min_ = float (filename_splitted[3])
                max_ = float (filename_splitted[4])
                print(("protein {} adj matrix missing...").format(p_name))
                protein_path = proteins_path+p_name+".pdb"
                atoms = readPDBFile(protein_path)
                coordinates = getResidueCoordinates(atoms)
                dict_residue_name = associateResidueName(coordinates)
                residue_names = np.array(list (dict_residue_name.items()))

                if comp_adj_fr is not None and window is not None:
                    computing_A_label = tk.Label(comp_adj_fr, text="Missing adjacency matrix for protein {}. Computing now... (This may take time)".format(p_name), bg = "dodger blue")
                    computing_A_label.pack()
                    comp_adj_fr.pack()
                    window.update()

                print("computing adjacency matrix with thresholds: min = {} and max = {} ... (This may take time)".format(min_, max_))
                output_path = os.path.abspath(os.path.join(adj_path, os.pardir))+add_slash_to_path

                #parallel computation, TODO: TEST PARALLEL COMPUTATION ON MAC AND THEN UNCOMMENT THIS
                #A =  adjacent_matrix(output_path, coordinates, p_name, min_, max_, comp_adj_fr, window)
                #non parallel computation
                A = adjacent_matrix_nonparallel(output_path, coordinates, p_name, min_, max_, comp_adj_fr, window)

            all_adj_file_exists = True

def readPDBFile(pdbFilePath):
    """
    Read the pdb file.
    Parameters:
        pdbFilePath: string, is the complete PDB file path to read.
    Returns:
        atoms: np.array, contains all the informations contained in the 'ATOM' key of the pdb files.
    """
    datetime_object = datetime.datetime.now()
    print('Start Reading PDB'+'\n')
    print(datetime_object)
    atoms = []
    with open(pdbFilePath) as pdbfile:
        for line in pdbfile:
            if line[:4] == 'ATOM':
                #split the line
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54], line[56:61], line[62:66]]
                atoms.append(splitted_line)
    datetime_object = datetime.datetime.now()
    print('End Reading PDB'+'\n')
    print(datetime_object)
    return np.array(atoms)

def getResidueCoordinates(atoms):
    """
    Compute the distances between the alpha carbons of the amino acids.
    Parameters:
        atoms: np.array, contains all the informations contained in the 'ATOM' key of the pdb files.
    Returns:
        coordinates: np.array, contains the euclidean distance between the alpha carbon of the amino acids of the proteins.
    """
    coordinates = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS',
                   'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                   'TRP', 'TYR']

    dropped = []

    last_residue_num = 0

    for i, atom in enumerate(atoms):

        residue_name = atom[3]
        residue_chain = atom[4]
        residue_num = int (atom[5])

        if residue_name in residues_list:

            if (residue_num != last_residue_num):

                if (atom[2].replace(" ","")=="CA"):
                    cord_C = atom[6:9]
                    coordinates.append([residue_name + str (residue_num) + " " + residue_chain, cord_C])
                    last_residue_num = residue_num

        else:
          dropped.append(residue_name)

    return np.array(coordinates, dtype=object)

def associateResidueName(coordinates):
    """
    Associate the protein residues names to the networkx node ID.
    Parameters:
        coordinates, np.array, contains the list of residues names
    Returns:
        dict_residue_name: dictionary {residue_nxID: residue_name}
    """
    dict_residue_name = dict()
    for i in range(coordinates.shape[0]):
        dict_residue_name[str (i)] = coordinates[i, 0]
    return dict_residue_name

def getResiduesSequence(pbdFilePath):
    """
    Read the amino acids sequence from the pdb file.
    Parameters:
        pdbFilePath: string, is the complete PDB file path to read.
    Returns:
        seq_res: np.array, contains all the informations contained in the 'SEQRES' key of the pdb files (sequence of the amino acids).
    """
    seq_res = []
    residues_list = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']

    with open(pbdFilePath) as pdbfile:
        for line in pdbfile:
            if line[:6]=='SEQRES':
                splitted_line = [line[19:22], line[23:26], line[27:30], line[31:34], line[35:38],
                                line[39:42], line[43:46], line[47:50], line[51:54], line[55:58],
                                line[59:62], line[63:66], line[67:70]]
                for residue in splitted_line:
                    if residue in residues_list:
                        seq_res.append(residue)

    return np.array(seq_res)

def read_adj_mat(adj_filepath, p, min_, max_):
    """
    Read the adjacency matrix file.
    Parameters:
        adj_filepath: string, is the complete adjacency matrix file path to read.
        p: string with len equals to 4, is the protein pdb code.
        min_: float, is the minimum threshold distance for extract only the non-covalent interactions between amino acids.
        max_: float, is the maximum threshold distance for extract only the significant interactions between amino acids.
    Returns: None
    """
    if (os.path.isfile("{}{}_adj_mat_{}_{}.txt".format(adj_filepath, p, min_, max_))):
        adj = np.loadtxt("{}{}_adj_mat_{}_{}.txt".format(adj_filepath, p, min_, max_))
        return adj
    else:
        raise Exception("Adj matrix for protein {} doesn't exists.".format(p))


def printProgressBar (iteration, total):

    length = 100
    fill = '█'
    printEnd = '\r'
    percent = ("{0:." + str(2) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f"\r |{bar}| Current progress: {percent}%", end = printEnd)

    if percent == 100:
        print()


def adjacent_matrix(output_path, coordinates, p, min_=4, max_=8, comp_adj_fr=None, window = None):
    from sklearn.metrics import pairwise_distances
    import time
    """
    Parallel computation of the adjacency matrix.
    Parameters:
        output_path: string, is the output file path.
        coordinates: np.array, contains the euclidean distance between the alpha carbon of the amino acids of the proteins.
        p: string with len equals to 4, is the protein pdb code.
        min_: float, is the minimum threshold distance for extract only the non covalent interactions between amino acids.
        max_: float, is the maximum threshold distance for extract only the significant interactions between amino acids.
        comp_adj_fr: tk.Frame, is the frame of the GUI that contains the progress bar.
        window: tk.Tk, is the window of the GUI.
    Returns: adj
    """
    start = time.time()
    n = coordinates.shape[0]
    cords = np.hstack(coordinates[:,1]).reshape(n,3)
    adj = np.zeros((n,n))
    #d = np.zeros((n,n), dtype=float)
    d = pairwise_distances(X = cords, n_jobs = -1)

    if comp_adj_fr is not None:
        pb = ttk.Progressbar(comp_adj_fr, orient="horizontal", mode="determinate", length=100)
        pb.pack()
        pb["value"] = 0
        label = tk.Label(comp_adj_fr, text="Current progress {}%".format(pb["value"]))
        label.pack()
        window.update()
    else:
        value = 0
        printProgressBar(value, n)

    for i in range(n):
        for j in range(i):
            if ((d[i][j]>min_) and (d[i][j]<max_)):
                adj[i][j] = 1
                adj[j][i] = 1

        if comp_adj_fr is not None:
            pb["value"] = round(((i + 1) / n) * 100, 2)
            label['text'] = "Current progress {}%".format(pb["value"])
            pb.pack()
            label.pack()
            window.update()
        else:
            printProgressBar(i + 1, n)

    if comp_adj_fr is not None:
        pb["value"] = round(((i + 1) / n) * 100, 2)
        label['text'] = "Current progress {}%".format(pb["value"])
        pb.pack()
        label.pack()
        window.update()
    else:
        printProgressBar(i + 1, n)

    end = time.time()
    print("Time for parallel PCN computation of protein {}: {} s".format(p, (end-start)))

    if not os.path.exists("{}Adj".format(output_path)):
        os.makedirs("{}Adj".format(output_path))
    np.savetxt("{}Adj{}{}_adj_mat_{}_{}.txt".format(output_path, add_slash_to_path, p, min_, max_), adj, fmt='%.2f')
    print("saved adj matrix")

    return adj

def adjacent_matrix_nonparallel(output_path, coordinates, p, min_=4, max_=8, comp_adj_fr=None, window = None):
    """
    Non parallel computation the adjacency matrix.
    Parameters:
        output_path: string, is the output file path.
        coordinates: np.array, contains the euclidean distance between the alpha carbon of the amino acids of the proteins.
        p: string with len equals to 4, is the protein pdb code.
        min_: float, is the minimum threshold distance for extract only the non-covalent interactions between amino acids.
        max_: float, is the maximum threshold distance for extract only the significant interactions between amino acids.
        comp_adj_fr: tk.Frame, is the frame of the GUI that contains the progress bar.
        window: tk.Tk, is the window of the GUI.
    Returns: None
    """
    n = coordinates.shape[0]
    adj = np.zeros((n,n))
    d = np.zeros((n,n), dtype=float)
    edge_list = []

    if comp_adj_fr is not None:
        pb = ttk.Progressbar(comp_adj_fr, orient="horizontal", mode = "determinate", length = 100)
        pb.pack()
        pb["value"] = 0
        label = tk.Label(comp_adj_fr, text = "Current progress {}%".format(pb["value"]))
        label.pack()
        window.update()
    else:

        value = 0
        printProgressBar(value, n)

    for i in range(n):
        for j in range(n):
            if (i!=j):
                p1 = np.array(coordinates[i, 1], dtype=float) #CORD_C
                p2 = np.array(coordinates[j, 1], dtype=float) #CORD_C
                d[i][j] = euclidean(p1, p2)
                if ((d[i][j]>min_) and (d[i][j]<max_)):
                    adj[i][j] = 1
                    adj[j][i] = 1
                    if (([j, i] not in edge_list) and ([i, j] not in edge_list)):
                        edge_list.append([i, j])

        if comp_adj_fr is not None:
            pb["value"]= round(((i+1)/n)*100, 2)
            label['text'] = "Current progress {}%".format(pb["value"])
            pb.pack()
            label.pack()
            window.update()
        else:
            printProgressBar(i + 1, n)

    if comp_adj_fr is not None:
        pb["value"]= round(((i+1)/n)*100, 2)
        label['text'] = "Current progress {}%".format(pb["value"])
        pb.pack()
        label.pack()
        window.update()
    else:
        printProgressBar(i + 1, n)


    #save Edgelists, Adj matrixs and Distances matrixs
    """
    if not os.path.exists("{}Distances".format(output_path)):
        os.makedirs("{}Distances".format(output_path))
    np.savetxt("{}Distances{}{}_distances.txt".format(output_path, add_slash_to_path, p), d, fmt='%.2f')
    print("saved distances matrix")

    if not os.path.exists("{}Edgelists".format(output_path)):
        os.makedirs("{}Edgelists".format(output_path))

    np.savetxt("{}Edgelists{}{}_edgelist_{}_{}.csv".format(output_path, add_slash_to_path, p, min_, max_), np.array(edge_list), fmt='%.2f')
    print("saved edge list")
    """
    if not os.path.exists("{}Adj".format(output_path)):
        os.makedirs("{}Adj".format(output_path))
    np.savetxt("{}Adj{}{}_adj_mat_{}_{}.txt".format(output_path, add_slash_to_path, p, min_, max_), adj, fmt='%.2f')
    print("saved adj matrix")

    return adj

def save_centralities(output_path, centralities, p_name, method = None):
    """
    Save the node centralities as txt file in the output directory.
    Parameters:
        output_path: string, path to use when save the centralities.
        centralities: dict {node: centrality}, node 'method' centralities.
        p_name: string, pdb code of the protein to study.
        method: string, the centrality measure algorithms used.
    Returns: None
    """

    centrality_output_path = "{}Centralities{}{}{}Txt{}".format(output_path, add_slash_to_path, method, add_slash_to_path, add_slash_to_path)

    if (not os.path.exists(centrality_output_path)):

        os.makedirs(centrality_output_path)

    f = open("{}{}_{}.txt".format(centrality_output_path, p_name, method),"w")
    f.write(str(centralities))
    f.close()


def save_labels(output_path, labels, residue_names, p_name, method=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Save the node centralities as txt file in the output directory.
    Parameters:
        output_path: string, path to use when save communities/clusters.
        labels : np.array, list of clusters/communities.
        residue_names: np.array, list of the residues names of the protein.
        p_name: string, pdb code of the protein to study.
        method: string, the spectral clustering / embedding+clustering / community detection algorithm used to create the partition.
        d: int, default None, dimension for the embedding
        beta: float, default None, decay factor for HOPE embedding
        walk_len: int, default None, length of the random walks used in node2vec embedding.
        num_walks: int, default None, number of random walks each node computed in node2vec embedding.
    Returns:
        dict_node_cluster_1: dict {node: label}, 'method' extracted clusters/communities.
    """
    supported_methods_clustering = ["unnorm_ssc", "norm_ssc", "unnorm_hsc", "norm_hsc", "hsc_shimalik", "ssc_shimalik", "skl_spectral_clustering"]
    supported_methods_embeddings = [
                                       "fuzzycmeans_hope", "kmeans_hope", "fuzzycmeans_laplacianeigenmaps", "kmeans_laplacianeigenmaps" ,
                                       "fuzzycmeans_node2vec", "kmeans_node2vec"
                                   ]
    supported_methods_communities = ["louvain", "leiden", "walktrap", "asyn_fluidc", "greedy_modularity", "infomap", "spinglass"]

    if platform == "win32":
        # Windows...
        supported_methods_communities.remove("infomap")

    if ((method in supported_methods_clustering)|(method in supported_methods_communities)|(method in supported_methods_embeddings)):

        if method in supported_methods_clustering:
            name = "Clusters"
        if method in supported_methods_communities:
            name = "Communities"
        if method in supported_methods_embeddings:
            name = "ClustersEmbeddings"

    residue_names_0 = np.array(residue_names[:, 0], dtype = int)
    residue_names_1 = np.array(residue_names[:, 1], dtype = str)

    labels_u, counts = np.unique(labels, return_counts=True)

    for i in range(len(labels_u)):

        print(labels_u[i], counts[i])

    dict_node_cluster_0 = dict ()
    for i, label in enumerate(labels):

      dict_node_cluster_0[str (residue_names_0[i])] = label

    residue_clusters = np.array(list (dict_node_cluster_0.items()), dtype=object)
    residue_names_cluster = dict ()

    k = int (max(labels)) + 1

    for label in range(k):

        temp = []

        for (resi, res_name) in residue_names:
            if ((residue_names[int (resi), 0]==residue_clusters[int (resi), 0]) and (str (int ((residue_clusters[int (resi), 1]))) == str (label))):

                temp.append(res_name)

        print(len(temp))
        residue_names_cluster[label] = temp
        print(f"{name} {label}: ", residue_names_cluster[label])

    summary_output_path = "{}{}{}{}".format(output_path, method, add_slash_to_path, "Summary")

    if (not os.path.exists(summary_output_path)):

        os.makedirs(summary_output_path)

    f = open("{}{}{}_{}_{}_k{}.txt".format(summary_output_path, add_slash_to_path, p_name, "{}_Summary".format(name), method, k),"w")

    for label in range(k):

        f.write("{} {}: {} ".format(name, label, residue_names_cluster[label]))
        f.write("\r\n")

    f.close()

    dict_node_cluster_1 = dict ()

    for i, label in enumerate(labels):

        dict_node_cluster_1[residue_names_1[i]] = int (label)

    method_output_path = "{}{}{}{}".format(output_path, method, add_slash_to_path, name)

    if (not os.path.exists(method_output_path)):

        os.makedirs(method_output_path)

    if (name == "Clusters"):

        f = open("{}{}{}_{}_{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, k),"w")
        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1

    elif (name == "ClustersEmbeddings"):

        if beta is not None:

            f = open("{}{}{}_{}_{}_d{}_beta{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, d, beta, k),"w")
        elif num_walks is not None and walk_len is not None:
            f = open("{}{}{}_{}_{}_d{}_wl{}_nw{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, d, walk_len, num_walks, k), "w")
        else:

            f = open("{}{}{}_{}_{}_d{}_k{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, d, k),"w")

        f.write(str (dict_node_cluster_1))
        f.close()
        return dict_node_cluster_1

    elif (name == "Communities"):

        f = open("{}{}{}_{}_{}_ncoms{}.txt".format(method_output_path, add_slash_to_path, p_name, name, method, k),"w")
        f.write(str (dict_node_cluster_1))
        f.close()

        return dict_node_cluster_1

    else:
        raise Exception ("method {} not supported".format(method))

def save_part_coef(output_path, part_coefs, p_name, method, k):
    """
    Save the nodes partecipation coefficients as txt file in the output directory.
    Parameters:
        output_path: string, path to use when save the node partecipation coefficients.
        part_coefs: dict {node: part_coef}, node partecipation coefficient.
        p_name: string, pdb code of the protein to study.
        method: string, the spectral clustering/ embedding+clusterin / community detection algorithm used to create the partition.
        k: int, number of clusters/communities extracted.
    Returns: None
    """
    method_output_path = "{}{}{}Part_coefs_txt".format(output_path, method, add_slash_to_path)
    if (not os.path.exists(method_output_path)):
        os.makedirs(method_output_path)
    f = open("{}{}{}_{}_part_coefs_k{}.txt".format(method_output_path, add_slash_to_path, p_name, method, k),"w")
    f.write(str (part_coefs))
    f.close()

#COMMUNITY EXTRACTION
def extract_labels_from_coms(num_nodes, coms, algorithm_name):

    if (algorithm_name != "Asyn FluidC"):
        n_coms = len(coms)
        print("number of {} communities: {}".format(algorithm_name, n_coms))

    labels = np.zeros((num_nodes, 1))
    dict_node_algorithm_com = dict ()

    for label, com in enumerate(coms):

        for i, node in enumerate(com):

            node =  int (float (node))
            dict_node_algorithm_com[node] = label
            labels[node] = label

    return labels

def louvain(G):
    """
    Implementation from the cdlib library for the community extraction Louvain algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start Louvain'+'\n')
    print(datetime_object)
    coms = louvain_cdlib(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Louvain")
    datetime_object = datetime.datetime.now()
    print('End Louvain'+'\n')
    print(datetime_object)
    return labels

def leiden(G):
    """
    Implementation from the cdlib library for the community extraction Leiden algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start Leiden'+'\n')
    print(datetime_object)
    coms = leiden_cdlib(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Leiden")
    datetime_object = datetime.datetime.now()
    print('End Leiden'+'\n')
    print(datetime_object)
    return labels

def walktrap(G):
    """
    Implementation from the cdlib library for the community extraction Walktrap algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start walktrap'+'\n')
    print(datetime_object)
    coms = walktrap_cdlib(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Walktrap")
    datetime_object = datetime.datetime.now()
    print('End walktrap'+'\n')
    print(datetime_object)
    return labels

def greedy_modularity(G):
    """
    Implementation from the cdlib library for the community extraction Greedy Modularity algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    coms = greedy_modularity_communities(G)
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Greedy Modularity")

    return labels

def infomap(G):
    """
    Implementation from the cdlib library for the community extraction Infomap algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start Infomap'+'\n')
    print(datetime_object)
    coms = infomap_cdlib(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Infomap")
    datetime_object = datetime.datetime.now()
    print('End Infomap'+'\n')
    print(datetime_object)
    return labels

def asyn_fluidc(G, k):
    """
    Implementation from the cdlib library for the community extraction Async FluidC algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start Async'+'\n')
    print(datetime_object)
    coms = asyn_fluidc_(G, k)
    num_nodes = G.number_of_nodes()
    print("number of Asyn FluidC communities: {}".format(k))
    labels = extract_labels_from_coms(num_nodes, coms, "Asyn FluidC")
    datetime_object = datetime.datetime.now()
    print('End Async'+'\n')
    print(datetime_object)
    return labels

def spinglass(G):
    """
    Implementation from the cdlib library for the community extraction Spinglass algorithm.

    Parameters:
        G: networkx.graph, the graph (the PCN) you want to extract the communities
    Returns:
        labels: numpy.array, the extracted communities
    """
    datetime_object = datetime.datetime.now()
    print('Start Spin Glass'+'\n')
    print(datetime_object)
    coms = spinglass_cdlib(G).communities
    num_nodes = G.number_of_nodes()
    labels = extract_labels_from_coms(num_nodes, coms, "Spin Glass")
    datetime_object = datetime.datetime.now()
    print('End Spin Glass'+'\n')
    print(datetime_object)

    return labels

#SPECTRAL CLUSTERING

def degree_matrix(A):
    """
    Compute the degree matrix of a graph
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
    Returns:
        D: numpy.array, the degree matrix of the graph.
    """
    n = A.shape[0]
    diag = np.zeros((n, n))
    rows, cols = A.nonzero()

    for row, col in zip(rows, cols):

        diag[row, row] += 1

    return diag

def compute_laplacian_matrix(A):
    """
    Compute the unnormalized laplacian matrix of a graph
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
    Returns:
        L: numpy.array, the unnormalized laplacian matrix of the graph.
    """
    D = degree_matrix(A)
    L = D-A
    return L

def compute_normalized_laplacian(A):
    """
    Compute the normalized laplacian matrix of a graph
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
    Returns:
        L: numpy.array, the normalized laplacian matrix of the graph.
    """
    D = degree_matrix(A)
    D_inv_sqrt = np.linalg.inv(np.sqrt(D))
    L_sym = 1 - np.dot(D_inv_sqrt, A).dot(D_inv_sqrt)
    return L_sym

def computeSortEigens(mat):
    """
    Compute and sort eigenvalues and eigenvectors of a given matrix.
    Eigenvalues are sorted in ascending order, then the associated eigenvectors are sorted following the eigenvalues sorting indices.
    Parameters:
        mat: numpy.array, is the matrix (adjacency or laplacian matrix) used for the computation of eigenvalues and eigenvectors.
    Returns:
        sortedEigenvalues: numpy.array, sorted eigenvalues of the matrix.
        sortedEigenvectors: numpy.array, sorted eigenvectors of the matrix.
    """
    eigenvalues, eigenvectors = eig(mat)

    idx = eigenvalues.argsort()
    sortedEigenvalues = eigenvalues[idx].real
    sortedEigenvectors = eigenvectors[:,idx].real

    return sortedEigenvalues, sortedEigenvectors

def computeBestK(eigenvalues, n_k = 1):
    """
    Compute the best number of clusters for a Spectral Clustering approach using the Max Eigengap Method.
    Parameters:
        eigenvalues: numpy.array or array like, is the list of eigenvalues of the Laplacian matrix.
        n_k: int, default equals to 1, is the number of best ks to return.
    Returns:
        best_ks: list of ints, is the list of number of best ks. Is length is equal to n_k.
    """
    sortedEigenvalues = sorted(eigenvalues)[::-1]
    max_eigengap = 0
    n = len(sortedEigenvalues)
    all_k = np.arange(0, n)
    best_k = 0
    all_eigengaps = np.zeros((n))

    for index, k in enumerate(all_k):

        if (k!=0):

            eigengap = abs(sortedEigenvalues[k]-sortedEigenvalues[k-1])
            all_eigengaps[k] = eigengap

            if (eigengap > max_eigengap):

                best_k = k
                max_eigengap = eigengap

        else: continue

    idx = all_eigengaps.argsort()[::-1]
    real_best_ks = all_k[idx]

    if ((best_k > 60)|(best_k==1)):

        best_k = 0
        for k in real_best_ks:

            if ((k<60)&(k>1)):

                best_k = k
                break

    best_ks = real_best_ks[(real_best_ks<60) & (real_best_ks>1)][:n_k]
    print("Best k: ", best_k)
    #print("Real Best {} ks: ".format(n_k), real_best_ks)
    print("Choosen Best {} ks: ".format(n_k), best_ks)

    return best_ks


def hardSpectralClustering(A, n_clusters = None, norm=False, embedding=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Hard Spectral Clustering Algorithm.
    Compute the eigenvectors and then use on them KMeans for clustering.
    In case of embedding, it uses the embedding on the graph and then use KMeans for extract the clusters.

    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
        n_clusters: int, default None. The number of clusters
        norm: boolean, default False. If is True the algorithm will compute the normalized laplacian matrix, otherwise the algorithm will compute the unnormalize laplacian matrix.
        embedding: string, default None. If None, no embeddings algorithms are applied. Outherwise, a supported embedding algorithm is applied on the eigenvectors of the laplacian matrix.
        d: int, default None, dimension for the embedding
        beta: float, default None, decay factor for HOPE embedding
        walk_len: int, default None, length of the random walks used in node2vec embedding.
        num_walks: int, default None, number of random walks each node computed in node2vec embedding.
    Returns:
        labels: extracted clusters
    """
    supported_embeddings = ["HOPE", "LaplacianEigenmaps", "Node2Vec"]

    if embedding is None:

        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)

        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)             #sorted eigenvalues/eigenvectors
        train = sortedEigenvectors[:, :n_clusters]

    else:

        if (embedding in supported_embeddings):

            G = from_numpy_array(A)
            if (embedding == "HOPE"): #embedding dimension (d) and decay factor (beta) as inputs
                model = HOPE(d=d, beta=beta)
            elif (embedding == "LaplacianEigenmaps"):
                model = LaplacianEigenmaps(d=d)
            elif (embedding == "Node2Vec"):
                model = Node2Vec(G, dimensions=d, walk_length=walk_len, num_walks=num_walks, workers = 1)

            if embedding == "Node2Vec":
                model = model.fit()
                train = model.wv.vectors
                train = np.array(train)
                if (train.shape[0] < A.shape[0]):
                    raise Exception("Parameters 'num_walks' and 'walk_length' for Node2Vec embedding are too small. Can't train the model.")
            else:
                train, t = model.learn_embedding(graph=G)
                train = np.array(train)

        else: raise Exception ("embedding {} not supported".format(embedding))

    km = KMeans(n_clusters=n_clusters).fit(train)
    labels = km.labels_

    return labels

def softSpectralClustering(A, n_clusters = None, norm=False, embedding = None,  d=None, beta=None, walk_len=None, num_walks=None):
    """
    Soft Spectral Clustering Algorithm.
    Compute the eigenvectors and then use on them Fuzzy C-Means for clustering.
    In case of embedding, it uses the embedding on the graph and then use Fuzzy C-Means for extract the clusters.

    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
        n_clusters: int, default None. The number of clusters
        norm: boolean, default False. If is True the algorithm will compute the normalized laplacian matrix, otherwise the algorithm will compute the unnormalize laplacian matrix.
        embedding: string, default None. If None, no embeddings algorithms are applied. Outherwise, a supported embedding algorithm is applied on the eigenvectors of the laplacian matrix.
        d: int, default None, dimension for the embedding
        beta: float, default None, decay factor for HOPE embedding
        walk_len: int, default None, length of the random walks used in node2vec embedding.
        num_walks: int, default None, number of random walks each node computed in node2vec embedding.
    Returns:
        labels: extracted clusters
    """
    supported_embeddings = ["HOPE", "LaplacianEigenmaps", "Node2Vec"]

    if embedding is None:

        if norm:
            L = compute_normalized_laplacian(A)
        else:
            L = compute_laplacian_matrix(A)

        sortedEigenvalues, sortedEigenvectors = computeSortEigens(L)            #sorted eigenvalues/eigenvectors
        train = sortedEigenvectors[:, :n_clusters]

    else:

        if (embedding in supported_embeddings):

            G = from_numpy_array(A)
            if (embedding == "HOPE"): #embedding dimension (d) and decay factor (beta) as inputs
                model = HOPE(d=d, beta=beta)
            elif (embedding == "LaplacianEigenmaps"):
                model = LaplacianEigenmaps(d=d)
            elif (embedding == "Node2Vec"):
                model = Node2Vec(G, dimensions=d, walk_length=walk_len, num_walks=num_walks, workers = 1)

            if embedding == "Node2Vec":
                model = model.fit()
                train = model.wv.vectors
                train = np.array(train)
                if (train.shape[0] < A.shape[0]):
                    raise Exception("Parameters 'num_walks' and 'walk_length' for Node2Vec embedding are too small. Can't train the model.")
            else:
                train, t = model.learn_embedding(graph=G)
                train = np.array(train)

        else: raise Exception ("embedding {} not supported".format(embedding))

    fcm = FCM(n_clusters=n_clusters)
    fcm.fit(train)
    labels = fcm.predict(train)
    return labels

def ssc_shimalik(A, n_clusters = None):
    """
    Soft Spectral Clustering with Shi Malik Approach.
    Compute the eigenvalues and eigenvectors of the laplacian matrix solving the generalized eigenvalue problem.
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
        n_clusters: int, default None. The number of clusters.

    Returns:
        labels: extracted clusters
    """

    L = compute_laplacian_matrix(A)
    D = degree_matrix(A)
    eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
    idx = eigenvalues.argsort()
    sortedEigenvalues = eigenvalues[idx].real
    sortedEigenvectors = eigenvectors[:,idx].real

    if n_clusters is None:
        n_clusters = computeBestK(sortedEigenvalues)

    train = sortedEigenvectors[:, :n_clusters]

    fcm = FCM(n_clusters=n_clusters)
    fcm.fit(train)
    labels = fcm.predict(train)

    return labels

def hsc_shimalik(A, n_clusters = None):
    """
    Hard Spectral Clustering with Shi Malik Approach.
    Compute the eigenvalues and eigenvectors of the laplacian matrix solving the generalized eigenvalue problem.
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
        n_clusters: int, default None. The number of clusters
    Returns:
        labels: extracted clusters
    """

    L = compute_laplacian_matrix(A)
    D = degree_matrix(A)
    eigenvalues, eigenvectors = eigh(L, D, eigvals_only=False)
    idx = eigenvalues.argsort()
    sortedEigenvalues = eigenvalues[idx].real
    sortedEigenvectors = eigenvectors[:,idx].real

    if n_clusters is None:
        n_clusters = computeBestK(sortedEigenvalues)

    train = sortedEigenvectors[:, :n_clusters]

    km = KMeans(n_clusters=n_clusters).fit(train)
    labels = km.labels_

    return labels

def skl_spectral_clustering(A, n_clusters = None):
    """
    Spectral clustering following the scikit learn implementation.
    Parameters:
        A: numpy.array, the adjacency matrix of the graph.
        n_clusters: int, default None. The number of clusters
    Returns:
        labels: np.array, extracted clusters
    """
    clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed').fit(A)
    return clustering.labels_


#all spectral clustering choices
def unnorm_hsc(A, n_clusters = None, norm=False, embedding=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Unnormalized Hard Spectral Clustering.
    Parameters: see hardSpectralClustering
    Returns: see hardSpectralClustering
    """
    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels

def norm_hsc(A, n_clusters = None, norm=True, embedding=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Normalized Hard Spectral Clustering.
    Parameters: see hardSpectralClustering
    Returns: see hardSpectralClustering
    """
    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels

def unnorm_ssc(A, n_clusters = None, norm=False, embedding=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Unnormalized Soft Spectral Clustering.
    Parameters: see softSpectralClustering
    Returns: see softSpectralClustering
    """
    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels

def norm_ssc(A, n_clusters = None, norm=True, embedding=None, d=None, beta=None, walk_len=None, num_walks=None):
    """
    Normalized Soft Spectral Clustering.
    Parameters: see softSpectralClustering
    Returns: see softSpectralClustering
    """
    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels

#EMBEDDING + NORMAL CLUSTERING ALGORITHMS

def fuzzycmeans_hope(A, n_clusters = None, norm=False, embedding="HOPE", d=2, beta=0.01, walk_len=None, num_walks=None):
    """
    Soft Clustering + HOPE Embedding.
    Parameters: see softSpectralClustering
    Returns: see softSpectralClustering
    """
    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels

def kmeans_hope(A, n_clusters = None, norm=False, embedding="HOPE", d=2, beta=0.01, walk_len=None, num_walks=None):
    """
    Hard Clustering + HOPE Embedding.
    Parameters: see hardSpectralClustering
    Returns: see hardSpectralClustering
    """
    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels


def kmeans_node2vec(A, n_clusters = None, norm=False, embedding="Node2Vec", d=2, beta=None, walk_len=30, num_walks=30):
    """
    Hard Clustering + node2vec Embedding.
    Parameters: see hardSpectralClustering
    Returns: see hardSpectralClustering
    """
    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels


def fuzzycmeans_node2vec(A, n_clusters = None, norm=False, embedding="Node2Vec", d=2, beta=None, walk_len=30, num_walks=30):
    """
    Soft Clustering + node2vec Embedding.
    Parameters: see softSpectralClustering
    Returns: see softSpectralClustering
    """
    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels


def kmeans_laplacianeigenmaps(A, n_clusters = None, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None, walk_len=None, num_walks=None):
    """
    Hard Clustering + LaplacianEigenmaps Embedding.
    Parameters: see hardSpectralClustering
    Returns:  see hardSpectralClustering
    """
    labels = hardSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels


def fuzzycmeans_laplacianeigenmaps(A, n_clusters = None, norm=False, embedding="LaplacianEigenmaps", d=2, beta=None, walk_len=None, num_walks=None):
    """
    Soft Clustering + LaplacianEigenmaps Embedding.
    Parameters: see softSpectralClustering
    Returns:  see softSpectralClustering
    """
    labels = softSpectralClustering(A, n_clusters, norm, embedding, d, beta, walk_len, num_walks)
    return labels


#CENTRALITY MEASURES

def betweenness(G, residue_names_1, n=10):
    """
    Compute the betweenness centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the betweenness centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: P_coef[node]}, for each node is linked its betweenness centrality
    """
    bc = betweenness_centrality(G)
    #bc= betweenness_centrality_parallel(G)
    bc = {int (float (k)):v for k,v in bc.items()}
    dict_node_centrality = dict ()

    for i, cent in bc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_bc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by betweenness centrality".format(n))
    for d in sorted_bc[:n]:
        print(d)

    return dict_node_centrality

def eigenvector_c(G, residue_names_1, n=10):
    """
    Compute the eigenvector centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the eigenvector centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: eigenvector_centrality[node]}, for each node is linked its eigenvector centrality
    """
    ec = eigenvector_centrality(G, max_iter=10000)
    ec = {int (float (k)):v for k,v in ec.items()}
    dict_node_centrality = dict ()

    for i, cent in ec.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_ec = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by eigenvector centrality".format(n))
    for d in sorted_ec[:n]:
        print(d)

    return dict_node_centrality

def degree_c(G, residue_names_1, n=10):
    """
    Compute the degree centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the degree centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: degree_centrality[node]}, for each node is linked its degree centrality
    """
    dc = degree_centrality(G)
    dc = {int (float (k)):v for k,v in dc.items()}
    dict_node_centrality = dict ()

    for i, cent in dc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_dc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by degree centrality".format(n))
    for d in sorted_dc[:n]:
        print(d)

    return dict_node_centrality

def closeness(G, residue_names_1, n=10):
    """
    Compute the closeness centrality of the nodes of the graph using the networkx implementation.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the closeness centrality.
        residue_names_1: np.array, list of the residues names of the protein.
        n: int, default equals to 10, number of best nodes with highest centrality to print.
    Returns:
        centralities: dict {node: closeness_centrality[node]}, for each node is linked its closeness centrality
    """
    cc = closeness_centrality(G)
    cc = {int (float (k)):v for k,v in cc.items()}
    dict_node_centrality = dict ()

    for i, cent in cc.items():

        dict_node_centrality[residue_names_1[i]] = cent

    sorted_cc = sorted(dict_node_centrality.items(), key=itemgetter(1), reverse=True)
    print("Top {} nodes by closeness_centrality".format(n))
    for d in sorted_cc[:n]:
        print(d)

    return dict_node_centrality


def participation_coefs(G, labels, residue_names_1):
    """
    Compute the participation coefficient of the nodes of the graph given a partition.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the closeness centrality.
        labels: np.array, extracted clusters/communities
        residue_names_1: np.array, list of the residues names of the protein.
    Returns:
        P: dict {node: P_coef[node]}, for each node is linked its participation coefficient
    """
    A = to_numpy_array(G)
    n = A.shape[0]
    P = dict()
    k_s = np.zeros((n))

    for i in range(n):
        k_i = np.sum(A[i,:])
        k_si = 0

        for j in range(n):
            if (i!=j):
                if ((labels[i] == labels[j]) and (A[i,j]!=0)):#se il nodo i e il nodo j sono dello stesso cluster e c'è un arco che li connette
                    k_si += A[i,j]

        k_s[i] = k_si
        P[residue_names_1[i]] = 1 - (k_s[i]/k_i)**2

    return P

#END


#HERE WE FOUND NOT USED ALGORITHMS
def z_score(G, labels, residue_names_1):
    """
    Compute the z score of the nodes of the graph.
    Parameters:
        G: networkx.graph, the graph (the PCN) you want to compute the closeness centrality.
        labels: np.array, extracted clusters/communities
        residue_names_1: np.array, list of the residues names of the protein.
    Returns:
        z_s: dict {node: z_score[node]}, for each node is linked its z-score
    """
    A = to_numpy_array(G)
    n = A.shape[0]
    z_s = dict()
    k_s = np.zeros((n))

    for i in range(n):
        k_i = np.sum(A[i,:])
        k_si = 0

        for j in range(n):
            if (i!=j):
                if ((clusters[i] == clusters[j]) and (A[i,j]!=0)):
                    k_si += A[i,j]

        k_s[i] = k_si

        mean_k_si = np.mean(k_s)
        sigma_k_si = np.std(k_s)
        for i in range(n):
            k_i = np.sum(A[i,:])
            z_s[residue_names_1[i]] = (k_i - mean_k_si) / sigma_k_si

    return z_s

def plot_z_p(G, clusters):
    """
    Dentist Chair plot: x axis z-score and y axis Part_coef.
    """
    part_coefs = participation_coefs(G, clusters)
    z_scores = z_score(G, clusters)
    plt.scatter(np.array(list(part_coefs.values())).astype(float), np.array(list(z_scores.values())).astype(float), facecolors='none', edgecolors='b')
    plt.xlabel("part_coefs")
    plt.ylabel("z_scores")

def color_map_clustering(clusters):
    """
    Clustering color map plot
    """
    n = len(clusters)
    cluster_map = np.zeros((n,n), dtype = int)
    labels_unique, counts = np.unique(clusters, return_counts=True)

    for i, label in enumerate(clusters):
        for j, label2 in enumerate(clusters):
            if ((cluster_map[i][j] == 0) and (label==label2)):
                cluster_map[i][j] = label+1

    print(cluster_map)
    plt.figure(figsize=(8,8))
    plt.matshow(cluster_map, vmin=0, vmax=max(clusters)+1, cmap='inferno', fignum=1)
    plt.legend(labels=labels_unique)
    plt.show()

def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x


def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * 4
    node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.starmap(
        betweenness_centrality_subset,
        zip(
            [G] * num_chunks,
            node_chunks,
            [list(G)] * num_chunks,
            [True] * num_chunks,
            [None] * num_chunks,
        ),
    )

    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c

