import pandas as pd 
import ast 
import os 
import numpy as np 
from sklearn.metrics import accuracy_score

cwd = os.getcwd()
sep = os.sep

results_path = cwd + sep + "results"

methods = ["Hierarchical clustering Average Linkage", "DBSCAN clustering"]
complex_types = ["monomer", "tetramer"]

if __name__ == "__main__":

    for method in methods:
        
        print(method)
        df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
        df_consequences.columns = ["Mutation", "Consequence"]

        # Remove unwanted rows
        df_consequences = df_consequences[~df_consequences['Consequence'].isin(['Not found', 'Variant of uncertain significance', 'Unknown'])]

        # Replace values
        df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic', 'Likely benign': 'Benign'})
        
        data_to_append = {"Mutation":"WT", "Consequence":"Benign"}
        df_consequences_append = pd.DataFrame(data_to_append, columns = ["Mutation", "Consequence"], index=[0])
        df_consequences = pd.concat([df_consequences, df_consequences_append], ignore_index=True)
        df_consequences.set_index("Mutation", inplace=True)
        consequences_dict = df_consequences.to_dict(orient="index")
        consequences_dict = {mutation: consequence["Consequence"] for mutation, consequence in consequences_dict.items()}

        for complex_type in complex_types:
            filepath = results_path + sep + "clusters_{}_{}.txt".format(method, complex_type)
            content = open(filepath, 'r').read()
            clusters_dict = ast.literal_eval(content)
            counts = {}
            temp = {}
            for cluster_idx, cluster_labels in clusters_dict.items():
                if "wt" in cluster_labels:
                    cluster_idx = "Benign"
                else:
                    counts[cluster_idx] = len(cluster_labels)
                temp[cluster_idx] = cluster_labels
            
            clusters_dict = temp
            if len(list(clusters_dict.keys())) == 1:
                continue
            elif len(list(clusters_dict.keys())) > 2:
                counts_v = np.array(list(counts.values()))
                where_min = np.argmin(counts_v)      
                cluster_idx_delete = str(list(counts.keys())[where_min])
                del clusters_dict[cluster_idx_delete]
            
            path_cluster_idx = [cluster_idx for cluster_idx in clusters_dict.keys() if cluster_idx != "Benign"][0]
            clusters_dict["Pathogenic"] = clusters_dict[path_cluster_idx]
            del clusters_dict[path_cluster_idx]
            #print(clusters_dict)
            map_mutant_cluster_dict = {
                cluster_label.upper(): cluster_idx 
                for cluster_idx, cluster_labels in clusters_dict.items() 
                for cluster_label in cluster_labels
            }
            temp = {}
            for cluster_label, cluster_idx in map_mutant_cluster_dict.items():
                if cluster_label != "WT":
                    w = cluster_label[0]
                    m = cluster_label[-1]
                    resi = int(cluster_label[1:-1])
                    new_resi = resi - 20
                    new_cluster_label = w + str(new_resi) + m
                    temp[new_cluster_label] = cluster_idx
                else:
                    temp[cluster_label] = cluster_idx

            map_mutant_cluster_dict = temp

            map_consequence = {"Benign": 0, "Pathogenic": 1}
            predictions = []
            consequences = []
            for mutation in df_consequences.index:
                if mutation in map_mutant_cluster_dict:
                    predictions.append(map_consequence[map_mutant_cluster_dict[mutation]])  
                    consequences.append(map_consequence[df_consequences["Consequence"][mutation]])
                else:
                    print("Can't predict mutant {}".format(mutation))
            corr = np.corrcoef(consequences, predictions)[0, 1]
            print("{} correlation between predictions and real consequences: {}".format(complex_type, corr))
            acc = accuracy_score(consequences, predictions)
            print("{} accuracy {}".format(complex_type, acc))