#starting with the ESM2 embeddings we do UMAP projections and then classify as benign the mutations that are near the wild type projection 

import requests
from pathlib import Path
import os 
from umap import UMAP
import torch 
import pandas as pd
import numpy as np
from scipy.spatial import distance
from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score

cwd = os.getcwd()
sep = os.sep

parent_dir = os.path.dirname(cwd)
embeddings_path = parent_dir + sep + 'esm2_embeddings'
fasta_sequences_path = parent_dir + sep + "fasta_sequences"
results_path = parent_dir + sep + "results"
if not os.path.exists(embeddings_path): 
    os.makedirs(embeddings_path)

df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["Mutation", "Consequence"]
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely benign': 'Benign'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic'})

if not os.path.exists(embeddings_path): 
    os.makedirs(embeddings_path)

def getembfromnpz(npz_file, filename="embeddings.npy"):
    """
    Extracts and returns the contents of 'embeddings.npy' from a given .npz file.
    
    Args:
        npz_file (str): Path to the .npz file.
        filename (str, optional): Name of the file inside the .npz to extract. Default is "embeddings.npy".
    
    Returns:
        numpy.ndarray: The extracted embeddings as a NumPy array.
    """
    with np.load(npz_file) as data:
        if filename in data:
            return data[filename]
        else:
            raise KeyError(f"'{filename}' not found in the NPZ file.")


def get_embeddings(fasta_file):
    
    key = "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"#os.getenv("NGC_API_KEY") or input("Paste the Run Key: ")

    with open(fasta_sequences_path + sep + fasta_file) as file:
        lines = file.readlines()
        variant = lines[0][1:]
        SEQ = lines[1].replace(":", "").strip()
    print("SEQUENCE {}: {}".format(variant, SEQ))
    EMB_FORMAT = "npz"

    response = requests.post(
        url="https://health.api.nvidia.com/v1/biology/meta/esm2-650m",
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {key}"
        },
        json={
            "sequences": [SEQ],
            "format": EMB_FORMAT,
        },
    )

    if response.status_code == 200:
        print("SUCCESS")
        ext = "zip" if response.headers["Content-Type"] == "application/zip" else EMB_FORMAT
        with open(Path(f"{embeddings_path}{sep}{fasta_file.split('.')[0]}_embeddings.{ext}"), "wb") as fb:
            fb.write(response.content)
    else:
        print("FAILED")



def logits_to_probs(
        logits: torch.Tensor, tokens: dict
) -> torch.Tensor:
    """Convert token logits to probabilities

    Args:
        logits (torch.Tensor): logits tensor with the [batch, sequence, hidden] dimensions
        tokens (torch.Tensor): ESM2 tokens

    Returns:
        probabilities (torch.Tensor): probability tensor with [batch, sequence, tokenizer.vocab_size]
    """
    aa_tokens = ['L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C']
    extra_indices = [i for i, token in enumerate(tokens) if token not in aa_tokens]

    aa_logits = logits[..., :33]  # filter out the 95 paddings and only keep 33 vocab positions
    aa_logits[..., extra_indices] = - torch.inf  # force non-amino acid token probs to zero
    return torch.softmax(aa_logits, dim=-1)


def predict(benign_variants, other_variants, threshold=0.1):

    y_pred = [] 
    y_true = []
    map_classes = {"Benign": 0, "Pathogenic": 1}

    for variant in other_variants:
        x_variant = df[df["Variant"] == variant][["UMAP-1", "UMAP-2"]].values.flatten()
        distances = []
        for variant_b in benign_variants:
            x_benign = df[df["Variant"] == variant_b][["UMAP-1", "UMAP-2"]].values.flatten()
            dist = distance.euclidean(x_variant, x_benign)
            #print("Distance between {} and {}: {:.4f}".format(variant, variant_b, dist))
            distances.append(dist)

        if len(distances) > 0:
            # Find the minimum distance to any benign variant
            dist = min(distances)

        if dist < threshold:
            prediction = 0 # Benign
        else: 
            prediction = 1 # Pathogenic

        y_pred.append(prediction)
        true = df_consequences[df_consequences["Mutation"] == variant]["Consequence"].values[0]
        true = map_classes[true]
        #print(f"Variant: {variant}, Prediction: {prediction}, True: {true}")
        y_true.append(true)
        
    return distances, y_pred, y_true

if __name__ == "__main__":
    
    embeddings_umap_path = cwd + sep + "embeddings_umap.csv"
    if os.path.exists(embeddings_umap_path):
        print("UMAP embeddings already exist, reading...")
        df = pd.read_csv(embeddings_umap_path)

    else:
        print("Generating UMAP embeddings...")
        fastas = os.listdir(fasta_sequences_path)
        npz_files = os.listdir(embeddings_path)
        
        model_name = 'esm2_t33_650M_UR50D'

        for fasta_file in fastas:
            
            variant = fasta_file.split(".")[0]
            npz_file = variant + "_" + "embeddings.npz"
            if npz_file in npz_files:
                print("Embeddings for {} already exist".format(variant))
            else:
                get_embeddings(fasta_file)

        labels = []
        variants = []
        embeddings = []
        for npz_file in npz_files:
            if "WT" in npz_file.upper():
                label = "Benign"
                variant = "wt"
            else:
                #print(df_consequences)
                variant = npz_file.split(".")[0].split("_")[1].upper()
                label = df_consequences[df_consequences["Mutation"] == variant]["Consequence"].values[0]

            variants.append(variant)
            labels.append(label)

            npz_filepath = embeddings_path + sep + npz_file
            embedding = getembfromnpz(npz_filepath)
            embedding = embedding.flatten()
            embeddings.append(embedding)


        embeddings = np.row_stack(embeddings)
        emb2d = UMAP(n_components=2, n_neighbors=3, min_dist=0.0, metric = "euclidean", random_state = 42).fit_transform(embeddings)
        df = pd.DataFrame({"UMAP-1": emb2d[:,0], "UMAP-2": emb2d[:,1], "Consequence": labels, "Variant": variants})
        #df.to_csv(cwd + sep + "embeddings_umap.csv", index=False)

    benign_variants = ["wt"]
    #print("Benign variants: ", benign_variants)

    df_consequences_unknown = df_consequences[df_consequences["Consequence"] == "Unknown"]
    df_consequences = df_consequences[df_consequences["Consequence"] != "Unknown"]
    variants = df_consequences["Mutation"].values
    
    thresholds = np.arange(0.05, 1, 0.025)  # Adjust the range and step size as needed
    best_rocauc = 0.0
    best_threshold = 0.0
    best_preds = None
    
    for threshold in thresholds:
        distances, y_pred, y_true = predict(benign_variants, variants, threshold=threshold)
        if len(np.unique(y_pred)) < 2:
            print(f"Skipping threshold {threshold:.2f} due to insufficient classes in y_pred.")
            continue
        if len(np.unique(y_true)) < 2:
            print(f"Skipping threshold {threshold:.2f} due to insufficient classes in y_true.")
            continue
        rocauc = roc_auc_score(y_true, y_pred)
        if rocauc > best_rocauc:
            print(f"New best threshold found: {threshold:.2f} with ROCAUC: {rocauc:.4f}")
            best_preds = y_pred
            best_rocauc = rocauc
            best_threshold = threshold
        print(f"Threshold: {threshold:.2f}, ROCAUC: {rocauc:.4f}")
    
    print(f"Best threshold: {best_threshold:.2f}, Best ROCAUC: {best_rocauc:.4f}")
    
    print(y_true)
    
    cf_matrix = confusion_matrix(y_true, best_preds, labels=[0, 1])
    print("Confusion Matrix:")
    print(cf_matrix)  

    print("Classification Report:")
    print(classification_report(y_true, best_preds, labels=[0, 1], target_names=["Benign", "Pathogenic"]))

    #Save the best predictions
    
    #For each variant extract the distance from WT 
    distances_wt = []
    for variant in variants:   
        x_variant = df[df["Variant"] == variant][["UMAP-1", "UMAP-2"]].values.flatten()
        x_wt = df[df["Variant"] == "wt"][["UMAP-1", "UMAP-2"]].values.flatten()
        dist = distance.euclidean(x_variant, x_wt)
        distances_wt.append(dist)
    distances = np.array(distances_wt)
    df_pred = pd.DataFrame({"Mutation": variants, "Pred_class": best_preds, "Distance_wt": distances})
    df_pred.to_csv(cwd + sep + "esm2_umap_predictions.csv", index=False)



    print("Unknown Consequences Predictions:")
    for variant in df_consequences_unknown["Mutation"].values:
        x_variant = df[df["Variant"] == variant][["UMAP-1", "UMAP-2"]].values.flatten()
        distances = []
        for variant_b in benign_variants:
            x_benign = df[df["Variant"] == variant_b][["UMAP-1", "UMAP-2"]].values.flatten()
            dist = distance.euclidean(x_variant, x_benign)
            distances.append(dist)

        if len(distances) > 0:
            # Find the minimum distance to any benign variant
            dist = min(distances)

        if dist < best_threshold:
            prediction = "Benign"
        else: 
            prediction = "Pathogenic"

        print(f"Variant: {variant}, Prediction: {prediction}")
