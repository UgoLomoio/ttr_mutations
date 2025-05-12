#!/usr/bin/env python3
import requests
from pathlib import Path
import os 
from umap import UMAP
import torch 
import matplotlib.pyplot as plt 

from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay, f1_score, make_scorer
import seaborn as sns
import numpy as np
import plotly.express as px
import pandas as pd
from sklearn.preprocessing import label_binarize
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    roc_auc_score, classification_report
)
import plotly.io as pio
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
   
pio.kaleido.scope.mathjax = None

# Create a custom scorer for F1 score with zero_division handling
f1_scorer = make_scorer(f1_score, pos_label=0, zero_division=1)  # Set zero_division to 0 or 1 based on preference

def silentremove(filename):
    import os, errno
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred


def plot_roc_curve(y_true, y_probs, filename="roc_curve.pdf"):
    """
    Plots the ROC curve for a binary classification problem and saves it as a .pdf file.

    Parameters:
        y_true: array-like of shape (n_samples,) - True binary labels (0 or 1).
        y_probs: array-like of shape (n_samples, 2) - Predicted probabilities for each class.
        filename: str - The file name to save the ROC curve plot (default: "roc_curve.pdf").
    """

    # Compute ROC curve and AUC
    fpr, tpr, _ = roc_curve(y_true, y_probs[:, 0])  
    roc_auc = auc(fpr, tpr)

    # Initialize plot
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC Curve (AUC = {roc_auc:.2f})')

    # Plot the diagonal (chance line)
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', lw=2)

    # Set plot labels and title
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', fontsize=14)
    plt.ylabel('True Positive Rate', fontsize=14)
    plt.title('Receiver Operating Characteristic (ROC) Curve', fontsize=16)
    plt.legend(loc='lower right', fontsize=12)

    # Save and show plot
    silentremove(os.path.join(results_path, filename))
    plt.savefig(os.path.join(results_path, filename))
    plt.show()




def plot_roc_curve_multiclass(y_true, y_probs, n_classes, filename="roc_curve_multiclass.pdf"):
    """
    Plots the ROC curve for a multi-class problem and saves it as a .pdf file.
    
    Parameters:
        y_true: true labels (array-like)
        y_probs: predicted probabilities for each class (array-like, shape = [n_samples, n_classes])
        n_classes: number of classes
        filename: the file name to save the ROC curve plot (default: "roc_curve_multiclass.pdf")
    """
    # Binarize the true labels for multi-class ROC calculation
    y_true_bin = label_binarize(y_true, classes=range(n_classes))
    
    # Initialize plot
    plt.figure(figsize=(8, 6))
    
    # Compute ROC curve and ROC AUC for each class
    roc_auc = dict()
    for i in range(n_classes):
        if np.sum(y_true_bin[:, i]) == 0:  # No positive samples for this class
            print(f"Warning: No positive samples for class {i}, skipping...")
            continue

        fpr, tpr, _ = roc_curve(y_true_bin[:, i], y_probs[:, i])
        roc_auc[i] = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=2, label=f'Class {i} (AUC = {roc_auc[i]:.2f})')
    
    # Plot the diagonal line (chance line)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    
    # Set plot parameters
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve - Multi-Class')
    
    # Add a legend
    plt.legend(loc='lower right')
    
    # Save the ROC curve plot
    silentremove(results_path + sep + filename)
    plt.savefig(results_path + sep + filename)
    plt.show()

def plot_confusion_matrix(y_true, y_pred, filename="confusion_matrix.pdf"):
    """
    Plots the confusion matrix and saves it as a .pdf file.
    
    Parameters:
    - y_true: True labels
    - y_pred: Predicted labels
    - filename: The file name to save the confusion matrix plot (default: "confusion_matrix.pdf")
    """
    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred)
    
    # Check if the confusion matrix has all the classes
    unique_classes = np.unique(y_true)
    print(f"Unique classes in y_true: {unique_classes}")
    
    # Ensure labels are properly set
    disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=unique_classes)
    
    # Plot the confusion matrix
    fig, ax = plt.subplots(figsize=(8, 6))
    disp.plot(cmap=plt.cm.Blues, ax=ax, text_kw={"fontsize": 16})
    if len(unique_classes) == 3:
        labels = ['Pathogenic', 'Likely pathogenic', 'Benign']
    else:
        labels = ['Pathogenic', 'Benign']
    ax.xaxis.set_ticklabels(labels)
    ax.yaxis.set_ticklabels(labels)
    # Save the confusion matrix plot
    silentremove(results_path + sep + filename)
    plt.savefig(results_path + sep + filename)
    plt.show()

# Function to compute and print the F1 score
def print_f1_score(y_true, y_pred):
    """
    Computes and prints the F1 score.
    
    Parameters:
        y_true: true labels (array-like)
        y_pred: predicted labels (array-like)
    """
    f1 = f1_score(y_true, y_pred, pos_label=0, zero_division=0)
    print(f"F1 Score: {f1:.2f}")
    return f1

# Create a function to build the pipeline with each model
def build_pipeline(model):
    return Pipeline(
        steps=[
            #('pca', PCA(n_components=num_pca_components)),
            ('model', model)  # Pass the model here directly
        ]
    )


cwd = os.getcwd()
sep = os.sep

embeddings_path = cwd + sep + 'esm2_embeddings'
fasta_sequences_path = cwd + sep + 'fasta_sequences'
results_path = cwd + sep + "results"
if not os.path.exists(embeddings_path): 
    os.makedirs(embeddings_path)

df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["Mutation", "Consequence"]

df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely benign': 'Benign'})
#df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic'})    

knn_grid = [
    {
        'model': [KNeighborsClassifier()],
        'model__n_neighbors': [2, 5, 10, 20],
        'model__weights': ['uniform', 'distance'],
        'model__algorithm': ['ball_tree', 'kd_tree', 'brute'],
        'model__leaf_size' : [15, 30],
        'model__p' : [1, 2]
    }
    ]

svm_grid = [
    {
        'model': [SVC()],
        'model__C' : [0.1, 1.0, 2.0, 10.0],
        'model__kernel' : ['linear', 'poly', 'rbf', 'sigmoid'],
        'model__degree' : [3],
        'model__gamma': ['scale'],
        'model__class_weight': ['balanced']  # Add class_weight for imbalance
    }
]

rfr_grid = [
    {
        'model': [RandomForestClassifier()],
        'model__n_estimators' : [10, 20, 100, 200],
        'model__criterion': ['gini', 'entropy'],
        'model__max_features': ['sqrt', 'log2'],
        'model__min_samples_split' : [5, 10],
        'model__min_samples_leaf': [1, 4],
        'model__class_weight': ['balanced']  # Add class_weight for imbalance
    }
]


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
    
    key = os.getenv("NGC_API_KEY") or input("Paste the Run Key: ")

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

def print_all_metrics(y_true, y_pred, y_prob=None):
    """
    Computes and prints accuracy, F1 score, precision, recall, ROC AUC, and confusion matrix.
    
    Parameters:
    - y_true: True labels
    - y_pred: Predicted labels
    - y_prob: Predicted probabilities (needed for ROC AUC in binary or multi-class classification)
    """
    
    unique_classes = np.unique(y_true)
    n_classes = len(unique_classes)
    
    # Accuracy Score
    accuracy = accuracy_score(y_true, y_pred)
    print(f"Accuracy: {accuracy:.2f}")
    
    # Check if it's binary classification
    if n_classes == 2:
        # Binary classification: using pos_label=1
        f1 = f1_score(y_true, y_pred, zero_division=0, pos_label=1)
        precision = precision_score(y_true, y_pred, zero_division=0, pos_label=1)
        recall = recall_score(y_true, y_pred, pos_label=1)
        
        # For ROC AUC: if y_prob is 2D with probabilities for both classes, choose column 1
        if y_prob is not None:
            if y_prob.ndim == 2 and y_prob.shape[1] > 1:
                roc_auc = roc_auc_score(y_true, y_prob[:, 1])
            else:
                roc_auc = roc_auc_score(y_true, y_prob)
            print(f"ROC AUC: {roc_auc:.2f}")
        
        # Classification report with pos_label not needed here
        report = classification_report(y_true, y_pred, zero_division=0, target_names=[str(cls) for cls in unique_classes])
    else:
        # Multiclass classification: use macro averaging
        f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
        precision = precision_score(y_true, y_pred, average='macro', zero_division=0)
        recall = recall_score(y_true, y_pred, average='macro')
        
        # For ROC AUC in multiclass: using one-vs-rest (OvR) approach with macro averaging
        if y_prob is not None:
            try:
                roc_auc = roc_auc_score(y_true, y_prob, multi_class='ovr', average='macro')
                print(f"ROC AUC: {roc_auc:.2f}")
            except ValueError:
                print("ROC AUC Score: Not computable (check your probability estimates and y_true labels)")
        
        # Classification report for multiclass
        report = classification_report(y_true, y_pred, zero_division=0, target_names=[str(cls) for cls in unique_classes])
    
    print(f"F1 Score: {f1:.2f}")
    print(f"Precision: {precision:.2f}")
    print(f"Recall: {recall:.2f}")
    
    # Confusion Matrix
    cm = confusion_matrix(y_true, y_pred)
    print("Confusion Matrix:")
    print(cm)
    
    print("\nClassification Report:")
    print(report)

if __name__ == "__main__":

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
    print(embeddings.shape)
    emb2d = UMAP(n_components=2, n_neighbors=5, min_dist=0.0).fit_transform(embeddings)

        # Define the color mapping
    color_map = {
        "Pathogenic": "red",
        "Likely pathogenic": "orange",
        "Unknown": "black",
        "Benign": "green"
    }

    df = pd.DataFrame({"UMAP-1": emb2d[:,0], "UMAP-2": emb2d[:,1], "Consequence": labels, "Variant": variants})
    fig = px.scatter(df, x='UMAP-1', y='UMAP-2', color='Consequence', hover_data='Variant', color_discrete_map=color_map)

        # Set options common to all traces with fig.update_traces
    fig.update_traces(mode='markers', marker_line_width=2, marker_size=24)

    # Customize the hoverlabel
    fig.update_layout(hoverlabel=dict(bgcolor="white",
                                    font=dict(size=20, family="Arial", color="black"),
                                    align="left"))
    
    # Customize the legend
    fig.update_layout(legend = dict(font = dict(family = "Courier", size = 30, color = "black")),
                  legend_title = dict(text = "Clinical Significance", font = dict(family = "Courier", size = 32, color = "blue")))

    fig.update_layout(title=dict(text='UMAP projection of the ESM2 embeddings obtained for TTR variants'),
                    margin=dict(t=100) , # Increase top margin
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
                    width = 1920
    )
    #Save the figure as pdf with 600 DPI
    silentremove(results_path + sep + "umap_plot.png")
    pio.write_image(fig, results_path + sep + "umap_plot.pdf", format="pdf", scale=1)  # scale increases resolution

    fig.show()

    with open(fasta_sequences_path + sep + fasta_file) as file:
        lines = file.readlines()
        variant = lines[0][1:]
        seq = lines[1].replace(":", "").strip()

    start_pos = 0
    end_pos = len(seq)
    positions = np.arange(start_pos, end_pos)

    npz_filepath = embeddings_path + sep + "WT_embeddings.npz"
    logits = getembfromnpz(npz_filepath, "logits.npy")
    logits = logits.squeeze(axis = 0)
    logits = logits.transpose(0, 1)
    tokens = getembfromnpz(npz_filepath, "tokens.npy")

    logits = torch.from_numpy(logits)
    tokens = torch.from_numpy(tokens)

    probs = logits_to_probs(logits, tokens)

    aa_tokens = ['L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C']

    aa_indices = [i for i, token in enumerate(tokens) if token in aa_tokens]
    aa_indices = np.arange(4, 24, 1)
    probs = np.array(probs).T
    probs = probs[aa_indices, 1:-1]

    seq = list(seq)

    plt.figure(figsize=(16, 6))

    ax = sns.heatmap(
        probs,
        cmap='viridis',
        xticklabels=seq,
        yticklabels=aa_tokens,
        vmin=0, vmax=1,
        cbar_kws={
            'label': 'Probability', 
            'orientation': 'horizontal',
            'shrink': 0.5,  # Make the colorbar 50% of its default size
            'aspect': 25,   # Make it wider relative to height
            'pad': 0.1     # Reduce padding between plot and colorbar
        }
    )

    # Label the axes and the plot
    ax.set_xlabel("Position in Sequence")
    ax.set_ylabel("Token Labels")
    ax.set_title("Positional Token Probabilities")

    plt.xticks(rotation=0, ha='center')
    plt.tight_layout()
    silentremove(results_path + sep + "esm2_uncertainty.png")
    plt.savefig(results_path + sep + "esm2_uncertainty.pdf", dpi=600)
    #plt.show()

    if "Likely pathogenic" in labels:
        multiclass = True
        map_labels = {"Pathogenic": 0, "Benign": 1, "Likely pathogenic": 2, "Unknown": 3}
    else:
        multiclass = False
        map_labels = {"Pathogenic": 0, "Benign": 1, "Unknown": 3}
    map_labels_rev = {int_label: label for label, int_label in map_labels.items()}
    labels = [map_labels[label] for label in labels] 
    labels = np.array(labels)
    variants = np.array(variants)
    X = embeddings[labels != 3]
    unknown_idx = np.argwhere(labels == 3).flatten()
    y = labels[labels != 3]
    X_unknown = {variants[idx]: embeddings[idx] for idx in unknown_idx}
    variants = variants[labels != 3]

    train_size = 0.7
    X_train, X_test, y_train, y_test, variants_train, variants_test = train_test_split(X, y, variants, train_size=train_size, stratify = y)#, random_state=0)
    print(X_train.shape, y_train.shape, variants_train.shape)
    num_pca_components = 2
    pca = UMAP(n_components=2, n_neighbors=5, min_dist=0.0)
    X_train_pca = pca.fit_transform(X_train)

    df = pd.DataFrame(X_train_pca[:, :2], columns=['UMAP1', 'UMAP2'])
    df['Variant Effect'] = y_train
    df['Variant'] = variants_train

    # Create interactive scatter plot with Plotly
    fig = px.scatter(df, x='UMAP1', y='UMAP2', hover_data = 'Variant', color='Variant Effect', 
                    labels={'UMAP1': 'UMAP First Principal Component', 
                            'UMAP2': 'UMAP Second Principal Component'},
                    color_continuous_scale='Viridis')

            # Set options common to all traces with fig.update_traces
    fig.update_traces(mode='markers', marker_line_width=2, marker_size=24)

    # Customize the hoverlabel
    fig.update_layout(hoverlabel=dict(bgcolor="white",
                                    font=dict(size=20, family="Arial", color="black"),
                                    align="left"))
    
    # Customize the legend
    fig.update_layout(legend = dict(font = dict(family = "Courier", size = 30, color = "black")),
                  legend_title = dict(font = dict(family = "Courier", size = 50, color = "blue")))

    fig.update_layout(title=dict(text='UMAP projection of the ESM2 embeddings obtained for TTR variants'),
                    margin=dict(t=100) , # Increase top margin
                    yaxis_zeroline=False, xaxis_zeroline=False, 
                    xaxis_showgrid=False, yaxis_showgrid=False,
                    template='plotly_white',
                    #plot_bgcolor='rgba(0, 0, 0, 0)',
                    #paper_bgcolor='rgba(0, 0, 0, 0)',
                    title_font_size=50,
                    xaxis_title_font_size=40,
                    yaxis_title_font_size=40,
                    xaxis_tickfont_size=30,
                    yaxis_tickfont_size=30
    )
    # Show the plot
    fig.show()

    cls_list = [SVC, RandomForestClassifier, KNeighborsClassifier]
    param_grid_list = [knn_grid, svm_grid, rfr_grid]

    # Initialize lists to store results
    result_list = []
    grid_list = []

    # Loop through each classifier and its corresponding parameter grid
    for cls, param_grid in zip(cls_list, param_grid_list):
        print(f"Running GridSearchCV for {cls.__class__.__name__}")

        # Build pipeline with the current model
        pipe = build_pipeline(cls)
        
        # Initialize GridSearchCV with the pipeline
        grid = GridSearchCV(
            estimator=pipe,
            param_grid=param_grid,
            scoring='accuracy',#f1_scorer,
            verbose=1,
            n_jobs=-1  # Use all available cores
        )
        
        # Fit the grid search on the training data
        grid.fit(X_train, y_train)
        
        # Store the results and the grid search object
        result_list.append(pd.DataFrame.from_dict(grid.cv_results_))
        grid_list.append(grid)

    #for res in result_list:
    #    print(res.sort_values('rank_test_score').head(5))


    best_model = grid.best_estimator_
    print(best_model)
    print("Inference on unknown clinical significance mutations.")
    predictions_unknown = {}
    for mutation, x in X_unknown.items():
        x = x.reshape(1, -1)
        prediction = best_model.predict(x)[0]
        predictions_unknown[mutation] = map_labels_rev[prediction]
        print(mutation, map_labels_rev[prediction])

    with open(results_path + sep + "predictions_unknown.txt", "w") as file:
        file.write(str(predictions_unknown))

    print("Results on test set")
    # Make predictions
    y_pred = best_model.predict(X_test)
    y_probs = best_model.predict_proba(X_test) # Probabilities for ROC curve

    # Plot and save the ROC curve
    if multiclass:
        plot_roc_curve_multiclass(y_test, y_probs, n_classes = 3, filename="roc_curve.pdf")
    else:
        plot_roc_curve(y_test, y_probs, filename="roc_curve.pdf")
    # Plot and save the confusion matrix
    plot_confusion_matrix(y_test, y_pred, filename="confusion_matrix.pdf")
    print_all_metrics(y_test, y_pred, y_prob=None)

    # Make predictions
    y_pred = best_model.predict(X)
    predictions_all = {}
    for i, mutation in enumerate(variants): 
        x = X[i]
        x = x.reshape(1, -1)
        prediction = best_model.predict(x)[0]
        predictions_all[mutation] = map_labels_rev[prediction]
        #print(mutation, map_labels_rev[prediction])

    with open(results_path + sep + "predictions.txt", "w") as file:
        file.write(str(predictions_all))
