import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import os 
import pandas as pd 
import numpy as np 
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt 

from sklearn.metrics import roc_curve, auc, confusion_matrix, ConfusionMatrixDisplay, f1_score
import pandas as pd
from sklearn.preprocessing import label_binarize
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    roc_auc_score, classification_report
)

cwd = os.getcwd()
sep = os.sep 
embeddings_path = cwd + sep + 'esm2_embeddings'
results_path = cwd + sep + "results"


df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["Mutation", "Consequence"]

df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely benign': 'Benign'})
#df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic'})   

# Define the MLP model
class MLP(nn.Module):
    def __init__(self, input_size=1280, hidden_size=256, output_size=2):
        super(MLP, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, output_size),
            nn.Softmax(dim=1)
        )
    
    def forward(self, x):
        return self.model(x)

    def predict(self, x):
        # Ensure that the model is in evaluation mode
        self.eval()
        with torch.no_grad():  # No need to track gradients for predictions
            outputs = self(x)  # Get the model's raw outputs
            _, predicted_classes = torch.max(outputs, 1)  # Choose the class with the highest probability
        return predicted_classes.long()

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

# Generate synthetic dataset
def generate_data(samples=1000, input_size=1280):
    X = torch.randn(samples, input_size)  # Random input tensors
    y = torch.randint(0, 2, (samples,))  # Binary labels (0 or 1)
    return X, y

# Training loop
def train_model(model, train_loader, val_loader, epochs=10):
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        for X_batch, y_batch in train_loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            optimizer.zero_grad()
            outputs = model(X_batch)
            loss = criterion(outputs, y_batch)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        # Validation
        model.eval()
        correct, total = 0, 0
        with torch.no_grad():
            for X_batch, y_batch in val_loader:
                X_batch, y_batch = X_batch.to(device), y_batch.to(device)
                outputs = model(X_batch)
                _, predicted = torch.max(outputs, 1)
                total += y_batch.size(0)
                correct += (predicted == y_batch).sum().item()
        
        print(f"Epoch {epoch+1}/{epochs}, Loss: {total_loss/len(train_loader):.4f}, Accuracy: {correct/total:.4f}")


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


def plot_roc_curve(y_true, y_probs, filename="roc_curve.png"):
    """
    Plots the ROC curve for a binary classification problem and saves it as a .png file.

    Parameters:
        y_true: array-like of shape (n_samples,) - True binary labels (0 or 1).
        y_probs: array-like of shape (n_samples, 2) - Predicted probabilities for each class.
        filename: str - The file name to save the ROC curve plot (default: "roc_curve.png").
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
    plt.savefig(os.path.join(results_path, filename))
    plt.show()


def plot_roc_curve_multiclass(y_true, y_probs, n_classes, filename="roc_curve_multiclass.png"):
    """
    Plots the ROC curve for a multi-class problem and saves it as a .png file.
    
    Parameters:
        y_true: true labels (array-like)
        y_probs: predicted probabilities for each class (array-like, shape = [n_samples, n_classes])
        n_classes: number of classes
        filename: the file name to save the ROC curve plot (default: "roc_curve_multiclass.png")
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
    plt.savefig(results_path + sep + filename)
    plt.show()

def plot_confusion_matrix(y_true, y_pred, filename="confusion_matrix.png"):
    """
    Plots the confusion matrix and saves it as a .png file.
    
    Parameters:
    - y_true: True labels
    - y_pred: Predicted labels
    - filename: The file name to save the confusion matrix plot (default: "confusion_matrix.png")
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
if __name__ == "__main__":

    npz_files = os.listdir(embeddings_path)
        
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

    device = "cpu"#torch.device("cuda" if torch.cuda.is_available() else "cpu")
    train_size = 0.7
    X_train, X_test, y_train, y_test, variants_train, variants_test = train_test_split(X, y, variants, train_size=train_size, stratify = y)#, random_state=0)
    X_train = torch.from_numpy(X_train).to(device).float()
    X_test = torch.from_numpy(X_test).to(device).float()
    y_train = torch.from_numpy(y_train).to(device).long() 
    y_test = torch.from_numpy(y_test).to(device).long()
    dataset_train = TensorDataset(X_train, y_train)
    dataset_test = TensorDataset(X_test, y_test)
    train_loader = DataLoader(dataset_train, batch_size=16, shuffle=True)
    test_loader = DataLoader(dataset_test, batch_size=1)

    # Initialize model, loss function, and optimizer
    if multiclass:
        n_classes = 3
    else:
        n_classes = 2
    model = MLP(output_size=n_classes).to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=0.0005)
    train_model(model, train_loader, test_loader, epochs = 400)

    y_preds = []
    y_probs = []

    # Iterate over the test set (if using DataLoader)
    for inputs, _ in test_loader:  # Here, we assume labels are in the second part of the batch
        y_pred = model.predict(inputs).cpu().detach().numpy()
        y_preds.extend(y_pred)  # Add predictions to the list
        y_prob = model.forward(inputs).cpu().detach().numpy()
        y_probs.extend(y_prob)  # Add probs to the list
    y_probs = np.array(y_probs)
    print(y_probs.shape)
    # Plot and save the ROC curve
    if multiclass:
        plot_roc_curve_multiclass(y_test, y_probs, n_classes = 3, filename="roc_curve.png")
    else:
        plot_roc_curve(y_test, y_probs, filename="roc_curve.png")
    # Plot and save the confusion matrix
    plot_confusion_matrix(y_test, y_preds, filename="confusion_matrix.png")
    print_all_metrics(y_test, y_preds, y_prob=None)
