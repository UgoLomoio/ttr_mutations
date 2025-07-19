import os 
import pandas as pd
from sklearn.metrics import roc_auc_score, confusion_matrix, classification_report
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
import numpy as np

def use_hgcn(df):
    df_copy = df.copy()
    for mutation in df['variant']:
        aa_original = mutation[0]
        aa_mutated = mutation[-1]
        aa_index = mutation[1:-1]
        aa_index_hugo = int(aa_index) + 20
        mutation_hugo = f"{aa_original}{aa_index_hugo}{aa_mutated}"
        #print(f"Processing mutation: {mutation}, HUGO mutation: {mutation_hugo}")
        df_copy.loc[df['variant'] == mutation, 'variant'] = mutation_hugo
    return df_copy


def predict_threshold(df, threshold=0.5):

    y_preds = []
    for i, row in df.iterrows():
        y_pred = 1 if row['score'] >= threshold else 0
        y_preds.append(y_pred)
    return np.array(y_preds)

 
cwd = os.getcwd()
sep = os.sep

parent_dir = os.path.dirname(cwd)

results_path = parent_dir + sep + "results"
df_consequences = pd.read_csv(results_path + sep + "mutations_consequences.csv")
df_consequences.columns = ["variant", "Consequence"]
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Not found': 'Unknown', 'Variant of uncertain significance': 'Unknown'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely benign': 'Benign'})
df_consequences['Consequence'] = df_consequences['Consequence'].replace({'Likely pathogenic': 'Pathogenic'})
df_consequences = use_hgcn(df_consequences)

predictions_vesm_path = cwd + sep + "vesmvariant_score.csv"
predictions_df = pd.read_csv(predictions_vesm_path)
predictions_df = predictions_df[predictions_df['variant'].isin(df_consequences['variant'].values)] # Filter to keep only variants present in consequences

variant_predictions = predictions_df['variant'].values
variant_target = df_consequences['variant'].values

if __name__ == "__main__":

    map_classes = {
        'Pathogenic': 1,
        'Benign': 0
    }
    map_classes_reverse = {v: k for k, v in map_classes.items()}

    """
    dt = DecisionTreeClassifier(random_state=42, 
                                max_depth=2 , 
                                min_samples_split=2, 
                                min_samples_leaf=1, 
                                class_weight='balanced',
                                criterion='gini')
    """

    df_consequences_unknown = df_consequences[df_consequences['Consequence'] == 'Unknown']
    df_consequences = df_consequences[df_consequences['Consequence'] != 'Unknown']
    predictions_df_unknown = predictions_df[predictions_df['variant'].isin(df_consequences_unknown['variant'].values)]
    predictions_df = predictions_df[predictions_df['variant'].isin(df_consequences['variant'].values)]
    
    X = predictions_df.drop(columns=['variant']).values
    y = df_consequences['Consequence'].map(map_classes)
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=y, random_state=42)
    #dt.fit(X_train, y_train)
    #y_pred = dt.predict(X)
    
    thresholds = np.arange(0.0, 1.0, 0.01)
    best_auc = 0
    best_threshold = 0.0
    best_pred = None
    for threshold in thresholds:
        
        y_pred = predict_threshold(predictions_df, threshold=threshold)
        predictions_df['Classification'] = y_pred

        y_pred = []
        y_true = []
        y_pred_unknown = {}
        for variant in df_consequences['variant']:
            if variant not in predictions_df['variant'].values:
                print(f"Variant {variant} not found in VESM predictions.")
                continue
            classification = predictions_df.loc[predictions_df['variant'] == variant, 'Classification'].values[0]
            classification = map_classes_reverse[classification]
            consequence = df_consequences.loc[df_consequences['variant'] == variant, 'Consequence'].values[0]
            y_true.append(map_classes[consequence])
            y_pred.append(map_classes[classification])
            #print(f"Variant: {variant}, Classification: {classification}, Consequence: {consequence}")  

        y_true = np.array(y_true)
        y_pred = np.array(y_pred)
        auc = roc_auc_score(y_true, y_pred)
        if auc > best_auc:
            print(f"New best AUC: {auc:.4f} at threshold {threshold:.2f}")
            best_auc = auc
            best_threshold = threshold
            best_pred = y_pred
        
    auc = roc_auc_score(y_true, best_pred)
    print(f"ROC-AUC: {auc:.4f}")
    cf_matrix = confusion_matrix(y_true, best_pred)
    print("Confusion Matrix:")
    print(cf_matrix)
    print("Classification Report:")
    print(classification_report(y_true, best_pred, target_names=list(map_classes.keys())[::-1]))

    X_unknown = predictions_df_unknown.drop(columns=['variant']).values
    #predictions_df_unknown['Classification'] = dt.predict(X_unknown)
    predictions_df_unknown['Classification'] = predict_threshold(predictions_df_unknown, threshold=best_threshold)
    for variant in df_consequences_unknown['variant']:
        if variant not in predictions_df_unknown['variant'].values:
            print(f"Variant {variant} not found in VESM predictions.")
            continue
        classification = predictions_df_unknown.loc[predictions_df_unknown['variant'] == variant, 'Classification'].values[0]
        classification = map_classes_reverse[classification]
        y_pred_unknown[variant] = classification

    #save the predictions
    df_to_save = predictions_df[['variant', 'Classification', 'score']].copy()
    df_to_save = df_to_save.rename(columns={'variant': 'Mutation', 'Classification': 'Pred_class', 'score': 'Pathogenicity_probability'})
    df_to_save['Pred_class'] = best_pred
    df_to_save['Pred_class'] = [map_classes_reverse[pred_class] for pred_class in df_to_save['Pred_class'].values]
    df_to_save.to_csv(cwd + sep + "vesm_predictions.csv", index=False) 

    """
    print("Variants with unknown consequences:")
    for variant, classification in y_pred_unknown.items():
        print(f"Variant: {variant}, Classification: {classification}")
    """