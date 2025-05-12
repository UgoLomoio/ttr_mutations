import pandas as pd
import os 

cwd = os.getcwd()
sep = os.sep 
results_path = cwd + sep + "mutant_stability"
file = "Analysis-Mutations3A4D_cleaned_ChainA"
filename = results_path + sep + f"{file}.xlsx"
outfilename = results_path + sep + f"{file}_latex.txt"

def df2latex(filename, outfilename):

    if filename.split(".")[-1] == "xlsx":
        df = pd.read_excel(filename, index_col=0)
    elif filename.split(".")[-1] == "csv":
        df = pd.read_csv(filename, index_col=0)
    else: return None

    df.to_latex(outfilename)

if __name__ == "__main__":

    df2latex(filename, outfilename)