import requests
import os 
import time 
import pandas as pd 

cwd = os.getcwd()
sep = os.sep 

pdb_path = cwd + sep + "pdbs-alphafold"
results_path = cwd + sep + "results"
#pdb_path = cwd + sep + "pdbs-real"

mutations_list_dict = {
    "monomer": cwd + sep + "ttrmutations_try.txt",
    "dimer": cwd + sep + "ttrmutations_dimer_try.txt",
    "tetramer": cwd + sep + "ttrmutations_tetramer_try.txt"
}
pdb_file_dict = {
    "monomer": pdb_path + sep + "monomer" + sep + "wt.pdb",
    "dimer": pdb_path + sep + "dimer" + sep + "wt-dimer.pdb",
    "tetramer": pdb_path + sep + "tetramer" + sep + "wt-tetramer.pdb",
}

url_single = "https://biosig.lab.uq.edu.au/dynamut2/api/prediction_list"
url_multiple = "https://biosig.lab.uq.edu.au/dynamut2/api/prediction_mm"
headers = {"Content-Type": "application/json"}

complex_types = list(pdb_file_dict.keys())

def convert_columns_to_float(directory, float_columns):
    try:
        csv_files = [f for f in os.listdir(directory) if f.endswith('.csv') and "converted" not in f]
        
        for file in csv_files:
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path, sep=',')
            
            for col in float_columns:
                if col in df.columns:
                    df[col] = df[col].astype(str).str.replace('.', ',')#.astype(float)
                else:
                    print(f"Warning: Column '{col}' not found in {file}.")
            
            output_file = os.path.join(directory, file.replace('.csv', '_converted.csv'))
            df.to_csv(output_file, index=False, sep=';')
            #os.remove(file_path)
            print(f"Successfully converted columns and saved to {output_file}. Deleted original file: {file}")
    except Exception as e:
        print(f"Error processing files: {e}")

# Function to read a list from a file
def read_list_from_file(file_path):
    with open(file_path, 'r') as file:
        # Read lines and remove any leading/trailing whitespace characters (like newline)
        my_list = [line.strip().replace("A ", "") for line in file.readlines()]
    return my_list

if __name__ == "__main__":

    datas = {}
    for complex_type in complex_types:
        print(complex_type)
        pdb_file = pdb_file_dict[complex_type]
        mutations_list = mutations_list_dict[complex_type]
        files_to_submit = {"pdb_file": open(pdb_file,"rb"),"mutations_list": open(mutations_list,"rb")}
        if complex_type == "monomer":
            req = requests.post(url_single, files=files_to_submit, data={})
        else:
            req = requests.post(url_multiple, files=files_to_submit) 
       
        job_id = req.json()["job_id"]
        params = {
            "job_id": job_id,
        }
        datas[complex_type] = params

    results = {}
    complex_todo = complex_types.copy()
    mutations = read_list_from_file(cwd + sep + "ttrmutations.txt")
    df_cols = [mutations]
    columns = ["Mutation"]
    while len(complex_todo) > 0:
        for complex_type in complex_todo:
            print("Checking job status for: ", complex_type)
            params = datas[complex_type]
            if complex_type == "monomer":
                req = requests.get(url_single, data=params)
            else:
                req = requests.get(url_multiple, data=params)
            response = req.json()
            #print(response)

            if "status" in response.keys():
                if response["status"] == "ERROR":
                    job_id = params["job_id"]
                    raise Exception("DynaMUT2 APIs raised an error processing job {} related to {}. Please check PDB file and mutation .txt file.".format(job_id, complex_type))
                else:
                    continue
            else:
                print("Completed job: ", complex_type)
                complex_todo.remove(complex_type)
                del datas[complex_type]
                del response["results_page"]
                col = []
                for i, result_dict in response.items():
                    del result_dict["chain"]
                    mutation = result_dict["mutation"]
                    prediction = result_dict["prediction"]
                    column_name = "ΔΔG prediction " + complex_type 
                    col.append(prediction)
                df_cols.append(col)
                columns.append(column_name)
            time.sleep(5)
        
    df = pd.DataFrame(df_cols, columns = columns)
    df.to_csv(results_path + sep + "results_dynamut2.csv")