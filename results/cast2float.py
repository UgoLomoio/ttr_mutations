import pandas as pd
import os

def convert_columns_to_float(directory, float_columns):
    try:
        csv_files = [f for f in os.listdir(directory) if f.endswith('.csv') and "converted" not in f]
        
        for file in csv_files:
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path, sep=',')
            
            for col in float_columns:
                for col_df in df.columns:
                    if col in col_df:
                        df[col] = df[col].astype(str).str.replace('.', ',')#.astype(float)
                else:
                    print(f"Warning: Column '{col}' not found in {file}.")
            
            output_file = os.path.join(directory, file.replace('.csv', '_converted.csv'))
            df.to_csv(output_file, index=False, sep=';')
            #os.remove(file_path)
            print(f"Successfully converted columns and saved to {output_file}. Deleted original file: {file}")
    except Exception as e:
        print(f"Error processing files: {e}")

if __name__ == "__main__":
    
    cwd = os.getcwd()
    sep = os.sep 
    directory = cwd + sep + "results" + sep # Specify the directory containing CSV files
    float_columns = ["ΔΔG prediction"] # Specify column names to convert
    convert_columns_to_float(directory, float_columns)