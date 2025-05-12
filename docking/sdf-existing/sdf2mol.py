import os 

cwd = os.getcwd()
sep = os.sep 


parent_dir = os.path.dirname(cwd)
docking_dir = parent_dir + sep + "docking"
sdf_path = docking_dir + sep + "sdf-existing"
files = os.listdir(sdf_path)
print(f"Files in {sdf_path}: {files}")
files = [file for file in files if file.split(".")[-1] == "sdf"]

if __name__ == "__main__":
    for file in files:
        file_path = os.path.join(sdf_path, file)
        mol2_path = file_path.replace(".sdf", ".mol2")
        if os.path.exists(mol2_path):
            print(f"File {mol2_path} already exists, skipping conversion.")
            continue

        # Convert SDF to MOL2 using Open Babel
        os.system(f"obabel {file_path} -O {mol2_path}")
        print(f"Converted {file} to {mol2_path}")