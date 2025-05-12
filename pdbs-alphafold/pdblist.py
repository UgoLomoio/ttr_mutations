import os 

cwd = os.getcwd()
sep = os.sep

pdb_path = cwd + sep + "pdbs-alphafold"
complex_types = ["monomer", "dimer", "tetramer"]
filenames4pcn = cwd + sep + "filenames4pcn.txt"

if __name__ == "__main__":
	with open(filenames4pcn, "w") as file:
		for complex_type in complex_types:
			directory = pdb_path + sep + complex_type
			filenames = [file for file in os.listdir(directory)]
			filenames_str = ""
			for filename in filenames:
				if filename != filenames[-1]:
					filenames_str+= "{}, ".format(filename)
				else:
					filenames_str+= filename
			line = "{}: {}\n".format(complex_type, filenames_str)
			file.write(line)
			print(line)
	