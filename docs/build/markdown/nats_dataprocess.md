# nats_dataprocess module


### nats_dataprocess.pdb_occupancy()
Cleans the pdb files in the current directory by quickly replacing its fixed version and launches the scoring
Each cleaned pdb filename is appended to a list to later log the data.


### nats_dataprocess.read_pdb_chains()
Read the selected chains for the protein dataset. The data is parsed from pdb_chains.txt file in
../data/input/etc.
:return: Dictionary. Keys: structure pdb id, Values: selected chain letter
:rtype: dict


### nats_dataprocess.replace(file_path, pattern, subst)
Small helper function to replace a str pattern by another one in a given file.
:param file_path: File path of the file to consider
:param pattern: str pattern to replace
:param subst: str pattern for replacing


### nats_dataprocess.score_proteins(pdb_filename)
Main structure scoring function. Describes the scoring parameters and set the values of two global lists containing
the initial score before relaxation on one hand, and the final score after relaxation on the other.
:param pdb_filename: pdb file name
