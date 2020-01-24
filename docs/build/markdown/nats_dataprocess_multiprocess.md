# nats_dataprocess_multiprocess module


### nats_dataprocess_multiprocess.pdb_occupancy()
Cleans the pdb files in the current directory by quickly replacing with its fixed version.
Each cleaned pdb filename is appended to a list to later log the data.


### nats_dataprocess_multiprocess.read_pdb_chains()
Read the selected chains for the protein dataset. The data is parsed from pdb_chains.txt file in
../data/input/etc.
:return: Dictionary. Keys: structure pdb id, Values: selected chain letter
:rtype: dict


### nats_dataprocess_multiprocess.score_proteins(matches, pdb_filename)
Base structure scoring function. Here only the relaxation with FastRelax algorithm is performed and the
relaxed pose is scored and passed through the mp channel to handle replicates.
:param matches: shared list (channel)
:param pdb_filename: pdb file name
