# nats_dataprocess_mpi module


### nats_dataprocess_mpi.pdb_occupancy()
Cleans the pdb files in the current directory by quickly replacing with its fixed version and launches the scoring.
Each cleaned pdb filename is appended to a list to later log the data.


### nats_dataprocess_mpi.read_pdb_chains()
Read the selected chains for the protein dataset. The data is parsed from pdb_chains.txt file in
../data/input/etc.
:return: Dictionary. Keys: structure pdb id, Values: selected chain letter
:rtype: dict


### nats_dataprocess_mpi.score_proteins(pdb_filename)
Structure scoring function. Launches the main scoring function on global number of replicates
:param pdb_filename: pdb file name
