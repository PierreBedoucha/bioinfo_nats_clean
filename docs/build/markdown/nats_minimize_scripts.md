# nats_minimize_scripts module


### nats_minimize_scripts.read_pdb_ends()
Reads at which index each pdb sequence is starting from the pdb_ends.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: ending index
:rtype: dict


### nats_minimize_scripts.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: starting index
:rtype: dict


### nats_minimize_scripts.test_missing_res(filename, pdbfile_path)
Detects gaps in residue sequence numbers and split the structure pdb file in two different ones for further
handling withe minimization scripts (Charmm)
:param filename: Name for the pdb file being handled
:param pdbfile_path: Path for the pdb file
:return: Array of the wo split pdb file names
:rtype: str[]
