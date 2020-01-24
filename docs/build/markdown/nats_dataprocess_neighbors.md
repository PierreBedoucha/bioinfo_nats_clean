# nats_dataprocess_neighbors module


### nats_dataprocess_neighbors.neighbor_res_select(pdbfilename, sel_chain, cutoff)
Selection of the catalytic residue neighbours from the structure pdb files. It uses NeighborSearch method from
BioPDB. All the detected residues indices are returned in a list.
:param pdbfilename: Pdb file name
:type pdbfilename: str
:param sel_chain: Identification letter of the structure chain
:type sel_chain: str
:param cutoff: Cutoff value in Angstroms for neighbouring search
:type cutoff: float
:return: List of detected residue indices
:rtype: list


### nats_dataprocess_neighbors.read_msa_fasta()
Reads multiple structure alignment from MUSTANG. It determines the structurally aligned core of the proteins.
Note: here, only the aligned regions are of interest, gaps are removed.
:return: Dictionary. Keys: structure pdb id, Values: aligned indices
:rtype: dict


### nats_dataprocess_neighbors.read_pdb_catalytic()
Reads the list of catalytic residues listed in pdb_catalytic.txt file in ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: catalytic residue index
:rtype: dict


### nats_dataprocess_neighbors.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: starting index
:rtype: dict
