# nats_plot_residue_score module


### nats_plot_residue_score.read_msa_fasta()
Reads multiple structure alignment from MUSTANG. It determines the structurally aligned core of the proteins.
Note: here, only the aligned regions are of interest, gaps are removed.
:return: Dictionary. Keys: structure pdb id, Values: aligned indices
:rtype: dict


### nats_plot_residue_score.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: starting index
:rtype: dict
