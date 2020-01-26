# nats_plot_residue_score module

This script plot the energy score (from pyROSETTA API) over all the structure files and per aligned residues.

The structures scores are listed in the relaxed pdb files and the script will isolate it for all the residues.

The final results is summed through all the current directory structure files.


### nats_plot_residue_score.read_msa_fasta()
Reads multiple structure alignment from MUSTANG.
It determines the structurally aligned core of the proteins.

Note: here, only the aligned regions are of interest, gaps are removed.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: aligned indices



* **Return type**

    dict



### nats_plot_residue_score.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: starting index



* **Return type**

    dict
