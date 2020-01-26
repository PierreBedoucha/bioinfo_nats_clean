# nats_dataprocess_neighbors module

Analysis for active sites neigbors in NATs

This script summarizes the energy score of specific residues in the available structures (relaxed pdb file)
in the current directory.
The structures scores are listed in the relaxed pdb files and the script will isolate only the resiudes neighbouring
the catalytic sites. The resulting score is summed over this selection.

The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.


### nats_dataprocess_neighbors.neighbor_res_select(pdbfilename, sel_chain, cutoff)
Selection of the catalytic residue neighbours from the structure pdb files.

It uses NeighborSearch method from BioPDB. All the detected residues indices are returned in a list.


* **Parameters**

    
    * **pdbfilename** (*str*) – Pdb file name


    * **sel_chain** (*str*) – Identification letter of the structure chain


    * **cutoff** (*float*) – Cutoff value in Angstroms for neighbouring search



* **Returns**

    List of detected residue indices



* **Return type**

    list



### nats_dataprocess_neighbors.read_msa_fasta()
Reads multiple structure alignment from MUSTANG.

It determines the structurally aligned core of the proteins.
Note: here, only the aligned regions are of interest, gaps are removed.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: aligned indices



* **Return type**

    dict



### nats_dataprocess_neighbors.read_pdb_catalytic()
Reads the list of catalytic residues listed in pdb_catalytic.txt file in ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: catalytic residue index



* **Return type**

    dict



### nats_dataprocess_neighbors.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: starting index



* **Return type**

    dict
