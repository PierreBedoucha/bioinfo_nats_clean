# nats_dataprocess_multiprocess module

Energy scoring of the NATSs protein models. MULTIPROCESS

This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API
The scoring function parameters are detailed in the score_proteins function. The structures are scored before and
after a relaxation step done with the FastRelax algorithm.

The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.


### nats_dataprocess_multiprocess.pdb_occupancy()
Cleans the pdb files in the current directory by quickly replacing with its fixed version.

Each cleaned pdb filename is appended to a list to later log the data.


### nats_dataprocess_multiprocess.read_pdb_chains()
Read the selected chains for the protein dataset.

The data is parsed from pdb_chains.txt file in ../data/input/etc.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: selected chain letter



* **Return type**

    dict



### nats_dataprocess_multiprocess.score_proteins(matches, pdb_filename)
Base structure scoring function.

Here only the relaxation with FastRelax algorithm is performed and the relaxed pose is scored and passed through
the mp channel to handle replicates.


* **Parameters**

    
    * **matches** – shared list (channel)


    * **pdb_filename** – pdb file name
