# nats_dataprocess_mpi module

Energy scoring of the NATSs protein models. MPI TEST

This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API
The scoring function parameters are detailed in the _main function.

The structures are scored before and after a relaxation step done with the FastRelax algorithm.

The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.


### nats_dataprocess_mpi.pdb_occupancy()
Cleans the pdb files in the current directory by quickly replacing with its fixed version and launches the scoring.

Each cleaned pdb filename is appended to a list to later log the data.


### nats_dataprocess_mpi.read_pdb_chains()
Read the selected chains for the protein dataset.

The data is parsed from pdb_chains.txt file in ../data/input/etc.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: selected chain letter



* **Return type**

    dict



### nats_dataprocess_mpi.score_proteins(pdb_filename)
Structure scoring function.

Launches the main scoring function on global number of replicates


* **Parameters**

    **pdb_filename** – pdb file name
