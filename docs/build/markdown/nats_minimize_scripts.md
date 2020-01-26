# nats_minimize_scripts module


### nats_minimize_scripts.read_pdb_bounds()
Reads at which index each pdb sequence is starting and ending

from the pdb_starts.txt and pdb_ends.txt files from ../data/input/etc


* **Returns**

    Tuple. Dictionaries. Keys: structure pdb id, Values: starting or ending index



* **Return type**

    Tuple



### nats_minimize_scripts.test_missing_res(filename, pdbfile_path)
Detects gaps in residue sequence numbers and split the structure pdb file in two different ones for further

handling with minimization scripts (Charmm)


* **Parameters**

    
    * **filename** – Name for the pdb file being handled


    * **pdbfile_path** – Path for the pdb file



* **Returns**

    Array of the wo split pdb file names



* **Return type**

    str[]
