# subplot_rmsd_score_stats module

This script handles PROCHECK analyses for the structure files in current directory.

Be careful, reference pdb files (amplitude 0.00 for both native and relaxed structures) are needed for comparison.

The PROCHECK files are created and placed in the current directory for later use.


### subplot_rmsd_score_stats.calc_distance(pdbfile_line_1, pdbfile_line_2)
Calculate in line distance (Angstroms) between two atoms.

param pdbfile_line_1: Str line for atom line in first pdb file


* **Parameters**

    **pdbfile_line_2** – Str line for atom line in second pdb file



* **Returns**

    Distance array from distance.cdist method



* **Return type**

    float[][]



### subplot_rmsd_score_stats.compute_rmsd(pdb_path1, pdb_path2, \*\*kwargs)
Computes RMS distance between two pdb structures and only from start to end indices in their sequence


* **Parameters**

    
    * **pdb_path1** (*str*) – First pdb file name


    * **pdb_path2** (*str*) – First pdb file name


    * **kwargs** – Keyword arguments with optional start and end in the pdb sequence instead of obtaining the
    boundaries from reading the msa data.



* **Returns**

    Rmsd value between pdb structures



* **Return type**

    float



### subplot_rmsd_score_stats.f(x, y, z, \*\*kwargs)
Add annotation on seaborn plot


* **Parameters**

    
    * **x** – x value for the annotation position on plot


    * **y** – y value for the annotation poistion on plot


    * **z** – Annotated value


    * **kwargs** – keyword arguments



### subplot_rmsd_score_stats.facet_line()
Draw scatterplot with point colors from a faceted DataFrame columns.


### subplot_rmsd_score_stats.facet_scatter(x, y, c, \*\*kwargs)
Draw scatterplot with point colors from a faceted DataFrame columns.


* **Parameters**

    
    * **x** – x axis data


    * **y** – y axis data


    * **kwargs** – keyword arguments



### subplot_rmsd_score_stats.facet_stairs(x, y, \*\*kwargs)
Draw scatterplot with point colors from a faceted DataFrame columns.


* **Parameters**

    
    * **x** – x axis data


    * **y** – y axis data


    * **kwargs** – keyword arguments



### subplot_rmsd_score_stats.procheck(sc_pdbpath, pdb_resolution)
Main PROCHECK computation function. The full atom (sc) structure pb file is required as well as its corresponding
resolution (Angstroms). All the resulting files are created in the current directory for later use.


* **Parameters**

    
    * **sc_pdbpath** (*str*) – Full atom pdb file path


    * **pdb_resolution** (*float*) – Resolution of the structure (Angstroms)



* **Returns**

    Dictionary. Keys: PROCHECK measurement types, Values: measures



* **Return type**

    dict



### subplot_rmsd_score_stats.read_msa_fasta()
Reads multiple structure alignment from MUSTANG.
It determines the structurally aligned core of the proteins.

Note: here, only the aligned regions are of interest, gaps are removed.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: aligned indices



* **Return type**

    dict



### subplot_rmsd_score_stats.read_pdb_resolution()
Reads the resolution from data for the pdb files.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: resolution (float)



* **Return type**

    dict



### subplot_rmsd_score_stats.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: starting index



* **Return type**

    dict



### subplot_rmsd_score_stats.read_pfam_align()
Reads multiple sequence alignment profile from pfam.

It determines the envelope (start and end) of the GNAT fold

Reads values from pfam_env.txt file in ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: envelope indices (tuple)



* **Return type**

    dict
