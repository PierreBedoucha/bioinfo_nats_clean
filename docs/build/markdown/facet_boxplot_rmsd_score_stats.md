# facet_boxplot_rmsd_score_stats module

Plotting Energy Score and RMSD in Facet Boxplots

The script analyses the obtained Rosetta Energy Scores for all the structures of the dataset (12 pdbs) for the specific
case when 5 replicates have been computed for 4 negative amplitudes.

It produces a facet grid, for the first mode 7 only, of the boxplots per amplitudes.

The energy score is normalized per structure residue.


### facet_boxplot_rmsd_score_stats.calc_distance(pdbfile_line_1, pdbfile_line_2)
Calculate in line distance (Angstroms) between two atoms.


* **Parameters**

    
    * **pdbfile_line_1** – Str line for atom line in first pdb file


    * **pdbfile_line_2** – Str line for atom line in second pdb file



* **Returns**

    Distance array from distance.cdist method



* **Return type**

    float[][]



### facet_boxplot_rmsd_score_stats.compute_rmsd(pdb_path1, pdb_path2, \*\*kwargs)
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



### facet_boxplot_rmsd_score_stats.f(x, y, z, \*\*kwargs)
Add annotation on seaborn plot


* **Parameters**

    
    * **x** – x value for the annotation position on plot


    * **y** – y value for the annotation poistion on plot


    * **z** – Annotated value


    * **kwargs** – keyword arguments



### facet_boxplot_rmsd_score_stats.facet_scatter(x, y, \*\*kwargs)
Draw scatterplot with point colors from a faceted DataFrame columns.


* **Parameters**

    
    * **x** – x axis data


    * **y** – y axis data


    * **kwargs** – keyword arguments



### facet_boxplot_rmsd_score_stats.read_msa_fasta()
Reads multiple structure alignment from MUSTANG.
It determines the structurally aligned core of the proteins.

Note: here, only the aligned regions are of interest, gaps are removed.


* **Returns**

    Dictionary. Keys: structure pdb id, Values: aligned indices



* **Return type**

    dict



### facet_boxplot_rmsd_score_stats.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values: starting index



* **Return type**

    dict



### facet_boxplot_rmsd_score_stats.read_pfam_align()
Reads multiple sequence alignment profile from pfam.
It determines the envelope (start and end) of the GNAT fold.

Reads values from pfam_env.txt file in ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id, Values:



* **Return type**

    dict



### facet_boxplot_rmsd_score_stats.select_ca_align(pdb_path, start, end)
Get the list of Calpha atom indices in the input structure and between start and end sequence indices


* **Parameters**

    
    * **pdb_path** (*str*) – First pdb file name


    * **start** (*int*) – Start sequence index for rmsd computation


    * **end** (*int*) – End sequence index for rmsd computation



* **Returns**

    List of selected Calpha atom indices



* **Return type**

    list
