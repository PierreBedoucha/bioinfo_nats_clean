# facet_boxplot_rmsd_score_stats module


### facet_boxplot_rmsd_score_stats.calc_distance(pdbfile_line_1, pdbfile_line_2)
Calculate in line distance (Angstroms) between two atoms.
:param pdbfile_line_1: Str line for atom line in first pdb file
:param pdbfile_line_2: Str line for atom line in second pdb file
:return: Distance array from distance.cdist method
:rtype: float[][]


### facet_boxplot_rmsd_score_stats.compute_rmsd(pdb_path1, pdb_path2, \*\*kwargs)
Computes RMS distance between two pdb structures and only from start to end indices in their sequence
:param pdb_path1: First pdb file name
:type pdb_path1: str
:param pdb_path2: First pdb file name
:type pdb_path2: str
:param kwargs: Keyword arguments with optional start and end in the pdb sequence instead of obtaining the
boundaries from reading the msa data.
:return: Rmsd value between pdb structures
:rtype: float


### facet_boxplot_rmsd_score_stats.f(x, y, z, \*\*kwargs)
Add annotation on seaborn plot
:param x: x value for the annotation position on plot
:param y: y value for the annotation poistion on plot
:param z: Annotated value
:param kwargs: keyword arguments


### facet_boxplot_rmsd_score_stats.facet_scatter(x, y, \*\*kwargs)
Draw scatterplot with point colors from a faceted DataFrame columns.
:param x: x axis data
:param y: y axis data
:param kwargs: keyword arguments


### facet_boxplot_rmsd_score_stats.read_msa_fasta()
Reads multiple structure alignment from MUSTANG. It determines the structurally aligned core of the proteins.
Note: here, only the aligned regions are of interest, gaps are removed.
:return: Dictionary. Keys: structure pdb id, Values: aligned indices
:rtype: dict


### facet_boxplot_rmsd_score_stats.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values: starting index
:rtype: dict


### facet_boxplot_rmsd_score_stats.read_pfam_align()
Reads multiple sequence alignment profile from pfam. It determines the envelope (start and end) of the GNAT fold
Reads values from pfam_env.txt file in ../data/input/etc
:return: Dictionary. Keys: structure pdb id, Values:
:rtype: dict


### facet_boxplot_rmsd_score_stats.select_ca_align(pdb_path, start, end)
Get the list of Calpha atom indices in the input structure and between start and end sequence indices
:param pdb_path: First pdb file name
:type pdb_path: str
:param start: Start sequence index for rmsd computation
:type start: int
:param end: End sequence index for rmsd computation
:type end: int
:return: List of selected Calpha atom indices
:rtype: list
