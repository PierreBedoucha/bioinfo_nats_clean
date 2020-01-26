# NATs Bioinformatics 2019

This repository contains the data analysis part of a project on NATs. It comprises a set of python scripts, independent from one another, and bound to the pdb structure files of the proteins in the dataset.


### 1. Setup

Below are som instructions on how to set everything up.

* Install [Anaconda](https://www.anaconda.com/download/) on your laptop, choose the Python 3.6 version.
* Run ```$ conda env create -f environment.yml``` from the provided file in /docs folder.

### 2. Documentation

Below is the autodoc generated documentation from sphinx Markdown

#### intro


* PyRosetta_TACC_MPI module
        More information at the following link [Using PyRosetta on TACC](https://www.rosettacommons.org/docs/latest/build_documentation/TACC)

* facet_boxplot_rmsd_score_stats module

        Plotting Energy Score and RMSD in Facet Boxplots

        The script analyses the obtained Rosetta Energy Scores for all the structures of the dataset (12 pdbs) for the specific
        case when 5 replicates have been computed for 4 negative amplitudes.

        It produces a facet grid, for the first mode 7 only, of the boxplots per amplitudes.

        The energy score is normalized per structure residue.

* nats_dataprocess module
        
        Energy scoring of the NATSs protein models

        This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API

        The scoring function parameters are detailed in the score_proteins function. The structures are scored before and
        after a relaxation step done with the FastRelax algorithm.

        The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.

* nats_dataprocess_loops module
  
        Analysis for Fluctuations of the loops in NATs

        The script analyses the normalized squared fluctuation data for the first 6 normal modes (successively).

        It plots the fluctuation graphs for the 7 protein structures (see data), adds the mean value for baseline description, and highlights the 2 loop regions.

        The modes of interest are the ones with the highest fluctuation values for the 2 functional loops.

* nats_dataprocess_mpi module

        Energy scoring of the NATSs protein models. MPI TEST

        This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API
        The scoring function parameters are detailed in the _main function.

        The structures are scored before and after a relaxation step done with the FastRelax algorithm.

        The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.

* nats_dataprocess_multiprocess module

        Energy scoring of the NATSs protein models. MULTIPROCESS

        This script computes the energy score of available structures (Pose object) in the current directory using pyROSETTA API
        The scoring function parameters are detailed in the score_proteins function. The structures are scored before and
        after a relaxation step done with the FastRelax algorithm.

        The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.

* nats_dataprocess_neighbors module

        Analysis for active sites neigbors in NATs

        This script summarizes the energy score of specific residues in the available structures (relaxed pdb file)
        in the current directory.
        The structures scores are listed in the relaxed pdb files and the script will isolate only the resiudes neighbouring
        the catalytic sites. The resulting score is summed over this selection.

        The resulting data is logged in a csv file ‘pyrosetta_out.csv’ for further analyses.

* nats_minimize_scripts module

        Charmm scripts for dataset structure minimization

        This script handles the creation of a collection of shell scripts for file transfer to run minimization steps
        on server (doggpil).

        It inputs the selected structure files (pdb) from the current directory.

* nats_plot_residue_score module

        This script plot the energy score (from pyROSETTA API) over all the structure files and per aligned residues.

        The structures scores are listed in the relaxed pdb files and the script will isolate it for all the residues.

        The final results is summed through all the current directory structure files.

* subplot_rmsd_score module

        This script plots the lines benchmarking the energy score (pyROSETTA API) per structure

        in the current directory and after their relaxation and per amplitudes of deformation.

        The RMSD is overlaid as continuous color scale for each datapoint. The energy score is normalized per structure residue.

* subplot_rmsd_score_stats module

        This script handles PROCHECK analyses for the structure files in current directory.

        Be careful, reference pdb files (amplitude 0.00 for both native and relaxed structures) are needed for comparison.

        The PROCHECK files are created and placed in the current directory for later use.


