# nats_dataprocess_loops module

Analysis for Fluctuations of the loops in NATs

The script analyses the normalized squared fluctuation data for the first 6 normal modes (successively).

It plots the fluctuation graphs for the 7 protein structures (see data), adds the mean value for baseline description, and highlights
the 2 loop regions.

The modes of interest are the ones with the highest fluctuation values for the 2 functional loops.


### nats_dataprocess_loops.facet_span(pdbid, y)
Facet span function to apply to the matplotlib facetGrid

Plots the 2 loop areas and the median dotted-line in gray


* **Parameters**

    
    * **pdbid** (*str*) – pdbids column used as index in the dataframe


    * **y** (*float*) – fluctuation values as column in the dataframe



### nats_dataprocess_loops.is_outlier(points, thresh=3.5)
Returns a boolean array with True if points are outliers and False
otherwise.

Boris Iglewicz and David Hoaglin (1993), “Volume 16: How to Detect and
Handle Outliers”, The ASQC Basic References in Quality Control:
Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.


* **Parameters**

    
    * **points** – An numobservations by numdimensions array of observations


    * **thresh** – The modified z-score to use as a threshold. Observations with
    a modified z-score (based on the median absolute deviation) greater
    than this value will be classified as outliers.



* **Returns**

    A numobservations-length boolean array.



### nats_dataprocess_loops.read_pdb_loops()
Reading loop locations for the different 7 protein structures.

Reads pdb_loops.txt file from ../data/input/etc containing simple data obtained from pdb structures.


* **Returns**

    Dictionary. Keys: structure pdb id (str), Values: loop 1 start (int), loop 1 end (int),
    loop 2 start (int), loop 2 end (int)



* **Return type**

    dict



### nats_dataprocess_loops.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc


* **Returns**

    Dictionary. Keys: structure pdb id (str), Values: starting index (ind)



* **Return type**

    dict
