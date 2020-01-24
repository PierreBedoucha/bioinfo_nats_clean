# nats_dataprocess_loops module


### nats_dataprocess_loops.facet_span(pdbid, y)
Facet span function to apply to the matplotlib facetGrid
Plots the 2 loop areas and the median dotted-line in gray
:param pdbid: pdbids column used as index in the dataframe
:type pdbid: str
:param y: fluctuation values as column in the dataframe
:type y: float


### nats_dataprocess_loops.is_outlier(points, thresh=3.5)
Returns a boolean array with True if points are outliers and False
otherwise.
Boris Iglewicz and David Hoaglin (1993), “Volume 16: How to Detect and
Handle Outliers”, The ASQC Basic References in Quality Control:
Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
:param points: An numobservations by numdimensions array of observations
:param thresh: The modified z-score to use as a threshold. Observations with
a modified z-score (based on the median absolute deviation) greater
than this value will be classified as outliers.
:return: A numobservations-length boolean array.


### nats_dataprocess_loops.read_pdb_loops()
Reading loop locations for the different 7 protein structures.
Reads pdb_loops.txt file from ../data/input/etc containing simple data obtained from pdb structures.
:return: Dictionary. Keys: structure pdb id (str), Values: loop 1 start (int), loop 1 end (int),
loop 2 start (int), loop 2 end (int)
:rtype: dict


### nats_dataprocess_loops.read_pdb_starts()
Reads at which index each pdb sequence is starting from the pdb_starts.txt file from ../data/input/etc
:return: Dictionary. Keys: structure pdb id (str), Values: starting index (ind)
:rtype: dict
