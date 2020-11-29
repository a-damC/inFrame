# inFrame
*pipeline used to analyse phase variable data*

The inFrame pipeline is used in the following order:
1. [inFrame_extract.py](https://github.com/a-damC/inFrame/blob/main/inFrame_extract.py) is used to extract coordinates of a gene from a sequence file
2. [inFrane_protCheck.py](https://github.com/a-damC/inFrame/blob/main/inFrame_protCheck.py) is used to see if the extracts produce a viable protein
3. [pv_dist_matrix.py](https://github.com/a-damC/inFrame/blob/main/pv_dist_matrix.py) takes the PV state data of different isolates and compares them
4. [mannWhitney.py](https://github.com/a-damC/inFrame/blob/main/mannWhitney.py) is a rank sum statistical test to determine if the isolates cluster the same way in the PV data as they do in allelic data
