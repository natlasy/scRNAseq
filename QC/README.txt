## scRNAseq.plate.R :
This R script helps getting the overall quality information of the cells in a positional approach in the sorted 384-well plates for single cell RNAseq protocols. The aim is to get an idea on the bias that may happen in the cell sorting or processing steps of the protocol across several plates.
The input is a count matrix contains features (must contain gene symbol) in the rows and the single cells in the columns and a text file (index.plate.txt) contains the cell-barcode numbers in one column from the first to the last well of the plate in rowwise direction.
##
