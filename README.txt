README

This python file reads sequences of DNA and breaks them down into the sequence dependent deformability value of each base pair step of DNA according to the published values in Young et al. 2022 (DOI: 10.3390/life12050759). The code can then output those individual values as a CSV file or plot them as a bar plot, or you can choose to use k-mer averaging approach. The k-mer averaging approach can also output the averaged deformability scores as a CSV file, or plot them as a scatter plot or a circular heatmap.
-------------------------------------------------------------------------------------------------------------------------------------
Dependencies:

Python 3.11.7 and the required external libraries compatible with Python 3.11.7 (numpy, pycirclize, matplotlib).

the software has been tested on python 3.11.7 and Windows 11 and Windows 10.

No required non-standard hardware.

---------------------------------------------------------------------------------------------------------------------------------
Installation:

Download the python file and the required "newStepValuesJSON.json" file, and place them into the same directory. Typical installation time on a "normal" desktop computer is a minute.


-------------------------------------------------------------------------------------------------------------------------------
Demo Instructions:


The sample data within the sample dataset file is a short DNA sequence of a minicircle of DNA. To run the code on the sample data, set the variable "sequence" to equal the DNA sequence from the sample dataset file. As the sequence is circular, the variable "circular" should be set equal to True. 

For this demo, we will look at a bar plot of the individual sequence dependent deformability values of this minicircle compared to the average deformability value of a 
base pair step, and we will look at a circular heatmap of the deformability values when averaged with a k-mer size of 35. Additionally, bases 160-250 on this circular heatmap will be labeled with red.

Individual Base Pair Step Deformability Values Bar plot:

For this analysis, set the variables "circular" and "barchart" to True, and set the variables "averaging", "writeCSV",
"scatterPlotKD", "circularHeatmapPlot", and "chartLabels" to False. Set the variable "whatToPlot" equal to 1. Now run the code.
The bar chart drawn by matplotlib should look identical to the one in the sample dataset file.

35-mer deformability score circular heatmap with bases 160-250 labeled with red:

For this analysis, set the variables "circular", "averaging", "chartLabels", and "circularHeatmapPlot" to True and set the variables "barChart", "writeCSV", and "scatterPlotKD" to False.
Set the variable "sliceSize" to 35. Set the variable "plotLabels" to the list ["other", 0, 160, "black", "label", 160, 251, "red", "other", 251, len(sequence), "black"]. 
Now run the code. The circular heatmap create by pycirclize should look identical to the one in the sample dataset file.

The code should only take a few seconds to run for both of these demos.


----------------------------------------------------------------------------------------------------------------------

Instructions for use.

Set the variable "sequence" to your chosen DNA sequence for analysis, and set the variable "sequenceName" to anything of your choice.
You can then choose the settings for analysis using the variables "whatToPlot," "averaging," "sliceSize," "writeCSV," "circular," "barChart," "scatterPlotKD," "chartLabels",
and "circularHeatmapPlot." For descriptions of what these settings do, please see the comments within the code file itself. The variable "plotLabels" is used to 
label the charts produced by the code file, instructions for how to input the labels are within the code file. 

In order to reproduce the analysis in the paper "DNA deformability defines sequence-dependent capture of E. coli gyrase," you must set the variable "sequence" to be the DNA sequence "AGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTCTAGACCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCGGAGAGACAACTTAAAGAGACTTAAAAGATTAATTTAAAATTTATCAAAAAGAGTATTGACTTAAAGTCTAACCTATAGGATACTTACAGCCATAGAGAGGGATAAGGTGAAATAATAGAATGGTATAATTGCGGCCGAGATCTCCATGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTCCTGTCGTCATATCTACAAGCCATCCCCCCACAGATACGGTAAACTAGCCTCGTTTTTGCATCAGGAAAGCAGA", and use the following settings.

whatToPlot = 2, averaging = False, sliceSize = 35, writeCSV = True, circular = True, barChart = True, scatterPlotKD = False, chartLabels = False, circularHeatMapPlot = False
whatToPlot = 2, averaging = True, sliceSize = 35, writeCSV = True, circular = True, barChart = False, scatterPlotKD = True, chartLabels = False, circularHeatmapPlot = True
whatToPlot = 2, averaging = True, sliceSize = 49, writeCSV = True, circular = True, barChart = False, scatterPlotKD = True, chartLabels = False, circularHeatmapPlot = True

