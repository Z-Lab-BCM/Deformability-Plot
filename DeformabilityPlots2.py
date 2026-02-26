from pycirclize import Circos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch
import json


#set whatToPlot to 0 if you want to plot deformability with no alteration
#set whatToPlot to 1 if you want to plot deformability compared to average of all steps
#set whatToPlot to 2 if you want to plot deformability compared to sequence average deformability
whatToPlot = 2

sequenceName = "ha"

sequence = "AGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTCTAGACCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAATTCGAGCTCGGTACCGGAGAGACAACTTAAAGAGACTTAAAAGATTAATTTAAAATTTATCAAAAAGAGTATTGACTTAAAGTCTAACCTATAGGATACTTACAGCCATAGAGAGGGATAAGGTGAAATAATAGAATGGTATAATTGCGGCCGAGATCTCCATGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCCTGAGTAGGACAAATCCGCCGGGAGCGGATTTGAACGTTGCGAAGCAACGGCCCGGAGGGTGGCGGGCAGGACGCCCGCCATAAACTGCCAGGCATCAAATTAAGCAGAAGGCCATCCTGACGGATGGCCTTTTTGCGTTTCTACAAACTCTTCCTGTCGTCATATCTACAAGCCATCCCCCCACAGATACGGTAAACTAGCCTCGTTTTTGCATCAGGAAAGCAGA"

averaging = True # if this is true the data will be split into k-mers of your choosing
sliceSize = 35 #Size of the k-mer, if averaging is false, this is automatically set to 2
writeCSV = False #outputs a CSV with the results of the analysis.
circular = True #set to True if your sequence is circular
barChart = True #this will be displayed automatically when the code runs
scatterPlotKD = True #this will be displayed automatically when the code runs
chartLabels = False #see below for explanation of how to make the labels
circularHeatmapPlot = False #this will be added to the current directory as a png file
literatureAverageOfSteps =  3.6896263958334 #average of the deformability values for the 136 unique tetramers

plotLabels = []

'''
The above list "plotLabels" is how you create the bar plot and label segments of the bar plot with colors and names.
If you are not going to do any labeling with colors set chartLabels = False

If you are going to have labels you must follow this pattern ["Name of label", label start (int), label end (int), "color"]. 
Repeat this pattern for as many segments as you want, but it must always be in this pattern!

The end of one segment must equal the start of the next segment.
if you don't want a segment to be labeled, its name must be "other"

keep in mind that python indexing starts at 0, and that slicing in python is inclusive of the start value but exclusive of the end value.

the end value of the final segment MUST be len(sequence)

'''

sequence = sequence.upper()




if chartLabels == False:
    plotLabels = ["other", 0 , len(sequence), "black"]
sliceValues = []

if averaging == False:
    sliceSize = 2

if writeCSV == True:
    outfile = open(f"{sequenceName}slice{sliceSize}CSV.csv", "w")
    outfile.write("Center Base Pair Number,Sequence,Deformability Score\n")

stepvaluefile = open("newstepValuesJSON.json", "r")
stepValueDict = json.load(stepvaluefile)
stepvaluefile.close()

def circularSequence(sequence):
    if sliceSize%2==1:
        firsthalf = sequence[0:(sliceSize//2+1)]
        lasthalf = sequence[len(sequence) - (sliceSize // 2+1):len(sequence)]
    elif sliceSize%2==0:
        firsthalf = sequence[0:(sliceSize//2+1)]
        lasthalf = sequence[len(sequence)-(sliceSize//2):len(sequence)]
    circseq = lasthalf + sequence + firsthalf
    return circseq

def deformability(seq, nSkip=False):
    seq=str(seq).upper()
    sliceValues = []
    for i in range(1, len(seq)-2):
        tetramer = seq[i-1:i+3]
        if ("N" in tetramer or "M" in tetramer or "R" in tetramer) and nSkip == True:
            continue
        stepValue = stepValueDict[tetramer]
        sliceValues.append(stepValue)
        if writeCSV == True and averaging == False:
            outfile.write(f"{i},{tetramer},{stepValue}\n")
    return sliceValues

def slicing(seq):
    slicesAverages = []
    slicesCenters = []
    slices = []
    slidingSlice = []
    if circular:
        for i in range(0, len(seq)-sliceSize+1):
            slices.append(seq[i:i + sliceSize])
            if i == 0:
                slidingSlice = sliceValues[0:sliceSize-2]
                continue
            elif i == 1:
                slidingSlice.append(sliceValues[sliceSize-2])
            elif i == len(seq)-sliceSize:
                continue
            else:
                del slidingSlice[0]
                slidingSlice.append(sliceValues[sliceSize-3+i])
            slicesAverages.append(np.mean(slidingSlice))
            slicesCenters.append(i+sliceSize//2)
            if writeCSV:
                if sliceSize%2==1:
                    outfile.write(f"{i},{seq[i:i+sliceSize]},{np.mean(slidingSlice)}\n")
                elif sliceSize%2==0:
                    outfile.write(f"{i},{seq[i:i+sliceSize]},{np.mean(slidingSlice)}\n")
    if circular == False:
        for i in range(0, len(seq)-sliceSize+1):
            if i == 0:
                slidingSlice = sliceValues[0:sliceSize-2]
            elif i == 1:
                slidingSlice.append(sliceValues[sliceSize-2])
            elif i == len(seq)-sliceSize:
                del slidingSlice[0]
            else:
                del slidingSlice[0]
                slidingSlice.append(sliceValues[sliceSize-3+i])
            slicesAverages.append(np.mean(slidingSlice))
            slicesCenters.append(i+sliceSize//2)
            if writeCSV:
                if sliceSize%2==1:
                    outfile.write(f"{i+sliceSize//2},{seq[i:i+sliceSize]},{np.mean(slidingSlice)}\n")
                elif sliceSize%2==0:
                    outfile.write(f"{i+sliceSize//2-1},{seq[i:i+sliceSize]},{np.mean(slidingSlice)}\n")
    return slicesAverages, slicesCenters, slices

def yValuePlotBar(yvalsinput):
    yvalsarr = np.array(yvalsinput)
    if whatToPlot == 0:
        return yvalsarr
    elif whatToPlot == 1:
        return yvalsarr - literatureAverageOfSteps
    elif whatToPlot == 2:
        return yvalsarr - np.mean(yvalsarr)

def simpleBarChart(axis, xvalsbar, yvalsbar):
    for i in range(1, len(plotLabels), 4):
        length = plotLabels[i+1] - plotLabels[i]
        if whatToPlot == 0:
            axis.bar(xvalsbar[plotLabels[i]:plotLabels[i+1]], yvalsbar[plotLabels[i]:plotLabels[i+1]], width=1, color=plotLabels[i+2])
        if whatToPlot == 1:
            axis.bar(xvalsbar[plotLabels[i]:plotLabels[i+1]], yvalsbar[plotLabels[i]:plotLabels[i+1]], width=1, color=plotLabels[i+2], bottom=literatureAverageOfSteps)
        if whatToPlot == 2:
            print(np.mean(sliceValues))
            axis.bar(xvalsbar[plotLabels[i]:plotLabels[i+1]], yvalsbar[plotLabels[i]:plotLabels[i+1]], width=1, color=plotLabels[i+2], bottom=np.mean(sliceValues))
        if plotLabels[i-1] != "other" and plotLabels[i+2] != "red":
            axis.add_patch(patch.Rectangle((plotLabels[i]+1, -1), length, 1, facecolor=plotLabels[i+2], label=plotLabels[i-1]))
            axis.text(plotLabels[i]+1+length/2, -0.5, plotLabels[i-1], ha="center", va="center", fontsize=8)
    axis.legend(handles= patches, loc="upper left")


def simpleScatterPlot(axis):
    for i in range(1, len(plotLabels), 4):
        axis.plot(xvalsScat[plotLabels[i]:plotLabels[i+1]+1], slicesAverages[plotLabels[i]:plotLabels[i+1]+1], color=plotLabels[i+2], marker=".", linestyle="-")
    if whatToPlot == 2:
        axis.axhline(y=np.mean(sliceValues), color="blue", label="Average Deformability")
    axis.legend(handles=patches, loc="upper left")

def xValuePlot(input):
    xvalsarr = np.arange(1, len(input)+1)
    xvalsarr = xvalsarr.astype(float)
    xvalsarr += 0.5
    return xvalsarr

def labelPatches():
    patches = []
    for i in range(0, len(plotLabels), 4):
        if plotLabels[i] == "other":
            continue
        patchy = patch.Patch(color=plotLabels[i + 3], label=plotLabels[i])
        patches.append(patchy)
    return patches

def circularHeatmap(slicesAverages):
    sectors = {"A": len(sequence)}
    circos = Circos(sectors, space=0)
    patches = labelPatches()
    for sector in circos.sectors:
        sectorName = f"{sequenceName} Slice Size {sliceSize} Deformability Heat Map"
        data = slicesAverages
        heatmap_track = sector.add_track((75, 100))
        heatmap_track.xticks_by_interval(50, show_endlabel=False, label_size=15)
        heatmap_track.xticks_by_interval(10, tick_length=1, show_label=False)
        heatmap_track.axis(fc="none")
        heatmap_track.heatmap(data, cmap="cividis")#, vmin=2.2463, vmax=4.833)
    circos.colorbar(bounds=(1, 0, 0.02, 0.3), vmin=min(slicesAverages), vmax=max(slicesAverages), cmap="cividis", label=f"{sliceSize}-mer Deformability Score")
    for i in range(0, len(plotLabels), 4):
        degStart = (plotLabels[i+1]/len(sequence))*360
        degEnd = (plotLabels[i+2]/len(sequence))*360
        if plotLabels[i] == "other":
            continue
        circos.rect(r_lim=(70,74), deg_lim=(degStart,degEnd), fc=plotLabels[i+3])

    def degstartstop(a, b):
        start = (a/len(sequence))*360
        stop = (b/len(sequence))*360
        return start,stop

    circos.plotfig()
    circos.ax.legend(handles=patches, bbox_to_anchor=(0.5, 0.5), loc="center", fontsize=12)
    plt.savefig(f"{sequenceName}slice{sliceSize}heatmap.svg", format="svg")



patches = labelPatches()

if circular:
    circSeq = circularSequence(sequence)
    sliceValues = deformability(circSeq)
    if averaging:
        slicesAverages, slicesCenters, slices = slicing(circSeq)
if circular == False:
    sliceValues = deformability(sequence)
    if averaging:
        slicesAverages, slicesCenters, slices = slicing(sequence)

if writeCSV == True:
    outfile.close()

if averaging == False and barChart == True:
    xvalsbar = xValuePlot(sliceValues)
    yvalsbar = yValuePlotBar(sliceValues)
    stepsGraph, ax = plt.subplots()
    ax.set_xlabel("Base Pair Step", fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.margins(x=0.01)
    simpleBarChart(ax, xvalsbar, yvalsbar)
    if whatToPlot == 0:
        ax.set_ylabel("Step Deformability", fontsize=15)
    if whatToPlot == 1:
        ax.set_ylabel("Step Deformability Compared to Average Deformability of a Base Pair Step", fontsize=15)
    if whatToPlot == 2:
        ax.set_ylabel("Step Deformability Compared to Average Deformability\nof Sequence", fontsize=12)
    ax.tick_params(axis="both", labelsize=5)
    stepsGraph.savefig(f"{sequenceName}_deformabilityProfile.svg", format="svg", bbox_inches="tight")
    plt.show()

if averaging == True and scatterPlotKD == True:
    scatterGraph = plt.figure()
    scatAx = scatterGraph.add_subplot()
    if circular:
        xvalsScat = xValuePlot(slicesAverages)
    elif circular == False:
        xvalsScat = np.array(slicesCenters)
    scatAx.set_xlabel("Base Pair Step", fontsize=20)
    simpleScatterPlot(scatAx)
    scatAx.tick_params(axis="both", which="major", labelsize=16)
    if whatToPlot == 0:
        scatAx.set_ylabel("Deformability Score")
    elif whatToPlot == 1:
        scatAx.set_ylabel("Deformability Score Compared to Average Deformability of a Base Pair Step")
    elif whatToPlot == 2:
        scatAx.set_ylabel("Deformability Score Compared to Average Deformability \nof Sequence", fontsize=12)
    plt.savefig(f"{sequenceName}slice{sliceSize}_deformabilityProfile.svg", format="svg")
    plt.show()


if averaging == True and circular == True and circularHeatmapPlot == True:
    circularHeatmap(slicesAverages)
