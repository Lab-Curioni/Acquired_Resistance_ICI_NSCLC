import sys

# version adapted from ngs-pipe script, authored by Jochen Singer
inFileName = sys.argv[1]
outFileName = sys.argv[2]
acceptMissingSites = sys.argv[3] # y/n (missing sites, genotype=./. for at least one sample, either take site into account by storing genotype as ./., or discard site by storing NA -> then sites are not used for R heatmap)

acceptMissing = True
if "n" in acceptMissingSites:
	acceptMissing = False


outFile = open(outFileName, 'w')
inFile = open(inFileName, 'r')

for line in inFile:
    if line.startswith("##"):
        continue

    lineSplit = line.strip().split("\t")
    if line.startswith("#"):
        for i in range(9, len(lineSplit)):
            outFile.write(lineSplit[i] + "\t")
        outFile.write("\n")
        continue

    formatSplit = lineSplit[8].split(":")
    try:
        gtIndex = formatSplit.index("GT")
    except:
        print("No GT field for position: " + lineSplit[0] + "\t" + lineSplit[1])
        continue

    for i in range(9, len(lineSplit)):
        sampleSplit = lineSplit[i].split(":")
        gt = sampleSplit[gtIndex]
	if gt == "./." and not acceptMissing:
		gt = "NA"
        outFile.write(gt + "\t")
    outFile.write("\n")


