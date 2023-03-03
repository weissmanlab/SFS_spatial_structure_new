from os import path

inputVal = "mixed_highR.txt"
outputStart = "mixed_highR_"
mapLoci = {}
lociList = []
flagged = set()
jitter = 0.00001

def writeOutput(position, mutationCount, outputFileName, nbase):
    if path.exists(outputFileName) == False:
        output = open(outputFileName, "w")
        output.write("position\t" + "x\t" + "n\t" + "folded\n")
        output.close()
    output = open(outputFileName, "a")
    output.write(str(position) + "\t" + str(mutationCount) + "\t" + str(nbase) + "\t" + "0\n")
    output.close()

#Compiles allele frequency data from input file.
x = 0
with open(inputVal, "r") as file:
    listLocation = 0
    for line in file:
        #Identify all the segregating sites
        if line[0:10] == "positions:":
            splitLine = line[12:]
            currentNum = ""
            for character in splitLine:
                if character != " " and character != "\n" and character != "":
                    currentNum += character
                else:
                    if currentNum == "":
                        break
                    cast = float(currentNum)
                    currentNum = ""
                    lociList.append(cast)
        #Identify the counts at the i-th segsite
        elif line[0] == "1" or line[0] == "0":
            for character in line:
                if character == "\n":
                    listLocation = 0
                    flagged = set()
                elif character == "1":
                    position = lociList[listLocation]
                    while position in flagged:
                        position += jitter
                    count = mapLoci.get(position)
                    if count is not None:
                        mapLoci[position] = mapLoci[position] + 1
                    else:
                        mapLoci[position] = 1
                    listLocation += 1
                    flagged.add(position)
                else:
                    listLocation += 1
        #Write output and clear data structures for next trial
        elif line[0] == "/":
            if x == 0:
                x += 1
                continue
            output = outputStart + str(x) + ".txt"
            for location in mapLoci:
                
                writeOutput(int(location * 1000000), mapLoci[location], output, 100)
            lociList = []
            mapLoci = {}
            flagged = set()
            x += 1
        else:
            continue
        
#If there is any trial data left, process it.
if len(mapLoci) > 0:
    output = outputStart + str(x) + ".txt"
    for location in mapLoci:
        writeOutput(location, mapLoci[location], output, 100)
