import os
import sys

def create_lines_specific_hitran_file(hitFileIn,
                                      freqs):
    """
    Pick out the closest lines to inputted frequencies, and write them out to
    a new file.
    """

    #Check input types.
    if not isinstance(hitFileIn,str):
        raise TypeError("the inputted hitran file (" + repr(hitFileIn) +
                            ") must be a string.\n")
    if not isinstance(freqs,list) and not isinstance(freqs,set):
        raise TypeError("the inputted frequencies must be a list or a set.\n")
    if len(freqs) < 1:
        raise TypeError("the inputted frequencies list is empty.\n")
    for freq in freqs:
        if not isinstance(freq,float) and not isinstance(freq,int):
            raise TypeError("The inputted frequency (" + repr(freq) +
                                ") must be a float or an int.")

    #Sort the freqs list.
    if isinstance(freqs,list):
        tmpSet = set(freqs)
        tmpList = list(tmpSet)
    else:
        tmpList = list(freqs)
    sortedFreqs = sorted(tmpList)

    #Initialize local variables.
    newFileLines = []
    prevLine = ""
    freqCounter = 1
    currentFreq = sortedFreqs[0]
    currentDiff = 1.e12

    #Open the inputted file.
    f = open(hitFileIn,
             "r")

    #Get the minimum frequency in the file by reading the first line.
    line = f.readline()
    lineBuf = line[3:]
    vals = lineBuf.split()
    minFileFreq = float(vals[0])

    #Check that there are not any inputted frequencies smaller than the
    #smallest frequency in the inputted hitran file.
    if sortedFreqs[0] < minFileFreq:
        raise ValueError("the inputted frequency (" + str(sortedFreqs[0]) +
                             ") is smaller than the smallest frequency (" +
                             str(minFileFreq) + ") in the HITRAN file (" +
                             hitFileIn + ".\n")

    #Rewind to the beginning of the file.
    f.seek(0)

    #Read through the file line by line.
    for line in f:

        #Get the line frequency from the file.
        lineBuf = line[3:]
        vals = lineBuf.split()
        fileFreq = float(vals[0])

        #Calculate the difference between the current frequency and the
        #frequency on the line of the inputted file.
        fDiff = abs(currentFreq - fileFreq)

        if fDiff <= currentDiff:

            #Update the difference.
            currentDiff = fDiff
        else:

            #Make sure that the difference from the first line of the file
            #does not exceed the default currentDiff value.
            if prevLine == "":
                raise ValueError("previous line is an empty string.\n")

            #Store the previous line.
            newFileLines.append(prevLine)

            #Iterate the frequency counter and check if all frequencies
            #have been found.
            freqCounter += 1
            if freqCounter > len(sortedFreqs):
                break

            #Set the current frequency to be the next frequency in the
            #inputted frequency list or else break if there are no
            #frequencies left in the list.
            lowBound = freqCounter - 1
            for i in range(lowBound,len(sortedFreqs)):
                currentFreq = sortedFreqs[i]
                prevDiff = abs(prevFileFreq - currentFreq)
                fDiff = abs(fileFreq - currentFreq)
                if prevDiff < fDiff:
                    newFileLines.append(prevLine)
                    freqCounter += 1
                else:
                    break
            currentDiff = fDiff

        #Store a copy of the first line.  This is needed because in order
        #to determine which line is closest to the desired frequency, we
        #must read one line past the line we ultimately want.
        prevLine = line
        prevFileFreq = fileFreq

    #Close the inputted file.
    f.close()

    #Make sure that all inputted lines were found.
    if len(newFileLines) != len(sortedFreqs):
        raise ValueError("found " + str(len(newFileLines)) + " out of " +
                             str(len(sortedFreqs)) + " inputted frequencies" +
                             " in the hitran file (" + hitFileIn + ").  " +
                             "One or more of the inputted frequencies is" +
                             " larger than the maximum frequency in the" +
                             " hitran file.")

    #Remove duplicates from newFileLines.
    tmpSet = set(newFileLines)
    tmpList = list(tmpSet)
    tmpList = sorted(tmpList)

    #Open a new file and write out the necessary lines.
    newHitFile = hitFileIn + ".lines_specific"
    f = open(newHitFile,
             "w")
    for line in tmpList:
        f.write(line)
    f.close()

    return newHitFile

