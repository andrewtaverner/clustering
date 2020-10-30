'''
Created on Jun 10, 2015

New one

@author: PiTav
'''
import numpy
import matplotlib.pyplot as plt
from scipy import stats
import itertools # Maybe I should fix this so all the groupby calls are either just "groupby" or "itertools.groupby"
from Bio import Seq, PDB, SeqIO, BiopythonWarning
from itertools import groupby
import random
import math
import time
import gzip
import warnings

translationTable = {'CTT': 'L', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 
                    'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 
                    'AGC': 'S', 'AAG': 'K', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 
                    'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 
                    'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CAA': 'Q', 
                    'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 
                    'CAG': 'Q', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 
                    'TGC': 'C', 'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TGG': 'W', 
                    'GGC': 'G', 'TAC': 'Y', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 
                    'TTG': 'L', 'TCC': 'S', 'GAA': 'E', 'TAA': '*', 'GCA': 'A', 
                    'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 
                    'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'TGA': '*', 'GAC': 'D', 
                    'CGT': 'R', 'TCA': 'S', 'ATG': 'M', 'CGC': 'R'}

def filterAlignment(firstSeq,secondSeq,thirdSeq,lenMultipleOfThree = True, minAlnLen = 300, 
                    eachGapMultipleOfThree = True, totalGapLengthMultipleOfThree = True,
                    gapContentThresh = .2, DNContentThresh = .2, DSContentThresh = 1, 
                    prematureStopCodon = True):
    # firstSeq, secondSeq, and thirdSeq should be Seq objects, not Seq records or strings
    # assigning 0 to minAlnLen means this filter is turned off (makes sense since False is normally 0)
    # First gap length filter is stronger than second, if first true, second will be (line 2)
    # Percent thresholds must be set to 1 to turn them off
    # DNContentThresh looks at all pairwise combinations of species to filter by DN amount
    # Same length filter is manditory
    seqLength = len(firstSeq)
    if len(firstSeq) != len(secondSeq) or len(firstSeq) != len(thirdSeq) or len(secondSeq) != len(thirdSeq):
        return (-1, "Different lengths")
    if lenMultipleOfThree and seqLength%3 > 0:
        return (-1, "Seq length not multiple of three")
    if minAlnLen > 0 and seqLength < minAlnLen:
        return (-1, "Seq length too short")
    if (eachGapMultipleOfThree and 
        (len([len(list(g)) for k, g in groupby(firstSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0 or
        len([len(list(g)) for k, g in groupby(secondSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0 or
        len([len(list(g)) for k, g in groupby(thirdSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0)):
        return (-1, "Each gap not multiple of three")
    # Only need to do this if the "eachGapMultipleOfThree" test is not used:
    if not eachGapMultipleOfThree and totalGapLengthMultipleOfThree and (firstSeq.count("-")%3 > 0 or 
        secondSeq.count("-")%3 > 0 or thirdSeq.count("-")%3 > 0):
        return (-1, "Total gap length does not sum to three")
    # gapContentThresh filter:
    if (firstSeq.count("-") > gapContentThresh * seqLength or secondSeq.count("-") > gapContentThresh * seqLength or
        thirdSeq.count("-") > gapContentThresh * seqLength):
        return (-1, "Greater than 20% gap")
    result = variantsNotPolarized(firstSeq,secondSeq)
    firstSecondDNAmount = result['DN Count']
    firstSecondDSAmount = result['DS Count']
    firstSecondStopCount = result['Stop Codon Count']
    result = variantsNotPolarized(firstSeq,thirdSeq)
    firstThirdDNAmount = result['DN Count']
    firstThirdDSAmount = result['DS Count']
    firstThirdStopCount = result['Stop Codon Count']
    result = variantsNotPolarized(secondSeq,thirdSeq)
    secondThirdDNAmount = result['DN Count']
    secondThirdDSAmount = result['DS Count']
    secondThirdStopCount = result['Stop Codon Count']
    # DNContentThresh filter
    if (firstSecondDNAmount > DNContentThresh * seqLength or firstThirdDNAmount > DNContentThresh * seqLength or 
        secondThirdDNAmount > DNContentThresh * seqLength):
        return (-1, "Too many amino acid substitutions")
    # DSContentThresh filter
    if (firstSecondDSAmount > DSContentThresh * seqLength or firstThirdDSAmount > DSContentThresh * seqLength or 
        secondThirdDSAmount > DSContentThresh * seqLength):
        return (-1, "Too many synonymous substitutions")
    if prematureStopCodon and (firstSecondStopCount > 0 or firstThirdStopCount > 0 or secondThirdStopCount > 0):
        return (-1, "Alignment contains stop codons")
    return (1, "Alignment passes filters")
    
def clusteringSameMutation(mutation_positions,clustering_calc_length = 500):
    mutation_positions = numpy.array(mutation_positions)//3
    clusteringDistances = mutation_positions - mutation_positions[:,None]
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances > 0,
        clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    return clusteringCount

def clusteringSameMutationIntrons(mutation_positions,intron_positions,geneLength,num_simulations = 100,clustering_calc_length = 500,seed = None):
    # NB: mutation_positions are in 0-indexed nucl space, and are all %3 = 0, intron_positions are in 0-indexed AA space
    # ALSO, geneLength is already in AA length
    # NB 2: This also calculates/simulates the null distribution of clustering... I'm hoping to "cheat" a bit and speed things up by reusing the intron mask
    # Turn mutation position into AA space coordinates:
    mutation_positions = numpy.array(mutation_positions)//3
    mutation_positions = numpy.setdiff1d(mutation_positions,intron_positions) # We don't want the mutations that are in a codon with an intron, so why not just remove them
    ## For testing:
    #originalIntronPos = intron_positions
    # First form the intron filtering mask:
    tempLogicalMatrix = numpy.zeros(shape = (len(mutation_positions), len(mutation_positions)),dtype = bool)
    mutation1 = numpy.tile(mutation_positions,(len(mutation_positions),1))
    mutation2 = numpy.transpose(mutation1)
    for currIntronPosition in intron_positions:
        tempLogicalMatrix = numpy.logical_or(tempLogicalMatrix,numpy.logical_and(currIntronPosition < mutation1, currIntronPosition > mutation2))
    # Then calculate clustering and filter for intron and length:
    clusteringDistances = mutation_positions - mutation_positions[:,None]
    clusteringDistances = clusteringDistances[tempLogicalMatrix]
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances > 0,
                                                                clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    # Now simulate:
    if seed != None: random.seed(seed)
    mutationsPerIntron = numpy.zeros((len(intron_positions)+1),dtype = int)
    #simulationHolder = numpy.zeros(len(clusteringCount)) # commented for testing
    simulationHolder = numpy.zeros((num_simulations,len(clusteringCount)))
    # Gene length should already be length in AA!
    intron_positions = numpy.append(numpy.append(-1,intron_positions),geneLength+1)
    for i in range(len(intron_positions)-1):
        mutationsPerIntron[i] = sum(numpy.logical_and(mutation_positions > intron_positions[i],mutation_positions < intron_positions[i+1]))
    cumulativeMutationCount = numpy.append(0,numpy.cumsum(mutationsPerIntron))
    for i in range(num_simulations):
        simMut = numpy.zeros(sum(mutationsPerIntron),dtype = int)
        for j in range(len(mutationsPerIntron)):
            simMut[cumulativeMutationCount[j]:cumulativeMutationCount[j+1]] = random.sample(range(int(math.floor(intron_positions[j]+1)),int(math.ceil(intron_positions[j+1]))),mutationsPerIntron[j])
        simMut = numpy.sort(simMut)
        ## Testing to see whether cheating gives the same result:
        #testingHolder = numpy.zeros(shape = (len(mutation_positions), len(mutation_positions)),dtype = bool)
        #simMut1Array = numpy.tile(simMut,(len(simMut),1))
        #simMut2Array = numpy.transpose(simMut1Array)
        #for currIntronPosition in originalIntronPos:
        #    testingHolder = numpy.logical_or(testingHolder,numpy.logical_and(currIntronPosition < simMut1Array, currIntronPosition > simMut2Array))
        #print(numpy.all(tempLogicalMatrix == testingHolder))
        simulationDistances = simMut - simMut[:,None]
        simulationDistances = simulationDistances[tempLogicalMatrix]
        simulationDistances = simulationDistances[numpy.logical_and(simulationDistances > 0, simulationDistances <= clustering_calc_length)]
        simulationHolder[i,:] = numpy.histogram(simulationDistances,bins = range(clustering_calc_length + 2))[0]
        #simulationHolder += numpy.histogram(simulationDistances,bins = range(clustering_calc_length + 2))[0] # commented for testing
    #simulationHolder = simulationHolder/num_simulations # commented for testing
    return {'data':clusteringCount,'norm':simulationHolder}

def clusteringDifferentMutation(mutation_positions_1,mutation_positions_2,clustering_calc_length = 500):
    mutation_positions_1 = numpy.array(mutation_positions_1)//3
    mutation_positions_2 = numpy.array(mutation_positions_2)//3
    clusteringDistances = numpy.abs(mutation_positions_1 - mutation_positions_2[:,None])
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances > 0,
        clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    return clusteringCount

def clusteringDifferentMutationIntrons(mutation_positions_1,mutation_positions_2,intron_positions,geneLength,num_simulations = 100,clustering_calc_length = 500,seed = None):
    ## For testing:
    #originalIntronPos = intron_positions
    # Turn mutation position into AA space coordinates:
    mutation_positions_1 = numpy.array(mutation_positions_1)//3
    mutation_positions_2 = numpy.array(mutation_positions_2)//3
    mutation_positions_1 = numpy.setdiff1d(mutation_positions_1,intron_positions)
    mutation_positions_2 = numpy.setdiff1d(mutation_positions_2,intron_positions)
    # First form the intron filtering mask:
    tempLogicalMatrix = numpy.zeros(shape = (len(mutation_positions_1),len(mutation_positions_2)),dtype = bool)
    mutation1 = numpy.tile(mutation_positions_1[:,None],(1,len(mutation_positions_2)))
    mutation2 = numpy.tile(mutation_positions_2,(len(mutation_positions_1),1))
    for currIntronPosition in intron_positions:
        tempLogicalMatrix = numpy.logical_or(tempLogicalMatrix,numpy.logical_or(numpy.logical_and(currIntronPosition < mutation1, currIntronPosition > mutation2),
                                                                                numpy.logical_and(currIntronPosition > mutation1, currIntronPosition < mutation2)
                                                                                ))
    clusteringDistances = numpy.abs(mutation_positions_1[:,None] - mutation_positions_2)
    clusteringDistances = clusteringDistances[tempLogicalMatrix]
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances != 0, clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    # Now simulate:
    if seed != None: random.seed(seed)
    mutation1PerIntron = numpy.zeros(len(intron_positions) + 1, dtype = int)
    mutation2PerIntron = numpy.zeros(len(intron_positions) + 1, dtype = int)
    simulationHolder = numpy.zeros(len(clusteringCount))
    intron_positions = numpy.append(numpy.append(-1,intron_positions),geneLength+1)
    for i in range(len(intron_positions) - 1):
        mutation1PerIntron[i] = sum(numpy.logical_and(mutation_positions_1 > intron_positions[i],mutation_positions_1 < intron_positions[i+1]))
        mutation2PerIntron[i] = sum(numpy.logical_and(mutation_positions_2 > intron_positions[i],mutation_positions_2 < intron_positions[i+1]))
    mutationsPerIntron = mutation1PerIntron + mutation2PerIntron
    cumulativeMutation1Count = numpy.append(0,numpy.cumsum(mutation1PerIntron))
    cumulativeMutation2Count = numpy.append(0,numpy.cumsum(mutation2PerIntron))
    for i in range(num_simulations):
        simMut1 = numpy.zeros(sum(mutation1PerIntron),dtype = int)
        simMut2 = numpy.zeros(sum(mutation2PerIntron),dtype = int)
        for j in range(len(mutationsPerIntron)):
            tempHolder = random.sample(range(int(math.floor(intron_positions[j]+1)),int(round(intron_positions[j+1]))),mutationsPerIntron[j])
            simMut1[cumulativeMutation1Count[j]:cumulativeMutation1Count[j+1]] = random.sample(tempHolder,mutation1PerIntron[j])
            simMut2[cumulativeMutation2Count[j]:cumulativeMutation2Count[j+1]] = numpy.setdiff1d(tempHolder,simMut1[cumulativeMutation1Count[j]:cumulativeMutation1Count[j+1]])
        ## testing to see whether cheating gives the same results:
        #testingHolder = numpy.zeros(shape = (len(mutation_positions_1),len(mutation_positions_2)),dtype = bool)
        #simMut1Array = numpy.tile(simMut1[:,None],(1,len(simMut2)))
        #simMut2Array = numpy.tile(simMut2,(len(simMut1),1))
        #for currIntronPosition in originalIntronPos:
        #    testingHolder = numpy.logical_or(testingHolder,numpy.logical_or(numpy.logical_and(currIntronPosition < simMut1Array, currIntronPosition > simMut2Array),
        #                                                                    numpy.logical_and(currIntronPosition > simMut1Array, currIntronPosition < simMut2Array)
        #                                                                    ))
        #print(numpy.all(tempLogicalMatrix == testingHolder))
        simulationDistances = numpy.abs(simMut1[:,None] - simMut2)
        simulationDistances = simulationDistances[tempLogicalMatrix]
        simulationDistances = simulationDistances[numpy.logical_and(simulationDistances != 0, simulationDistances <= clustering_calc_length)]
        simulationHolder += numpy.histogram(simulationDistances,bins = range(clustering_calc_length + 2))[0]
    simulationHolder = simulationHolder/num_simulations
    return {'data':clusteringCount,'norm':simulationHolder}

def analyticalNormSameMut(numMutations, geneLength, clustering_calc_length = 500):
    codonRange = numpy.arange(clustering_calc_length + 1)
    numMutations = float(numMutations)
    tempHolder = (numMutations*(numMutations-1))/(geneLength-1) - (numMutations*(numMutations-1))/(geneLength*(geneLength - 1))*codonRange
    #tempHolder = -((float(numMutations)**2-numMutations)/geneLength**2)*codonRange + (float(numMutations)**2-numMutations)/geneLength
    tempHolder[tempHolder < 0] = 0
    return tempHolder

def analyticalNormDiffMut(numMutations1, numMutations2, geneLength, clustering_calc_length = 500):
    codonRange = numpy.arange(clustering_calc_length + 1)
    numMutations1 = float(numMutations1)
    numMutations2 = float(numMutations2)
    tempHolder = ((2*numMutations1*numMutations2)/(geneLength-1)) - ((2*numMutations1*numMutations2)/(geneLength*(geneLength-1)))*codonRange
    #tempHolder = -(2*float(numMutations1*numMutations2)/geneLength**2)*codonRange + 2*float(numMutations1*numMutations2)/geneLength
    #tempHolder = -(2*(float(numMutations1)*numMutations2)/geneLength**2)*codonRange + 2*(float(numMutations1)*numMutations2)/geneLength
    tempHolder[tempHolder < 0] = 0
    return tempHolder

def propertyClusteringNew(firstSeq,secondSeq,outSeq,speciesOneDN,speciesTwoDN,lookupTable):
    if '-' in firstSeq:
        firstSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in firstSeq]))
    if '-' in secondSeq:
        secondSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in secondSeq]))
    if '-' in outSeq:
        outSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in outSeq]))
    firstSeqTrans = firstSeq.translate()
    secondSeqTrans = secondSeq.translate()
    outSeqTrans = outSeq.translate()
    speciesOneInc = [x for x in speciesOneDN if (lookupTable[firstSeqTrans[x//3]] - lookupTable[outSeqTrans[x//3]]) > 0]
    speciesOneDec = [x for x in speciesOneDN if (lookupTable[firstSeqTrans[x//3]] - lookupTable[outSeqTrans[x//3]]) < 0]
    speciesTwoInc = [x for x in speciesTwoDN if (lookupTable[secondSeqTrans[x//3]] - lookupTable[outSeqTrans[x//3]]) > 0]
    speciesTwoDec = [x for x in speciesTwoDN if (lookupTable[secondSeqTrans[x//3]] - lookupTable[outSeqTrans[x//3]]) < 0]
    speciesOneRein = clusteringSameMutation(speciesOneInc) + clusteringSameMutation(speciesOneDec)
    speciesOneComp = clusteringDifferentMutation(speciesOneInc,speciesOneDec)
    speciesTwoRein = clusteringSameMutation(speciesTwoInc) + clusteringSameMutation(speciesTwoDec)
    speciesTwoComp = clusteringDifferentMutation(speciesTwoInc,speciesTwoDec)
    speciesCrossRein = (clusteringDifferentMutation(speciesOneInc,speciesTwoInc) + 
                        clusteringDifferentMutation(speciesOneDec,speciesTwoDec))
    speciesCrossComp = (clusteringDifferentMutation(speciesOneInc,speciesTwoDec) + 
                        clusteringDifferentMutation(speciesOneDec,speciesTwoInc))
    # The first four things output by this function are more for debugging
    return {'Species One Dec':speciesOneDec,'Species One Inc':speciesOneInc,
            'Species Two Dec':speciesTwoDec,'Species Two Inc':speciesTwoInc,
            'Species One Comp':speciesOneComp,'Species One Rein':speciesOneRein,
            'Species Two Comp':speciesTwoComp,'Species Two Rein':speciesTwoRein,
            'Cross Species Comp':speciesCrossComp,'Cross Species Rein':speciesCrossRein}

def clusteringMetric(variant,normalization,clusteringLowerLimit=1,
                     clusteringUpperLimit=100,normLowerLimit=300,
                     normUpperLimit=350):
    # Restricting the length of the vector should elimintate numpy divide by zero warnings
    # if the selected range to calulate the "clustering metric" is ok
    clustering = variant[:normUpperLimit].astype(float)/normalization[:normUpperLimit]

    return sum(clustering[clusteringLowerLimit:clusteringUpperLimit] - numpy.mean(clustering[normLowerLimit:normUpperLimit]))

def pdbMatchHelper(startingPositionGuess, endingPositionGuess, currPeptide, proteinRef):
    # If the peptide sequence has a gap, skip this peptide:
    if len(currPeptide) != endingPositionGuess - startingPositionGuess:
        return -1
    numMatches = numpy.sum([x == y for x,y in zip(list(currPeptide.get_sequence()), list(proteinRef[startingPositionGuess:endingPositionGuess]))])
    peptideLength = len(currPeptide)
    if numMatches/peptideLength < 0.9:
        print("Finding 'alignment'")
        peptideSeq = currPeptide.get_sequence()
        peptideFragment = peptideSeq[0:10] if peptideLength > 10 else peptideSeq
        peptideFragmentLen = len(peptideFragment)
        tempList = []
        for i in range(len(proteinRef) - peptideFragmentLen + 1):
            tempList.append(numpy.sum([x[0] == x[1] for x in zip(proteinRef[i:i+peptideFragmentLen],peptideFragment)]))
        tempList = numpy.array(tempList)
        possibleStartingPositions = numpy.where(tempList >= 9)[0]
        if len(possibleStartingPositions) == 0:
            return -1
        if len(possibleStartingPositions) > 1 and peptideLength <= 10: 
            print("Multiple matches for short peptide, aborting")
            return -1
        candidateMatchVector = []
        for currStartingPos in possibleStartingPositions:
            candidateMatchVector.append([numpy.sum([x[0] == x[1] for x in zip(proteinRef[currStartingPos:currStartingPos + peptideLength],peptideSeq)])][0])
        candidateMatchVector = numpy.array(candidateMatchVector)
        if len(candidateMatchVector) > 1 and sum(candidateMatchVector == max(candidateMatchVector)) > 1:
            # This is so unlikely to happen, it might have a bug and I'd never know
            print("Multiple matches for long peptide, aborting")
            return -1
        startingPosition = possibleStartingPositions[candidateMatchVector == max(candidateMatchVector)][0]
        endingPosition = startingPosition + peptideLength
        if max(candidateMatchVector)/peptideLength < 0.9:
            print("Unable to find a good match, aborting")
            return -1
        return (startingPosition,endingPosition,max(candidateMatchVector))
    else:
        refProteinFragment = proteinRef.rstrip('*')[startingPositionGuess:endingPositionGuess]
        return (startingPositionGuess,endingPositionGuess,numMatches)

def clusteringPDB(sequenceHolderFileName, geneName, speciesOneName, speciesTwoName, speciesOutName, PDBList, pdbPath,sequenceHolder = None):
    warnings.simplefilter('ignore',BiopythonWarning)
    
    startTime = time.time()
    
    print("Processing Gene: " + geneName)
    
    parser = PDB.MMCIFParser()
    ppb = PDB.PPBuilder()
    
    resultHolder = dict()
    resultHolder['Gene name'] = geneName
    try:    
        if sequenceHolder == None:
            sequenceHolder = SeqIO.to_dict(SeqIO.parse(open(sequenceHolderFileName),'fasta'))
        
        currSeqName = dict([(x,y) for y in sequenceHolder.keys() for x in [speciesOneName,speciesTwoName,speciesOutName] if x in y])
        if len(currSeqName) < 3:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " beacuse missing one or more of the three species"
            return resultHolder
        
        speciesOneName = currSeqName[speciesOneName]
        speciesTwoName = currSeqName[speciesTwoName]
        speciesOutName = currSeqName[speciesOutName]
        
        try:
            firstSeq = str(sequenceHolder[speciesOneName].seq)
            secondSeq = str(sequenceHolder[speciesTwoName].seq)
            outSeq = str(sequenceHolder[speciesOutName].seq)
        except:
            firstSeq = sequenceHolder[speciesOneName]
            secondSeq = sequenceHolder[speciesTwoName]
            outSeq = sequenceHolder[speciesOutName]
        
        if len(firstSeq) != len(secondSeq) or \
            len(firstSeq) != len(outSeq):
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because lengths don't match"
            return resultHolder
        if len(firstSeq) < 300:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because sequence too short"
            return resultHolder
        if len(firstSeq) % 3 > 0:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because sequence length not multiple of three"
            return resultHolder
        if firstSeq.count('-') > .2 * len(firstSeq) or secondSeq.count('-') > .2 * len(secondSeq) or outSeq.count('-') > len(outSeq):
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because sequence contained too many gaps"
            return resultHolder
        if len([len(list(g)) for k, g in itertools.groupby(firstSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species one gaps not multiple of three"
            return resultHolder
        if len([len(list(g)) for k, g in itertools.groupby(secondSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species two gaps not multiple of three"
            return resultHolder
        if len([len(list(g)) for k, g in itertools.groupby(outSeq) if (k == '-' and len(list(g))%3 != 0)]) > 0:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species out gaps not multiple of three"
            return resultHolder
        if '*' in Seq.Seq(''.join([x if x != '-' else 'N' for x in firstSeq])).translate()[:-1]:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species one contains a premature stop codon"
            return resultHolder
        if '*' in Seq.Seq(''.join([x if x != '-' else 'N' for x in secondSeq])).translate()[:-1]:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species two contains a premature stop codon"
            return resultHolder
        if '*' in Seq.Seq(''.join([x if x != '-' else 'N' for x in outSeq])).translate()[:-1]:
            resultHolder['ReasonForFailure'] = "Skipped " + geneName + " because species out contains a premature stop codon"
            return resultHolder
        
        firstSeqProtein = Seq.Seq(''.join([x if x != "-" else "N" for x in firstSeq])).translate()
        
        # PDB info:
        lengthOfMatchDict = dict()
        locationOfCarbonDict = dict()
        lengthOfPDBDict = dict()
        startingPositionDict = dict()
        endingPositionDict = dict()
        peptideIndexDict = dict()
        
        for currPDB in PDBList:
            for currChainID in PDBList[currPDB]:
                lengthOfMatchDict[(currPDB,currChainID)] = 0
                lengthOfPDBDict[(currPDB,currChainID)] = 0
                locationOfCarbonDict[(currPDB,currChainID)] = []
                peptideIndexDict[(currPDB,currChainID)] = []
                startingPositionDict[(currPDB,currChainID)] = []
                endingPositionDict[(currPDB,currChainID)] = []
        
        for pdbInd, currPDB in enumerate(PDBList):
            print('    Processing PDB ' + str(pdbInd + 1) + ' of ' + str(len(PDBList)))
            try:
                fileReader = gzip.open(pdbPath + currPDB + '.cif.gz','rt')
                structure = parser.get_structure(currPDB,fileReader)
                peptides = ppb.build_peptides(structure)
            except Exception:
                print('Unable to read PDB file')
                continue
            if len(peptides) == 0:
                continue
            for currChainID in PDBList[currPDB]:
                # Get the peptide index within the PDB file:
                peptideIndex = [tempInd for tempInd,x in enumerate(peptides) if x[0].full_id[2] == currChainID]
                if len(peptideIndex) == 0:
                    continue
                offsetGuess = 0
                for currPeptideIndex in peptideIndex:
                    startingPosition = peptides[currPeptideIndex][0].full_id[3][1] - 1 + offsetGuess
                    endingPosition = peptides[currPeptideIndex][-1].full_id[3][1] + offsetGuess
                    lengthOfPDBDict[(currPDB,currChainID)] += len(peptides[currPeptideIndex])    
                    result = pdbMatchHelper(startingPosition,endingPosition,peptides[currPeptideIndex],firstSeqProtein)
                    if result != -1:
                        startingPositionDict[(currPDB,currChainID)].append(result[0])
                        endingPositionDict[(currPDB,currChainID)].append(result[1])
                        lengthOfMatchDict[(currPDB,currChainID)] += result[2]
                        peptideIndexDict[(currPDB,currChainID)].append(currPeptideIndex)
                        locationOfCarbonDict[(currPDB,currChainID)].append(['CA' in x for x,y in zip(peptides[currPeptideIndex],firstSeqProtein[result[0]:result[1]])])
                        offsetGuess = result[0] - startingPosition
        
        pdbToRemove = [x for x in lengthOfMatchDict if lengthOfPDBDict[x] == 0 or lengthOfMatchDict[x]/lengthOfPDBDict[x] < .9]
        #print(pdbToRemove)
        
        for currKey in pdbToRemove:
            del lengthOfMatchDict[currKey] # we only need to remove this one because we'll get out the key corresponding to the longest value
        
        if len(lengthOfMatchDict) == 0:
            resultHolder['ReasonForFailure'] = 'Skipped ' + geneName + " because PDB sequence didn't match sequence"
            return resultHolder
        
        maxLength = max(lengthOfMatchDict.values())
        pdbToUse = [x for x in lengthOfMatchDict if lengthOfMatchDict[x] == maxLength][0]
        
        fileReader = gzip.open(pdbPath + pdbToUse[0] + '.cif.gz','rt')
        structure = parser.get_structure(pdbToUse[0],fileReader)
        peptides = ppb.build_peptides(structure)
        
        editedFirstSeq = []
        editedSecondSeq = []
        editedOutSeq = []
        ca_atomList = []
        proteinSeq = []
        for ind in range(len(peptideIndexDict[pdbToUse])):
            startingPosition = startingPositionDict[pdbToUse][ind]
            endingPosition = endingPositionDict[pdbToUse][ind]
            positionsToMask = numpy.where(numpy.array(locationOfCarbonDict[pdbToUse][ind]) == False)[0]
            tempHolderOne = firstSeq[startingPosition*3:endingPosition*3]
            tempHolderTwo = secondSeq[startingPosition*3:endingPosition*3]
            tempHolderThree = outSeq[startingPosition*3:endingPosition*3]
            if len(positionsToMask) > 0:
                tempHolderOne = bytearray(tempHolderOne,'utf8')
                tempHolderTwo = bytearray(tempHolderTwo,'utf8')
                tempHolderThree = bytearray(tempHolderThree,'utf8')
                # I just treat M like it's a mask; really doesn't matter, as long as all three match
                for seqInd in positionsToMask:
                    tempHolderOne[seqInd*3:(seqInd*3 + 3)] = bytearray(b'ATG')
                    tempHolderTwo[seqInd*3:(seqInd*3 + 3)] = bytearray(b'ATG')
                    tempHolderThree[seqInd*3:(seqInd*3 + 3)] = bytearray(b'ATG')
                tempHolderOne = tempHolderOne.decode()
                tempHolderTwo = tempHolderTwo.decode()
                tempHolderThree = tempHolderThree.decode()
            editedFirstSeq.append(tempHolderOne)
            editedSecondSeq.append(tempHolderTwo)
            editedOutSeq.append(tempHolderThree)
            [ca_atomList.append(x['CA']) if y else ca_atomList.append(None) for x,y in zip(peptides[peptideIndexDict[pdbToUse][ind]],locationOfCarbonDict[pdbToUse][ind])]
            proteinSeq.append(''.join([PDB.Polypeptide.three_to_one(x.get_resname()) if y else 'M' for x,y in zip(peptides[peptideIndexDict[pdbToUse][ind]],locationOfCarbonDict[pdbToUse][ind])]))
            assert(len(''.join(editedFirstSeq))/3 == len(ca_atomList))
            #[test.append(x.resname) if y else test.append(None) for x,y in zip(peptides[currPeptideIndex],locationOfCarbonDict[pdbToUse][ind])]
            if ind + 1 < len(peptideIndexDict[pdbToUse]):
                startingPosition = endingPosition
                endingPosition = startingPositionDict[pdbToUse][ind+1]
                editedFirstSeq.append('ATG'*(endingPosition - startingPosition))
                editedSecondSeq.append('ATG'*(endingPosition - startingPosition))
                editedOutSeq.append('ATG'*(endingPosition - startingPosition))
                proteinSeq.append('M'*(endingPosition - startingPosition))
                [ca_atomList.append(None) for x in range(endingPosition - startingPosition)]
                #[test.append(None) for x in range(endingPosition - startingPosition)]
        
        editedFirstSeq = ''.join(editedFirstSeq)
        editedSecondSeq = ''.join(editedSecondSeq)
        editedOutSeq = ''.join(editedOutSeq)
        proteinSeq = ''.join(proteinSeq)
        
        if len(editedFirstSeq)/3 != len(ca_atomList):
            resultHolder['ReasonForFailure'] = "Skipping " + geneName + " because Ca location vector is different length than expected"
            return resultHolder
        
        numMatches = sum([x[0] == x[1] for x in zip(proteinSeq,str(Seq.Seq(''.join([x if x != '-' else 'N' for x in editedFirstSeq])).translate()))])
        if numMatches/float(len(proteinSeq)) < 0.95:
            resultHolder['ReasonForFailure'] = "Skipping " + geneName + " unable to process PDB to match gene sequence"
            return resultHolder
        
        results = variantsNotPolarized(editedFirstSeq,editedSecondSeq)
        if results['DS Count'] == 0 and results['DN Count'] == 0:
            resultHolder['ReasonForFailure'] = 'Skipping ' + geneName + ' trimmed protein does not contain any substitutions'
            return resultHolder
        
        DNList = numpy.array(results['DN List'])//3
        DSList = numpy.array(results['DS List'])//3
        
        results = variantsByParsimony(editedFirstSeq,editedSecondSeq,editedOutSeq)
        
        DNListSpeciesOne = numpy.array(results['Species One DN'])//3
        DNListSpeciesTwo = numpy.array(results['Species Two DN'])//3
        
        DNNormTempHolder = numpy.array([])
        DSNormTempHolder = numpy.array([])
        DNNormSpeciesOneTempHolder = numpy.array([])
        DNNormSpeciesTwoTempHolder = numpy.array([])
        DNDNSpeciesOneTempHolder = numpy.array([])
        DNDNSpeciesTwoTempHolder = numpy.array([])
        DNDNBetweenTempHolder = numpy.array([])
        DNDNTempHolder = numpy.array([])
        DSDSTempHolder = numpy.array([])
        DNDSTempHolder = numpy.array([])
        
        try:
            ca_atomList = numpy.array(ca_atomList)
            # DNDN species one
            if len(DNListSpeciesOne) > 1:
                DNDNSpeciesOneTempHolder = ca_atomList[DNListSpeciesOne] - ca_atomList[DNListSpeciesOne,None]
                DNDNSpeciesOneTempHolder = DNDNSpeciesOneTempHolder[numpy.triu_indices(len(DNListSpeciesOne),1)]
            # DNDN species two
            if len(DNListSpeciesTwo) > 1:
                DNDNSpeciesTwoTempHolder = ca_atomList[DNListSpeciesTwo] - ca_atomList[DNListSpeciesTwo,None]
                DNDNSpeciesTwoTempHolder = DNDNSpeciesTwoTempHolder[numpy.triu_indices(len(DNListSpeciesTwo),1)]
            # DNDN between species
            if len(DNListSpeciesOne) > 0 and len(DNListSpeciesTwo) > 0:
                DNDNBetweenTempHolder = ca_atomList[DNListSpeciesOne] - ca_atomList[DNListSpeciesTwo,None]
                DNDNBetweenTempHolder = DNDNBetweenTempHolder.flatten()
            # DNDN
            if len(DNList) > 1:
                DNDNTempHolder = ca_atomList[DNList] - ca_atomList[DNList,None]
                DNDNTempHolder = DNDNTempHolder[numpy.triu_indices(len(DNList),1)]
            # DSDS
            if len(DSList) > 1:
                DSDSTempHolder = ca_atomList[DSList] - ca_atomList[DSList,None]
                DSDSTempHolder = DSDSTempHolder[numpy.triu_indices(len(DSList),1)]
            # DNDS
            if len(DNList) > 0 and len(DSList) > 0:
                DNDSTempHolder = ca_atomList[DNList] - ca_atomList[DSList,None]
                DNDSTempHolder = DNDSTempHolder.flatten()
            # Old normalization (which is not only stupid, it's time consuming to calculate)
            for i in DNListSpeciesOne:
                tempHolder = ca_atomList[(i+1):]
                tempHolder = tempHolder[tempHolder != None]
                DNNormSpeciesOneTempHolder = numpy.append(DNNormSpeciesOneTempHolder,numpy.array(ca_atomList[i]) - tempHolder)
            for i in DNListSpeciesTwo:
                tempHolder = ca_atomList[(i+1):]
                tempHolder = tempHolder[tempHolder != None]
                DNNormSpeciesTwoTempHolder = numpy.append(DNNormSpeciesTwoTempHolder,numpy.array(ca_atomList[i]) - tempHolder)
            for i in DNList:
                tempHolder = ca_atomList[(i+1):]
                tempHolder = tempHolder[tempHolder != None]
                DNNormTempHolder = numpy.append(DNNormTempHolder,numpy.array(ca_atomList[i]) - tempHolder)
            for i in DSList:
                tempHolder = ca_atomList[(i+1):]
                tempHolder = tempHolder[tempHolder != None]
                DSNormTempHolder = numpy.append(DSNormTempHolder,numpy.array(ca_atomList[i]) - tempHolder)
            # New normalization attempt
            pseudoPrimaryPos = numpy.arange(len(ca_atomList))
            pseudoPrimaryPos = pseudoPrimaryPos[ca_atomList != None] # remove sites where we don't have information about where the alpha carbon is
            ca_atomList = ca_atomList[ca_atomList != None] # same as above comment, and we don't need this thing anymore at this point...
            tempNormHolder = ca_atomList - ca_atomList[:,None]
            primarySeqLen = pseudoPrimaryPos - pseudoPrimaryPos[:,None]
            tempNormHolder = tempNormHolder[numpy.triu_indices(len(ca_atomList),1)]
            primarySeqLen = primarySeqLen[numpy.triu_indices(len(ca_atomList),1)]
        except TypeError:
            resultHolder['ReasonForFailure'] = 'Something went wrong masking None from ca_atomList'
            return resultHolder
            
        
        timeElapsed = time.time() - startTime
        print("Processed " + geneName + " in " + str(timeElapsed) + " seconds")
        
        resultHolder = {'DNDN Lengths Species One':DNDNSpeciesOneTempHolder,'DNDN Lengths Species Two':DNDNSpeciesTwoTempHolder,
                        'DNDN Lengths Between':DNDNBetweenTempHolder,
                        'DNDN Lengths':DNDNTempHolder,'DNDS Lengths':DNDSTempHolder,'DSDS Lengths':DSDSTempHolder,
                        'DN Lengths Species One Norm':DNNormSpeciesOneTempHolder,
                        'DN Lengths Species Two Norm':DNNormSpeciesTwoTempHolder,
                        'New Norm':tempNormHolder,'DS Count':len(DSList),'DN Count':len(DNList),
                        'Species One DN Count':len(DNListSpeciesOne),'Species Two DN Count':len(DNListSpeciesTwo),
                        'PDBMatch':pdbToUse,'Time elapsed':timeElapsed,'Primary Seq Len':primarySeqLen, 'Gene name': geneName}
        return resultHolder
    except:
        resultHolder['ReasonForFailure'] = "God only knows what went wrong with this gene"
        return resultHolder

def readFasta(fileName,lengthQC = True):
    sequenceHolder = dict()
    fileHandle = open(fileName)
    currStr = fileHandle.readline()
    while currStr != '':
        if currStr[0] == '>':
            currOrg = currStr[1:]
            sequenceHolder[currOrg] = list()
            currStr = fileHandle.readline()
            while currStr != '' and currStr[0] != '>':
                sequenceHolder[currOrg].append(currStr.rstrip('\n'))
                currStr = fileHandle.readline()
    fileHandle.close()
    for currKey in sequenceHolder.keys():
        sequenceHolder[currKey] = ''.join(sequenceHolder[currKey])
    seqLength = len(sequenceHolder[sequenceHolder.keys()[0]])
    if lengthQC:
        for currKey in sequenceHolder.keys():
            if len(sequenceHolder[currKey]) != seqLength:
                return -1
    return sequenceHolder

def variantsNotPolarized(firstSeq, secondSeq):
    ambiguousGenotypeOrBreak = 0
    stopCodons = 0
    translationProblems = 0
    numDS = 0
    numDN = 0
    DSList = []
    DNList = []

    firstSeq = firstSeq.upper()
    secondSeq = secondSeq.upper()

    if not (len(firstSeq) == len(secondSeq)): return -1
    seqLength = len(firstSeq)
    for index in range(0,seqLength-3,3):
        codonStart = index
        codonEnd = index + 3
        if ('-' in firstSeq[codonStart:codonEnd] or 'N' in firstSeq[codonStart:codonEnd] or 'X' in firstSeq[codonStart:codonEnd] or
            '-' in secondSeq[codonStart:codonEnd] or 'N' in secondSeq[codonStart:codonEnd] or 'X' in secondSeq[codonStart:codonEnd]):
            ambiguousGenotypeOrBreak = ambiguousGenotypeOrBreak + 1
            continue
        # Count stop codons
        try:
            if (translationTable[str(firstSeq[codonStart:codonEnd])] == '*' or
                translationTable[str(secondSeq[codonStart:codonEnd])] == '*'):
                stopCodons = stopCodons + 1 
        except KeyError:
            pass

        # same
        if str(firstSeq[codonStart:codonEnd]) == str(secondSeq[codonStart:codonEnd]): continue
        try:
            if translationTable[str(firstSeq[codonStart:codonEnd])] == translationTable[str(secondSeq[codonStart:codonEnd])]:
                numDS = numDS + 1
                DSList.append(index)
            else:
                numDN = numDN + 1
                DNList.append(index)
        except KeyError:
            translationProblems = translationProblems + 1

    return {'Ambiguous Codon': ambiguousGenotypeOrBreak, 'Translation Problems':translationProblems,
            'DS Count':numDS,'DN Count':numDN,'DS List':DSList,'DN List':DNList,'Stop Codon Count':stopCodons,
            'Sequence Length':seqLength}

def gapClustering(firstSeq,secondSeq,windowLength):
    breakClusteringHolderRev = numpy.zeros([1,windowLength//3],dtype=int)
    opportunityRev = numpy.zeros([1,windowLength//3],dtype=int)
    breakClusteringHolderFor = numpy.zeros([1,windowLength//3],dtype=int)
    opportunityFor = numpy.zeros([1,windowLength//3],dtype=int)
    gapLocationTemp = numpy.diff([1 if x == '-' else 0 for x in "A" + secondSeq + "A"])
    gapStarts = numpy.where(gapLocationTemp == 1)[0]
    gapEnds = numpy.where(gapLocationTemp == -1)[0]
    # Corrects any probelems with codon starting time
    gapStarts = gapStarts//3*3
    gapEnds = gapEnds//3*3
    for currStart in gapStarts:
        if currStart == 0: continue
        tempSecondSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in secondSeq[max(currStart-windowLength,0):currStart]])).translate()
        tempFirstSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in firstSeq[max(currStart-windowLength,0):currStart]])).translate()
        tempHolder = numpy.array([1 if (x[0] != x[1] and x[0] != 'X' and x[1] != 'X') else 0 for x in zip(tempSecondSeq,tempFirstSeq)])
        tempHolderOpportunity = numpy.array([1 if (x[0] != 'X' and x[1] != 'X') else 0 for x in zip(tempSecondSeq,tempFirstSeq)])
        if len(tempHolder) < windowLength//3:
            tempHolder = numpy.concatenate((numpy.zeros([1,windowLength//3-len(tempHolder)],dtype=int)[0],tempHolder))
            tempHolderOpportunity = numpy.concatenate((numpy.zeros([1,windowLength//3-len(tempHolderOpportunity)],dtype=int)[0],tempHolderOpportunity))
        breakClusteringHolderRev += tempHolder
        opportunityRev += tempHolderOpportunity
    for currEnd in gapEnds:
        if currEnd == len(firstSeq): continue
        tempSecondSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in secondSeq[currEnd:currEnd+windowLength]])).translate()
        tempFirstSeq = Seq.Seq(''.join([x if x != "-" else "N" for x in firstSeq[currEnd:currEnd+windowLength]])).translate()
        tempHolder = numpy.array([1 if (x[0] != x[1] and x[0] != 'X' and x[1] != 'X') else 0 for x in zip(tempSecondSeq,tempFirstSeq)])
        tempHolderOpportunity = numpy.array([1 if (x[0] != 'X' and x[1] != 'X') else 0 for x in zip(tempSecondSeq,tempFirstSeq)])
        if len(tempHolder) < windowLength//3:
            tempHolder = numpy.concatenate((tempHolder,numpy.zeros([1,windowLength//3-len(tempHolder)],dtype=int)[0]))
            tempHolderOpportunity = numpy.concatenate((tempHolderOpportunity,numpy.zeros([1,windowLength//3-len(tempHolderOpportunity)],dtype=int)[0]))
        breakClusteringHolderFor += tempHolder
        opportunityFor += tempHolderOpportunity
    return {'revClustering':breakClusteringHolderRev,'opportunityRev':opportunityRev,
            'forClustering':breakClusteringHolderFor,'opportunityFor':opportunityFor}

def variantsByParsimony(firstSeq, secondSeq, outSeq):
    ambiguousGenotypeOrBreak = 0
    stopCodons = 0
    translationProblems = 0
    nonPolarizableDS = 0
    nonPolarizableDN = 0
    numSpeciesOneDS = 0
    numSpeciesOneDN = 0
    numSpeciesTwoDS = 0
    numSpeciesTwoDN = 0
    speciesOneDN = []
    speciesOneDS = []
    speciesTwoDN = []
    speciesTwoDS = []

    firstSeq = firstSeq.upper()
    secondSeq = secondSeq.upper()
    outSeq = outSeq.upper()

    if not (len(outSeq) == len(firstSeq) and len(firstSeq) == len(secondSeq)): return -1
    seqLength = len(outSeq)
    for index in range(0,seqLength-3,3):
        codonStart = index
        codonEnd = index + 3
        if ('-' in outSeq[codonStart:codonEnd] or 'N' in outSeq[codonStart:codonEnd] or 'X' in outSeq[codonStart:codonEnd] or
            '-' in firstSeq[codonStart:codonEnd] or 'N' in firstSeq[codonStart:codonEnd] or 'X' in firstSeq[codonStart:codonEnd] or
            '-' in secondSeq[codonStart:codonEnd] or 'N' in secondSeq[codonStart:codonEnd] or 'X' in secondSeq[codonStart:codonEnd]):
            ambiguousGenotypeOrBreak = ambiguousGenotypeOrBreak + 1
            continue

        # Count stop codons
        try:
            if (translationTable[str(firstSeq[codonStart:codonEnd])] == '*' or
                translationTable[str(secondSeq[codonStart:codonEnd])] == '*' or
                translationTable[str(outSeq[codonStart:codonEnd])] == '*'):
                stopCodons = stopCodons + 1 
        except KeyError:
            pass
        
        # There have been some changes to this logic. Reference the "finalizedForPub" version to see the original
        # exactly same
        if (str(firstSeq[codonStart:codonEnd]) == str(secondSeq[codonStart:codonEnd]) 
            and str(firstSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd])): continue
        # same (though differing from outgroup) 
        elif str(firstSeq[codonStart:codonEnd]) == str(secondSeq[codonStart:codonEnd]): continue
        # Difference, but not polarizable 
        elif (str(firstSeq[codonStart:codonEnd].translate()) != str(outSeq[codonStart:codonEnd].translate()) and
              str(secondSeq[codonStart:codonEnd].translate()) != str(outSeq[codonStart:codonEnd].translate())):
                nonPolarizableDN = nonPolarizableDN + 1
        # Difference in first seq
        elif str(secondSeq[codonStart:codonEnd].translate()) == str(outSeq[codonStart:codonEnd].translate()) and \
            str(firstSeq[codonStart:codonEnd].translate()) != str(outSeq[codonStart:codonEnd].translate()):
            speciesOneDN.append(index)
            numSpeciesOneDN = numSpeciesOneDN + 1            
        # Difference in secondSeq
        elif str(firstSeq[codonStart:codonEnd].translate()) == str(outSeq[codonStart:codonEnd].translate()) and \
            str(secondSeq[codonStart:codonEnd].translate()) != str(outSeq[codonStart:codonEnd].translate()):
            speciesTwoDN.append(index)
            numSpeciesTwoDN = numSpeciesTwoDN + 1
        # Now deal with DS:
        elif (str(firstSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd]) and
              str(secondSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd])):
                nonPolarizableDS = nonPolarizableDS + 1
        elif str(secondSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd]) and \
            str(firstSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd]):
            speciesOneDS.append(index)
            numSpeciesOneDS += 1
        elif str(firstSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd]) and \
            str(secondSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd]):
            speciesTwoDS.append(index)
            numSpeciesTwoDS += 1
        else:
            assert False
    
    return {'Ambiguous Codon':ambiguousGenotypeOrBreak,'Translation Problems':translationProblems,
            'Non-polarizable DS':nonPolarizableDS, 'Non-polarizable DN':nonPolarizableDN,
            'Species One DS Count':numSpeciesOneDS,'Species One DN Count':numSpeciesOneDN,
            'Species Two DS Count':numSpeciesTwoDS,'Species Two DN Count':numSpeciesTwoDN,
            'Species One DS':speciesOneDS,'Species One DN':speciesOneDN,
            'Species Two DS':speciesTwoDS,'Species Two DN':speciesTwoDN,
            'Sequence Length':seqLength,'Stop Codon Count':stopCodons}

# This is not used, but probably still works
def variantsByAncestor(firstSeq, secondSeq, ancSeq):
    ambiguousGenotypeOrBreak = 0
    translationProblems = 0
    stopCodons = 0
    numSpeciesOneDS = 0
    numSpeciesOneDN = 0
    numSpeciesTwoDS = 0
    numSpeciesTwoDN = 0
    speciesOneTotalCodon = 0
    speciesTwoTotalCodon = 0
    speciesOneDN = []
    speciesOneDS = []
    speciesTwoDN = []
    speciesTwoDS = []
    
    firstSeq = firstSeq.upper()
    secondSeq = secondSeq.upper()
    ancSeq = ancSeq.upper()

    if not (len(ancSeq) == len(firstSeq) and len(firstSeq) == len(secondSeq)): return -1
    seqLength = len(ancSeq)
    for index in range(0,seqLength-3,3):
        codonStart = index
        codonEnd = index + 3
        if ('-' in ancSeq[codonStart:codonEnd] or 'N' in ancSeq[codonStart:codonEnd] or 'X' in ancSeq[codonStart:codonEnd]):
            ambiguousGenotypeOrBreak = ambiguousGenotypeOrBreak + 1
            continue

        # Count stop codons, except if it occurs at the end of the sequence
        try:
            if (translationTable[firstSeq[codonStart:codonEnd]] == '*' or
                translationTable[secondSeq[codonStart:codonEnd]] == '*' or
                translationTable[ancSeq[codonStart:codonEnd]] == '*'):
                stopCodons = stopCodons + 1 
        except KeyError:
            pass

        if not ('-' in firstSeq[codonStart:codonEnd] or 'N' in firstSeq[codonStart:codonEnd] or 'X' in firstSeq[codonStart:codonEnd]):
            speciesOneTotalCodon = speciesOneTotalCodon + 1
            # Difference in first seq
            if firstSeq[codonStart:codonEnd] != ancSeq[codonStart:codonEnd]:
                try:
                    if translationTable[firstSeq[codonStart:codonEnd]] == translationTable[ancSeq[codonStart:codonEnd]]:
                        speciesOneDS.append(index)
                        numSpeciesOneDS = numSpeciesOneDS + 1
                    else:
                        speciesOneDN.append(index)
                        numSpeciesOneDN = numSpeciesOneDN + 1
                except KeyError:
                    translationProblems = translationProblems + 1
        if not ('-' in secondSeq[codonStart:codonEnd] or 'N' in secondSeq[codonStart:codonEnd] or 'X' in secondSeq[codonStart:codonEnd]):
            speciesTwoTotalCodon = speciesTwoTotalCodon + 1
            # Difference in secondSeq
            if secondSeq[codonStart:codonEnd] != ancSeq[codonStart:codonEnd]:
                try:
                    if translationTable[secondSeq[codonStart:codonEnd]] == translationTable[ancSeq[codonStart:codonEnd]]:
                        speciesTwoDS.append(index)
                        numSpeciesTwoDS = numSpeciesTwoDS + 1
                    else:
                        speciesTwoDN.append(index)
                        numSpeciesTwoDN = numSpeciesTwoDN + 1
                except KeyError:
                    translationProblems = translationProblems + 1
                
    return {'Ambiguous Codon':ambiguousGenotypeOrBreak,'Translation Problems':translationProblems,
            'Species One DS Count':numSpeciesOneDS,'Species One DN Count':numSpeciesOneDN,
            'Species Two DS Count':numSpeciesTwoDS,'Species Two DN Count':numSpeciesTwoDN,
            'Species One DS':speciesOneDS,'Species One DN':speciesOneDN,'Stop Codon Count':stopCodons,
            'Species Two DS':speciesTwoDS,'Species Two DN':speciesTwoDN,'Sequence Length':seqLength,
            'Species One Total Codon':speciesOneTotalCodon,'Species Two Total Codon':speciesTwoTotalCodon}

# not used (keeping for error checking)
def propertyVariants(firstSeq, secondSeq, outSeq, lookupTable):
    ambiguousGenotypeOrBreak = 0
    stopCodons = 0
    translationProblems = 0
    nonPolarizableDS = 0
    nonPolarizableDN = 0
    numSpeciesOneInc = 0
    numSpeciesOneDec = 0
    numSpeciesOneNeu = 0
    numSpeciesTwoInc = 0
    numSpeciesTwoDec = 0
    numSpeciesTwoNeu = 0
    speciesOneInc = []
    speciesOneDec = []
    speciesTwoInc = []
    speciesTwoDec = []

    firstSeq = firstSeq.upper()
    secondSeq = secondSeq.upper()
    outSeq = outSeq.upper()

    if not (len(outSeq) == len(firstSeq) and len(firstSeq) == len(secondSeq)): return -1
    seqLength = len(outSeq)
    for index in range(0,seqLength-3,3):
        codonStart = index
        codonEnd = index + 3
        if ('-' in outSeq[codonStart:codonEnd] or 'N' in outSeq[codonStart:codonEnd] or 'X' in outSeq[codonStart:codonEnd] or
            '-' in firstSeq[codonStart:codonEnd] or 'N' in firstSeq[codonStart:codonEnd] or 'X' in firstSeq[codonStart:codonEnd] or
            '-' in secondSeq[codonStart:codonEnd] or 'N' in secondSeq[codonStart:codonEnd] or 'X' in secondSeq[codonStart:codonEnd]):
            ambiguousGenotypeOrBreak = ambiguousGenotypeOrBreak + 1
            continue

        # Count stop codons, except if it occurs at the end of the sequence
        try:
            if (translationTable[firstSeq[codonStart:codonEnd]] == '*' or
                translationTable[secondSeq[codonStart:codonEnd]] == '*' or
                translationTable[outSeq[codonStart:codonEnd]] == '*'):
                stopCodons = stopCodons + 1 
        except KeyError:
            pass

        # exactly same
        if (firstSeq[codonStart:codonEnd] == secondSeq[codonStart:codonEnd] 
            and firstSeq[codonStart:codonEnd] == outSeq[codonStart:codonEnd]): continue
        # same (though differing from outgroup)
        elif firstSeq[codonStart:codonEnd] == secondSeq[codonStart:codonEnd]: continue
        # Difference, but not polarizable 
        elif (firstSeq[codonStart:codonEnd] != outSeq[codonStart:codonEnd] and
              secondSeq[codonStart:codonEnd] != outSeq[codonStart:codonEnd]):
            try:
                if translationTable[firstSeq[codonStart:codonEnd]] == translationTable[secondSeq[codonStart:codonEnd]]:
                    nonPolarizableDS = nonPolarizableDS + 1
                else:
                    nonPolarizableDN = nonPolarizableDN + 1
            except KeyError:
                translationProblems = translationProblems + 1
        # Difference in first seq
        elif secondSeq[codonStart:codonEnd] == outSeq[codonStart:codonEnd]:
            try:
                differentAA = translationTable[firstSeq[codonStart:codonEnd]]
                refAA = translationTable[outSeq[codonStart:codonEnd]]
                if differentAA != refAA:
                    if lookupTable[differentAA] < lookupTable[refAA]:
                        speciesOneDec.append(index)
                        numSpeciesOneDec = numSpeciesOneDec + 1
                    elif lookupTable[differentAA] > lookupTable[refAA]:
                        speciesOneInc.append(index)
                        numSpeciesOneInc = numSpeciesOneInc + 1
                    else: numSpeciesOneNeu = numSpeciesOneNeu + 1
            except KeyError:
                translationProblems = translationProblems + 1
        # Difference in secondSeq
        elif firstSeq[codonStart:codonEnd] == outSeq[codonStart:codonEnd]:
            try:
                differentAA = translationTable[secondSeq[codonStart:codonEnd]]
                refAA = translationTable[outSeq[codonStart:codonEnd]]
                if differentAA != refAA:
                    if lookupTable[differentAA] < lookupTable[refAA]:
                        speciesTwoDec.append(index)
                        numSpeciesTwoDec = numSpeciesTwoDec + 1
                    elif lookupTable[differentAA] > lookupTable[refAA]:
                        speciesTwoInc.append(index)
                        numSpeciesTwoInc = numSpeciesTwoInc + 1
                    else: numSpeciesTwoNeu = numSpeciesTwoNeu + 1
            except KeyError:
                translationProblems = translationProblems + 1
        else:
            assert False
    
    return {'Ambiguous Codon':ambiguousGenotypeOrBreak,'Translation Problems':translationProblems,
            'Non-polarizable DS':nonPolarizableDS, 'Non-polarizable DN':nonPolarizableDN,
            'Species One Dec Count':numSpeciesOneDec,'Species One Inc Count':numSpeciesOneInc,
            'Species One Neu Count':numSpeciesOneNeu,'Species Two Neu Count':numSpeciesTwoNeu,
            'Species Two Dec Count':numSpeciesTwoDec,'Species Two Inc Count':numSpeciesTwoInc,
            'Species One Dec':speciesOneDec,'Species One Inc':speciesOneInc,
            'Species Two Dec':speciesTwoDec,'Species Two Inc':speciesTwoInc,
            'Sequence Length':seqLength,'Stop Codon Count':stopCodons}

# not used
def clusteringNonPolarizedOld(DNList,DSList,seqLength):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    for indX, x in enumerate(DNList):
        normalizationHolder['DN'][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            variantHolder['DNDN'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            if DSList[indY] < x: indY = indY + 1; continue
            variantHolder['DNDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX, x in enumerate(DSList):
        normalizationHolder['DS'][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            variantHolder['DSDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1

    return {'Variant':variantHolder, 'Normalization': normalizationHolder}

# not used
def clusteringNonPolarized(DNList,DSList,seqLength):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 # subtracting 3 here might not be the right thing to do AMT 8 Nov 2017

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    for indX, x in enumerate(DNList):
        normalizationHolder['DN'][1:seqLength//3-x//3] += 1 # There might be a problem here also AMT 8 Nov 2017
        indY = indX + 1
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            variantHolder['DNDN'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            if DSList[indY] < x: indY = indY + 1; continue
            variantHolder['DNDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX, x in enumerate(DSList):
        normalizationHolder['DS'][1:seqLength//3-x//3] += 1 # There might be a problem here also
        indY = indX + 1
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            variantHolder['DSDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            if DNList[indY] < x: indY = indY + 1;continue
            variantHolder['DNDS'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1

    return {'Variant':variantHolder, 'Normalization': normalizationHolder}

# not used
def clusteringNonPolarizedNew(DNList,DSList,seqLength):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    for x in DNList:
        normalizationHolder['DN'][1:seqLength//3 - x//3] += 1 #This doesn't have a +1 at the end because seqLength is not 0 indexed and because the last "chance" for a pair would be one codon before the end
        normalizationHolder['DN'][1:x//3 + 1] += 1 #This has a +1 at the end because Python indexes variables like this: [start,end)
        for y in DNList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDN'][abs(y//3 - x//3)] += 1
        for y in DSList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDS'][abs(y//3 - x//3)] += 1

    for x in DSList:
        normalizationHolder['DS'][1:seqLength//3 - x//3] += 1
        normalizationHolder['DS'][1:x//3 + 1] += 1
        for y in DSList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DSDS'][abs(y//3 - x//3)] += 1
        for y in DNList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDS'][abs(y//3 - x//3)] += 1

    return {'Variant':variantHolder, 'Normalization': normalizationHolder}

# not used, but possibly still works
def clusteringAncestral(DNList,DSList,firstSeq,secondSeq):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = len(firstSeq)
    seqLength = seqLength - 3

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    for indX, x in enumerate(DNList):
        for y in range(x+3,seqLength,3):
            if y - x > 1500:
                break
            tempCodon1 = firstSeq[y:y+3]
            tempCodon2 = secondSeq[y:y+3]
            if ('-' in tempCodon1 or 'N' in tempCodon1 or 'X' in tempCodon1 or 
                '-' in tempCodon2 or 'N' in tempCodon2 or 'X' in tempCodon2):
                continue
            else:
                normalizationHolder['DN'][(y-x)//3] += 1
        #normalizationHolder['DN'][1:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            variantHolder['DNDN'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            if DSList[indY] < x: indY = indY + 1; continue
            variantHolder['DNDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX, x in enumerate(DSList):
        for y in range(x+3,seqLength,3):
            if y - x > 1500:
                break
            tempCodon1 = firstSeq[y:y+3]
            tempCodon2 = secondSeq[y:y+3]
            if ('-' in tempCodon1 or 'N' in tempCodon1 or 'X' in tempCodon1 or 
                '-' in tempCodon2 or 'N' in tempCodon2 or 'X' in tempCodon2):
                continue
            else:
                normalizationHolder['DS'][(y-x)//3] += 1
        # normalizationHolder['DS'][1:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            variantHolder['DSDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            if DNList[indY] < x: indY = indY + 1;continue
            variantHolder['DNDS'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1

    return {'Variant':variantHolder, 'Normalization': normalizationHolder}

# not used
def clusteringIntronsNotPolarizedOld(DNList, DSList, seqLength, intronPositions):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    # Clustering count
    # x is 0 indexed and always a multiple of 3
    # mutation cannot occur in a codon in which there is an intron
    for indX,x in enumerate(DNList):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DNList[indY]//3),
                                                 intronPositions/3 < DNList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DNDN'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            if DSList[indY] < x: indY = indY + 1; continue
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DSList[indY]//3),
                                                 intronPositions/3 < DSList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DNDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(DSList):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DSList[indY]//3),
                                                 intronPositions/3 < DSList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DSDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
            
    for x in DNList:
        if x >= max(intronPositions): break
        if not any(intronPositions/3 > x//3): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder['DN'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1
    for x in DSList:
        if x >= max(intronPositions): break
        if not any(intronPositions/3 > x//3): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder['DS'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1

    return {'Normalization':normalizationHolder,'Variant':variantHolder}
    #------------------------------------------------------------------------------------

# not used
def clusteringIntronsNotPolarized(DNList, DSList, seqLength, intronPositions):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])
    
    if len(intronPositions) == 0:
        return {'Normalization':normalizationHolder,'Variant':variantHolder}
    
    # Clustering count
    # x is 0 indexed and always a multiple of 3
    # mutation cannot occur in a codon in which there is an intron
    for indX,x in enumerate(DNList):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DNList[indY]//3),
                                                 intronPositions/3 < DNList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DNDN'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            if DSList[indY] < x: indY = indY + 1; continue
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DSList[indY]//3),
                                                 intronPositions/3 < DSList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DNDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(DSList):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(DSList) and DSList[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= DSList[indY]//3),
                                                 intronPositions/3 < DSList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DSDS'][DSList[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(DNList) and DNList[indY] - x <= 1500:
            if DNList[indY] < x: indY = indY + 1; continue
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions %3 == 0, intronPositions/3 <= DNList[indY]//3),
                                                 intronPositions/3 < DNList[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder['DNDS'][DNList[indY]//3 - x//3] += 1
            indY = indY + 1
            
    for x in DNList:
        if x >= max(intronPositions): break
        if not any(intronPositions/3 > x//3): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder['DN'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1
    for x in DSList:
        if x >= max(intronPositions): break
        if not any(intronPositions/3 > x//3): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder['DS'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1

    return {'Normalization':normalizationHolder,'Variant':variantHolder}
    #------------------------------------------------------------------------------------

# not used (no new function to replace it, but I should write one)
def clusteringIntronsNotPolarizedNew(DNList, DSList, seqLength, intronPositions):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    variantHolder['DNDN'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DNDS'] = numpy.array([0 for x in range(0,501)])
    variantHolder['DSDS'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DN'] = numpy.array([0 for x in range(0,501)])
    normalizationHolder['DS'] = numpy.array([0 for x in range(0,501)])

    # Clustering count
    # x is 0 indexed and always a multiple of 3
    # mutation cannot occur in a codon in which there is an intron
    for x in DNList:
        if any(intronPositions/3 > x//3):
            nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
            normalizationHolder['DN'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1
        if any(intronPositions/3 < x//3) or any(numpy.logical_and(intronPositions%3==0,intronPositions//3==x//3)):
            if any(numpy.logical_and(intronPositions%3==0,intronPositions//3==x//3)):
                nextIntron = x
            elif any(intronPositions/3 < x//3):
                nextIntron = intronPositions[numpy.argmin(intronPositions/3 < x//3)-1]
            normalizationHolder['DN'][x//3 - nextIntron//3 + 1:x//3 + 1] += 1
        lowerBoundXLogical = intronPositions/3 > x//3
        upperBoundXLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= x//3),
                                              intronPositions/3 < x//3)
        for y in DNList:
            if y < x and abs(y-x) <= 1500:
                lowerBoundLogical = intronPositions/3 > y//3
                if any(numpy.logical_and(lowerBoundLogical,upperBoundXLogical)):
                    variantHolder['DNDN'][x//3 - y//3] += 1
            if y > x and abs(y-x) <= 1500:
                upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= y//3),
                                                     intronPositions/3 < y//3)
                if any(numpy.logical_and(lowerBoundXLogical,upperBoundLogical)):
                    variantHolder['DNDN'][y//3 - x//3] += 1
        for y in DSList:
            if y < x and abs(y-x) <= 1500:
                lowerBoundLogical = intronPositions/3 > y//3
                if any(numpy.logical_and(lowerBoundLogical,upperBoundXLogical)):
                    variantHolder['DNDS'][x//3 - y//3] += 1
            if y > x and abs(y-x) <= 1500:
                upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= y//3),
                                                     intronPositions/3 < y//3)
                if any(numpy.logical_and(lowerBoundXLogical,upperBoundLogical)):
                    variantHolder['DNDS'][y//3 - x//3] += 1

    for x in DSList:
        if any(intronPositions/3 > x//3):
            nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
            normalizationHolder['DS'][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3+1] += 1
        if any(intronPositions/3 < x//3) or any(numpy.logical_and(intronPositions%3==0,intronPositions//3==x//3)):
            if any(numpy.logical_and(intronPositions%3==0,intronPositions//3==x//3)):
                nextIntron = x
            elif any(intronPositions/3 < x//3):
                nextIntron = intronPositions[numpy.argmin(intronPositions/3 < x//3)-1]
            normalizationHolder['DS'][x//3 - nextIntron//3 + 1:x//3 + 1] += 1
        lowerBoundXLogical = intronPositions/3 > x//3
        upperBoundXLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= x//3),
                                              intronPositions/3 < x//3)
        for y in DNList:
            if y < x and abs(y-x) <= 1500:
                lowerBoundLogical = intronPositions/3 > y//3
                if any(numpy.logical_and(lowerBoundLogical,upperBoundXLogical)):
                    variantHolder['DNDS'][x//3 - y//3] += 1
            if y > x and abs(y-x) <= 1500:
                upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= y//3),
                                                     intronPositions/3 < y//3)
                if any(numpy.logical_and(lowerBoundXLogical,upperBoundLogical)):
                    variantHolder['DNDS'][y//3 - x//3] += 1
        for y in DSList:
            if y < x and abs(y-x) <= 1500:
                lowerBoundLogical = intronPositions/3 > y//3
                if any(numpy.logical_and(lowerBoundLogical,upperBoundXLogical)):
                    variantHolder['DSDS'][x//3 - y//3] += 1
            if y > x and abs(y-x) <= 1500:
                upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= y//3),
                                                     intronPositions/3 < y//3)
                if any(numpy.logical_and(lowerBoundXLogical,upperBoundLogical)):
                    variantHolder['DSDS'][y//3 - x//3] += 1

    return {'Normalization':normalizationHolder,'Variant':variantHolder}
    #------------------------------------------------------------------------------------

# not used
def clusteringOld(speciesOneDN, speciesOneDS, speciesTwoDN, speciesTwoDS,
               speciesOneName, speciesTwoName, seqLength, speciesOutName = None):
    
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DSDS','DNDS','DNDN','DNDNcross']]:
        variantHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DS','DN','DNcross']]:
        normalizationHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    
    for indX,x in enumerate(speciesOneDN):
        normalizationHolder[(speciesOneName, speciesTwoName, speciesOutName,'DN',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneDN) and speciesOneDN[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDN',)][speciesOneDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            if speciesOneDS[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesOneDS):
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DS',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DSDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for indX,x in enumerate(speciesTwoDN):
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DN',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoDN) and speciesTwoDN[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDN',)][speciesTwoDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            if speciesTwoDS[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesTwoDS):
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DS',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DSDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for x in speciesTwoDN:
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNcross',)][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNcross',)][0:x//3] += 1
        for y in speciesOneDN: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
    for x in speciesOneDN:
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNcross',)][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNcross',)][0:x//3] += 1
        for y in speciesTwoDN: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
                
    return {'Normalization':normalizationHolder,'Variant':variantHolder}

# not used
def clustering(speciesOneDN, speciesOneDS, speciesTwoDN, speciesTwoDS,
               speciesOneName, speciesTwoName, seqLength, speciesOutName = None):
    
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DSDS','DNDS','DNDN','DNDNcross']]:
        variantHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DS','DN','DNcross']]:
        normalizationHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    
    for indX,x in enumerate(speciesOneDN):
        normalizationHolder[(speciesOneName, speciesTwoName, speciesOutName,'DN',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneDN) and speciesOneDN[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDN',)][speciesOneDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            if speciesOneDS[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesOneDS):
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DS',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DSDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneDN) and speciesOneDN[indY] - x <= 1500:
            if speciesOneDN[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS')][speciesOneDN[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for indX,x in enumerate(speciesTwoDN):
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DN',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoDN) and speciesTwoDN[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDN',)][speciesTwoDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            if speciesTwoDS[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesTwoDS):
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DS',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DSDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoDN) and speciesTwoDN[indY] - x <= 1500:
            if speciesTwoDN[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS')][speciesTwoDN[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for x in speciesTwoDN:
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNcross',)][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNcross',)][0:x//3] += 1
        for y in speciesOneDN: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
    for x in speciesOneDN:
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNcross',)][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNcross',)][0:x//3] += 1
        for y in speciesTwoDN: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
                
    return {'Normalization':normalizationHolder,'Variant':variantHolder}

# not used
def clusteringNew(speciesOneDN, speciesOneDS, speciesTwoDN, speciesTwoDS,
               speciesOneName, speciesTwoName, seqLength, speciesOutName = None):
    
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DSDS','DNDS','DNDN','DNDNcross']]:
        variantHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DS','DN']]:
        normalizationHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    
    for x in DNList:
        normalizationHolder['DN'][1:seqLength//3 - x//3] += 1 #This doesn't have a +1 at the end because seqLength is not 0 indexed and because the last "chance" for a pair would be one codon before the end
        normalizationHolder['DN'][1:x//3 + 1] += 1 #This has a +1 at the end because Python indexes variables like this: [start,end)
        for y in DNList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDN'][abs(y//3 - x//3)] += 1
        for y in DSList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDS'][abs(y//3 - x//3)] += 1
    for x in DSList:
        normalizationHolder['DS'][1:seqLength//3 - x//3] += 1
        normalizationHolder['DS'][1:x//3 + 1] += 1
        for y in DSList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DSDS'][abs(y//3 - x//3)] += 1
        for y in DNList:
            if abs(y-x) <= 1500 and not abs(y-x) == 0:
                variantHolder['DNDS'][abs(y//3 - x//3)] += 1

    for x in speciesOneDN:
        normalizationHolder[(speciesOneName, speciesTwoName, speciesOutName,'DN',)][1:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DN')][1:x//3 + 1] += 1
        for y in speciesOneDN:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDN',)][abs(y//3 - x//3)] += 1
        for y in speciesOneDS:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS',)][abs(y//3 - x//3)] += 1
        for y in speciesTwoDN:
            if abs(y-x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
    for x in speciesOneDS:
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DS',)][1:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoDN,speciesOutName,'DS')][1:x//3 + 1] += 1
        for y in speciesOneDS:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DSDS')][abs(y//3-x//3)] += 1
        for y in speciesOneDN:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS')][abs(y//3-x//3)] += 1

    for x in speciesTwoDN:
        normalizationHolder[(speciesTwoName, speciesOneName, speciesOutName,'DN',)][1:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DN')][1:x//3 + 1] += 1
        for y in speciesTwoDN:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDN',)][abs(y//3 - x//3)] += 1
        for y in speciesTwoDS:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS',)][abs(y//3 - x//3)] += 1
        for y in speciesOneDN:
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDNcross',)][abs(y//3 - x//3)] += 1
    for x in speciesTwoDS:
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DS',)][1:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneDN,speciesOutName,'DS')][1:x//3 + 1] += 1
        for y in speciesTwoDS:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DSDS')][abs(y//3-x//3)] += 1
        for y in speciesTwoDN:
            if abs(y-x) <= 1500 and not x == y:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS')][abs(y//3-x//3)] += 1
                
    return {'Normalization':normalizationHolder,'Variant':variantHolder}

# no idea what this is, probably not used
def clusteringIntrons(speciesOneDN, speciesOneDS, speciesTwoDN, speciesTwoDS,
                      speciesOneName, speciesTwoName, seqLength, intronPositions, 
                      speciesOutName = None):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength = seqLength - 3 #I don't look for mutations in the last codon so this adjusts things accordingly 

    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DSDS','DNDS','DNDN','DNDNcross']]:
        variantHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['DS','DN','DNcross']]:
        normalizationHolder[tupleKey] = numpy.array([0 for x in range(0,501)])

    # Clustering count
    # x is 0 indexed and always a multiple of 3
    # mutation cannot occur in a codon in which there is an intron
    for indX,x in enumerate(speciesOneDN):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(speciesOneDN) and speciesOneDN[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesOneDN[indY]//3),
                                                 intronPositions/3 < speciesOneDN[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDN',)][speciesOneDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            if speciesOneDS[indY] < x: indY = indY + 1; continue
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesOneDS[indY]//3),
                                                 intronPositions/3 < speciesOneDS[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DNDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesOneDS):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(speciesOneDS) and speciesOneDS[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesOneDS[indY]//3),
                                                 intronPositions/3 < speciesOneDS[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'DSDS',)][speciesOneDS[indY]//3 - x//3] += 1
            indY = indY + 1

    for indX,x in enumerate(speciesTwoDN):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(speciesTwoDN) and speciesTwoDN[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesTwoDN[indY]//3),
                                                 intronPositions/3 < speciesTwoDN[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDN',)][speciesTwoDN[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            if speciesTwoDS[indY] < x: indY = indY + 1; continue
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesTwoDS[indY]//3),
                                                 intronPositions/3 < speciesTwoDS[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DNDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesTwoDS):
        indY = indX + 1
        lowerBoundLogical = intronPositions/3 > x//3
        while indY < len(speciesTwoDS) and speciesTwoDS[indY] - x <= 1500:
            upperBoundLogical = numpy.logical_or(numpy.logical_and(intronPositions % 3 == 0,intronPositions/3 <= speciesTwoDS[indY]//3),
                                                 intronPositions/3 < speciesTwoDS[indY]//3)
            if any(numpy.logical_and(lowerBoundLogical,upperBoundLogical)):
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'DSDS',)][speciesTwoDS[indY]//3 - x//3] += 1
            indY = indY + 1
            
    for x in speciesOneDN:
        if x >= max(intronPositions): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DN',)][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3] += 1
    for x in speciesOneDS:
        if x >= max(intronPositions): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'DS',)][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3] += 1
    for x in speciesTwoDN:
        if x >= max(intronPositions): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DN',)][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3] += 1
    for x in speciesTwoDS:
        if x >= max(intronPositions): break
        nextIntron = intronPositions[numpy.argmax(intronPositions/3 > x//3)]
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'DS',)][nextIntron//3 - x//3 + (nextIntron%3 != 0):(seqLength-1)//3-x//3] += 1

    return {'Normalization':normalizationHolder,'Variant':variantHolder}

def propertyClustering(speciesOneInc,speciesOneDec,speciesTwoInc,speciesTwoDec,
                       speciesOneName,speciesTwoName,seqLength,speciesOutName = None):
    normalizationHolder = dict()
    variantHolder = dict()
    seqLength - seqLength - 3
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['comp','rein','compCross','reinCross']]:
        variantHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    for tupleKey in [x + (speciesOutName,) + (y,) for x in [(speciesOneName,speciesTwoName),(speciesTwoName,speciesOneName)] for y in ['within','cross']]:
        normalizationHolder[tupleKey] = numpy.array([0 for x in range(0,501)])
    
    for indX,x in enumerate(speciesOneInc):
        normalizationHolder[(speciesOneName, speciesTwoName, speciesOutName,'within',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneInc) and speciesOneInc[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'rein',)][speciesOneInc[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneDec) and speciesOneDec[indY] - x <= 1500:
            if speciesOneDec[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'comp',)][speciesOneDec[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesOneDec):
        normalizationHolder[(speciesOneName, speciesTwoName, speciesOutName,'within',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesOneDec) and speciesOneDec[indY] - x <= 1500:
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'rein',)][speciesOneDec[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesOneInc) and speciesOneInc[indY] - x <= 1500:
            if speciesOneInc[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'comp')][speciesOneInc[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for indX,x in enumerate(speciesTwoInc):
        normalizationHolder[(speciesTwoName, speciesOneName, speciesOutName,'within',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoInc) and speciesTwoInc[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'rein',)][speciesTwoInc[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoDec) and speciesTwoDec[indY] - x <= 1500:
            if speciesTwoDec[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'comp',)][speciesTwoDec[indY]//3 - x//3] += 1
            indY = indY + 1
    for indX,x in enumerate(speciesTwoDec):
        normalizationHolder[(speciesTwoName, speciesOneName, speciesOutName,'within',)][0:seqLength//3-x//3] += 1
        indY = indX + 1
        while indY < len(speciesTwoDec) and speciesTwoDec[indY] - x <= 1500:
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'rein',)][speciesTwoDec[indY]//3 - x//3] += 1
            indY = indY + 1
        indY = 0
        while indY < len(speciesTwoInc) and speciesTwoInc[indY] - x <= 1500:
            if speciesTwoInc[indY] < x: indY = indY + 1; continue
            variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'comp')][speciesTwoInc[indY]//3 - x//3] += 1
            indY = indY + 1
    
    for x in speciesOneInc:
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'cross')][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'cross')][0:x//3] += 1
        for y in speciesTwoInc: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'reinCross')][abs(y//3 - x//3)] += 1
        for y in speciesTwoDec:
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'compCross')][abs(y//3 - x//3)] += 1

    for x in speciesOneDec:
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'cross')][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesOneName,speciesTwoName,speciesOutName,'cross')][0:x//3] += 1
        for y in speciesTwoDec: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'reinCross')][abs(y//3 - x//3)] += 1
        for y in speciesTwoInc: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesOneName,speciesTwoName,speciesOutName,'compCross')][abs(y//3 - x//3)] += 1

    for x in speciesTwoInc:
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'cross')][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'cross')][0:x//3] += 1
        for y in speciesOneInc: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'reinCross')][abs(y//3 - x//3)] += 1
        for y in speciesOneDec:
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'compCross')][abs(y//3 - x//3)] += 1

    for x in speciesTwoDec:
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'cross')][0:seqLength//3-x//3] += 1
        normalizationHolder[(speciesTwoName,speciesOneName,speciesOutName,'cross')][0:x//3] += 1
        for y in speciesOneDec: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'reinCross')][abs(y//3 - x//3)] += 1
        for y in speciesOneInc: 
            if abs(y - x) <= 1500:
                variantHolder[(speciesTwoName,speciesOneName,speciesOutName,'compCross')][abs(y//3 - x//3)] += 1
    
    return {'Normalization':normalizationHolder,'Variant':variantHolder}

def graphPropertyClusters(variantProperty,normalizationProperty,variant,normalization,
                          speciesOneName,speciesTwoName,propertyName,speciesOutName = None
                          ,graphUpperLimit=400,normLowerLimit=300,normUpperLimit=350,save=False,savePath='',
                          fileSpecificName='',smoothing = False,window = 5,calculateP = False,expectation='clustering'):
    ''' Less intuitive graphing
    speciesOneComp = ((variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'comp')][1:500].astype(float)/
                       normalizationProperty[(speciesOneName,speciesTwoName,speciesOutName,'within')][1:500])/
                      (variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][1:500].astype(float)/
                       normalization[(speciesOneName,speciesTwoName,speciesOutName,'DN')][1:500]))
    speciesTwoComp = ((variantProperty[(speciesTwoName,speciesOneName,speciesOutName,'comp')][1:500].astype(float)/
                       normalizationProperty[(speciesTwoName,speciesOneName,speciesOutName,'within')][1:500])/
                      (variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][1:500].astype(float)/
                       normalization[(speciesTwoName,speciesOneName,speciesOutName,'DN')][1:500]))
    speciesOneRein = ((variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'rein')][1:500].astype(float)/
                       normalizationProperty[(speciesOneName,speciesTwoName,speciesOutName,'within')][1:500])/
                      (variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][1:500].astype(float)/
                       normalization[(speciesOneName,speciesTwoName,speciesOutName,'DN')][1:500]))
    speciesTwoRein = ((variantProperty[(speciesTwoName,speciesOneName,speciesOutName,'rein')][1:500].astype(float)/
                       normalizationProperty[(speciesTwoName,speciesOneName,speciesOutName,'within')][1:500])/
                      (variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][1:500].astype(float)/
                       normalization[(speciesTwoName,speciesOneName,speciesOutName,'DN')][1:500]))
    crossSpeciesComp = ((variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'compCross')][1:500].astype(float)/
                         (normalizationProperty[(speciesOneName,speciesTwoName,speciesOutName,'cross')][1:500] + normalizationProperty[(speciesTwoName,speciesOneName,speciesOutName,'cross')][1:500]))/
                        (variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')][1:500].astype(float)/
                         (normalization[(speciesOneName,speciesTwoName,speciesOutName,'DN')][1:500] + normalization[(speciesTwoName,speciesOneName,speciesOutName,'DN')][1:500])))
    crossSpeciesRein = ((variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'reinCross')][1:500].astype(float)/
                         (normalizationProperty[(speciesOneName,speciesTwoName,speciesOutName,'cross')][1:500] + normalizationProperty[(speciesTwoName,speciesOneName,speciesOutName,'cross')][1:500]))/
                        (variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][1:500].astype(float)/
                         (normalization[(speciesOneName,speciesTwoName,speciesOutName,'DN')][1:500] + normalization[(speciesTwoName,speciesOneName,speciesOutName,'DN')][1:500])))
    
    fig = plt.figure()
    fig.suptitle(speciesOneName + '-' + speciesTwoName + " Outgroup: " + speciesOutName + " "  + propertyName + " Clustering")
    ax = fig.add_subplot(111)
    ax.set_xlabel('Distance (Codons)')
    ax.set_ylabel('Fraction')
    ax.plot(range(1,graphUpperLimit + 1),speciesOneComp[1:graphUpperLimit + 1],'g-',label=speciesOneName + ' Compensatory')
    ax.plot(range(1,graphUpperLimit + 1),speciesOneRein[1:graphUpperLimit + 1],'g--',label=speciesOneName + ' Reinforcing')
    ax.plot(range(1,graphUpperLimit + 1),speciesTwoComp[1:graphUpperLimit + 1],'y-',label=speciesTwoName + ' Compensatory')
    ax.plot(range(1,graphUpperLimit + 1),speciesTwoRein[1:graphUpperLimit + 1],'y--',label=speciesTwoName + ' Reinforcing')
    ax.plot(range(1,graphUpperLimit + 1),crossSpeciesComp[1:graphUpperLimit + 1]*2,'b-',label='Cross Species Compensatory')
    ax.plot(range(1,graphUpperLimit + 1),crossSpeciesRein[1:graphUpperLimit + 1]*2,'b--',label='Cross Species Reinforcing')
    plt.ylim(ymax = (max([max(speciesOneComp),max(speciesOneRein),max(speciesTwoComp),max(speciesTwoRein),max(crossSpeciesComp)*2,max(crossSpeciesRein)*2])*1.5))
    plt.legend()
    if save:
        plt.savefig(savePath + propertyName + "_" + speciesOneName + '-' + speciesTwoName + '-' + speciesOutName + '_normalized.png',bbox_inches = 'tight')
        plt.close(fig)
    '''

    if calculateP:
        pvals = propertyClusteringSignificance(variantProperty,variant,speciesOneName,speciesTwoName,speciesOutName,
                                               expectation=expectation,upperLimit = 50)

    speciesOneComp = (variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'comp')][1:500].astype(float)/
                      variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][1:500].astype(float))
    speciesTwoComp = (variantProperty[(speciesTwoName,speciesOneName,speciesOutName,'comp')][1:500].astype(float)/
                      variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][1:500].astype(float))
    speciesOneRein = (variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'rein')][1:500].astype(float)/
                      variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][1:500].astype(float))
    speciesTwoRein = (variantProperty[(speciesTwoName,speciesOneName,speciesOutName,'rein')][1:500].astype(float)/
                      variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][1:500].astype(float))
    crossSpeciesComp = (variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'compCross')][1:500].astype(float)/
                        variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')][1:500].astype(float))
    crossSpeciesRein = (variantProperty[(speciesOneName,speciesTwoName,speciesOutName,'reinCross')][1:500].astype(float)/
                        variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][1:500].astype(float))

    # Smoothing here
    if smoothing:
        speciesOneComp = windowSmoothing(speciesOneComp,window)
        speciesTwoComp = windowSmoothing(speciesTwoComp,window)
        speciesOneRein = windowSmoothing(speciesOneRein,window)
        speciesTwoRein = windowSmoothing(speciesTwoRein,window)
        crossSpeciesComp = windowSmoothing(crossSpeciesComp,window)
        crossSpeciesRein = windowSmoothing(crossSpeciesRein,window)

    if calculateP:
        fig = plt.figure(figsize=(11,6))
    else:
        fig = plt.figure()
    fig.suptitle(speciesOneName + '-' + speciesTwoName + " Outgroup: " + speciesOutName + " "  + propertyName + " Clustering Fraction")
    ax = fig.add_subplot(111)
    ax.set_xlabel('Distance (Codons)')
    ax.set_ylabel('Fraction')
    if calculateP:
        ax.plot(range(1,graphUpperLimit + 1),speciesOneComp[1:graphUpperLimit + 1],'g-',label=speciesOneName + ' Compensatory p = ' + str(round(pvals[(speciesOneName,speciesTwoName,speciesOutName,'comp')],2)))
        ax.plot(range(1,graphUpperLimit + 1),speciesOneRein[1:graphUpperLimit + 1],'g--',label=speciesOneName + ' Reinforcing p = ' + str(round(pvals[(speciesOneName,speciesTwoName,speciesOutName,'rein')],2)))
        ax.plot(range(1,graphUpperLimit + 1),speciesTwoComp[1:graphUpperLimit + 1],'y-',label=speciesTwoName + ' Compensatory p = ' + str(round(pvals[(speciesTwoName,speciesOneName,speciesOutName,'comp')],2)))
        ax.plot(range(1,graphUpperLimit + 1),speciesTwoRein[1:graphUpperLimit + 1],'y--',label=speciesTwoName + ' Reinforcing p = ' + str(round(pvals[(speciesTwoName,speciesOneName,speciesOutName,'rein')],2)))
        ax.plot(range(1,graphUpperLimit + 1),crossSpeciesComp[1:graphUpperLimit + 1],'b-',label='Cross Species Compensatory p = ' + str(round(pvals[(speciesOneName,speciesTwoName,speciesOutName,'compCross')],2)))
        ax.plot(range(1,graphUpperLimit + 1),crossSpeciesRein[1:graphUpperLimit + 1],'b--',label='Cross Species Reinforcing p = ' + str(round(pvals[(speciesOneName,speciesTwoName,speciesOutName,'reinCross')],2)))
    else:
        ax.plot(range(1,graphUpperLimit + 1),speciesOneComp[1:graphUpperLimit + 1],'g-',label=speciesOneName + ' Compensatory')
        ax.plot(range(1,graphUpperLimit + 1),speciesOneRein[1:graphUpperLimit + 1],'g--',label=speciesOneName + ' Reinforcing')
        ax.plot(range(1,graphUpperLimit + 1),speciesTwoComp[1:graphUpperLimit + 1],'y-',label=speciesTwoName + ' Compensatory')
        ax.plot(range(1,graphUpperLimit + 1),speciesTwoRein[1:graphUpperLimit + 1],'y--',label=speciesTwoName + ' Reinforcing')
        ax.plot(range(1,graphUpperLimit + 1),crossSpeciesComp[1:graphUpperLimit + 1],'b-',label='Cross Species Compensatory')
        ax.plot(range(1,graphUpperLimit + 1),crossSpeciesRein[1:graphUpperLimit + 1],'b--',label='Cross Species Reinforcing')
    maxVal = max([max(speciesOneComp[1:graphUpperLimit + 1]),max(speciesOneRein[1:graphUpperLimit + 1]),
                  max(speciesTwoComp[1:graphUpperLimit + 1]),max(speciesTwoRein[1:graphUpperLimit + 1]),
                  max(crossSpeciesComp[1:graphUpperLimit + 1]),max(crossSpeciesRein[1:graphUpperLimit + 1])])
    minVal = min([min(speciesOneComp[1:graphUpperLimit + 1]),min(speciesOneRein[1:graphUpperLimit + 1]),
                  min(speciesTwoComp[1:graphUpperLimit + 1]),min(speciesTwoRein[1:graphUpperLimit + 1]),
                  min(crossSpeciesComp[1:graphUpperLimit + 1]),min(crossSpeciesRein[1:graphUpperLimit + 1])])
    plt.ylim(ymax = (maxVal + .3*(maxVal-minVal)))
    #plt.legend(loc=9,ncol=2)
    plt.legend(bbox_to_anchor=(0.,.898,1.,.102), loc=6,ncol=2, mode='expand',borderaxespad=0.)
    if save:
        if fileSpecificName != '':
            plt.savefig(savePath + '/' + fileSpecificName + '_' + propertyName + "_" + speciesOneName + '-' + speciesTwoName + '-' + speciesOutName + '_raw.png',bbox_inches = 'tight')
        else:
            plt.savefig(savePath + '/' + propertyName + "_" + speciesOneName + '-' + speciesTwoName + '-' + speciesOutName + '_raw.png',bbox_inches = 'tight')
        plt.close(fig)
    else:
        plt.show()
    return

#def propertyClusteringSignificance(propertyVariant,variant,speciesOneName,speciesTwoName,
#                                   speciesOutName=None,propertyName='',lowerLimit=1,upperLimit=501,
#                                   expectation = None, clusteringLowerLimit = 1, clusteringUpperLimit = 20):
#    listToReturn = dict()
#    listOfKeys = [((speciesOneName,speciesTwoName,speciesOutName,'comp'),(speciesOneName,speciesTwoName,speciesOutName,'DNDN')),
#                  ((speciesOneName,speciesTwoName,speciesOutName,'rein'),(speciesOneName,speciesTwoName,speciesOutName,'DNDN')),
#                  ((speciesTwoName,speciesOneName,speciesOutName,'comp'),(speciesTwoName,speciesOneName,speciesOutName,'DNDN')),
#                  ((speciesTwoName,speciesOneName,speciesOutName,'rein'),(speciesTwoName,speciesOneName,speciesOutName,'DNDN')),
#                  ((speciesOneName,speciesTwoName,speciesOutName,'compCross'),(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')),
#                  ((speciesOneName,speciesTwoName,speciesOutName,'reinCross'),(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross'))]
#    for key in listOfKeys:
#        scalingFactor = (sum(propertyVariant[key[0]][lowerLimit:upperLimit])/float(sum(variant[key[1]][lowerLimit:upperLimit])))
#        if expectation == None or expectation == 'Uniform':
#            null = scalingFactor * variant[key[1]][lowerLimit:upperLimit]
#            if sum(null < 5) > .2 * len(null):
#                print 'Warning, expectation values are low, error in chi-square may result'
#            listToReturn[key[0]] = stats.chisquare(propertyVariant[key[0]][lowerLimit:upperLimit],null)[1]
#        elif expectation == 'clustering':
#            null = scalingFactor * variant[key[1]]
#            clusteringExp = sum(null[clusteringLowerLimit:clusteringUpperLimit])
#            everythingExp = sum(null[clusteringUpperLimit+1:upperLimit])
#            clusteringObs = sum(propertyVariant[key[0]][clusteringLowerLimit:clusteringUpperLimit])
#            everythingObs = sum(propertyVariant[key[0]][clusteringUpperLimit+1:upperLimit])
#            listToReturn[key[0]] = stats.chisquare([clusteringObs,everythingObs],[clusteringExp,everythingExp])[1]
#        else:
#            raise ValueError("'" + expectation + "' is unsupported by this function")
#
#    return listToReturn

def clusteringSignificance(variant,normalization,speciesOneName,speciesTwoName,speciesOutName=None,
                           lowerLimit=1,upperLimit=501,options=7,expectation=None,clusteringLowerLimit = 1,
                           clusteringUpperLimit = 20):
    # options: ones place = species one p-val is calculated, twos place = species two p-val is calculated,
    # fours place = between species p-val is calculated
    # default null model is uniform distribution of observations over clustering range, other options
    # are 'longRange' or 'DSDS'. In the case of 'longRange', DOF=number categories, for all others,
    # DOF=number categories - 1. For within/between species clustering significance, the within species is the
    # null model. 
    listToCalc = []
    listToReturn = dict()
    if options == 1 or options == 3 or options == 5 or options == 7:
        listToCalc.append((speciesOneName,speciesTwoName,speciesOutName))
    if options == 2 or options == 3 or options == 6 or options == 7:
        listToCalc.append((speciesTwoName,speciesOneName,speciesOutName))

    # Prevents divide by zero erros, and the first value isn't used
    for key in normalization:
        normalization[key][0] = 1

    for key in listToCalc:
        if expectation == None or expectation == 'uniform':
            scalingFactor = (sum(variant[key + ('DNDN',)][lowerLimit:upperLimit])/
                             float(sum(normalization[key + ('DN',)][lowerLimit:upperLimit])))
            null = scalingFactor * normalization[key + ('DN',)][lowerLimit:upperLimit]
            DDOF = 0
            listToReturn[key + ('within',)] = stats.chisquare(variant[key + ('DNDN',)][lowerLimit:upperLimit],null,DDOF)[1]
        elif expectation == 'DSDS':
            scalingFactor = (sum(variant[key + ('DNDN',)][lowerLimit:upperLimit])/
                             float(sum(variant[key + ('DSDS',)][lowerLimit:upperLimit])))
            null = scalingFactor * variant[key + ('DSDS',)][lowerLimit:upperLimit]
            DDOF = 0
            listToReturn[key + ('within',)] = stats.chisquare(variant[key + ('DNDN',)][lowerLimit:upperLimit],null,DDOF)[1]
        elif expectation == 'clustering':
            scalingFactor = (sum(variant[key + ('DNDN',)][lowerLimit:upperLimit])/
                             float(sum(normalization[key + ('DN',)][lowerLimit:upperLimit])))
            null = scalingFactor * normalization[key + ('DN',)]
            clusteringExp = sum(null[clusteringLowerLimit:clusteringUpperLimit])
            everythingExp = sum(null[clusteringUpperLimit+1:upperLimit])
            clusteringObs = sum(variant[key + ('DNDN',)][clusteringLowerLimit:clusteringUpperLimit])
            everythingObs = sum(variant[key + ('DNDN',)][clusteringUpperLimit+1:upperLimit])
            listToReturn[key + ('within',)] = stats.chisquare([clusteringObs,everythingObs],[clusteringExp,everythingExp])[1]
        else:
            raise ValueError("'" + expectation + "' is unsupported by this function")

    if options == 4 or options == 5 or options == 6 or options == 7:
        scalingFactor = (sum(variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')][lowerLimit:upperLimit])/
                         float(sum(variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][lowerLimit:upperLimit])))
        null = scalingFactor * variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')][lowerLimit:upperLimit]
        listToReturn[(speciesOneName,speciesTwoName,speciesOutName,'between',)] = stats.chisquare(variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][lowerLimit:upperLimit],null)[1]

        scalingFactor = (sum(variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')][lowerLimit:upperLimit])/
                         float(sum(variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][lowerLimit:upperLimit])))
        null = scalingFactor * variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')][lowerLimit:upperLimit]
        listToReturn[(speciesTwoName,speciesOneName,speciesOutName,'between',)] = stats.chisquare(variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross',)][lowerLimit:upperLimit],null)[1]

    return listToReturn

def graphClusters(variant,normalization,speciesOneName,speciesTwoName,speciesOutName=None,
                  graphingOptions=7,graphUpperLimit=400,normLowerLimit=300,normUpperLimit=350,save=False,savePath='',
                  fileSpecificName='',smoothing = False, window = 5, calculateP = False, expectation = 'clustering'):
    # graphing options, ones place = species one is graphed, twos place = species two is graphed, fours place = between species graphed
    graphList = []
    if graphingOptions == 1 or graphingOptions == 3 or graphingOptions == 5 or graphingOptions == 7:
        graphList.append((speciesOneName,speciesTwoName,speciesOutName))
    if graphingOptions == 2 or graphingOptions == 3 or graphingOptions == 6 or graphingOptions == 7:
        graphList.append((speciesTwoName,speciesOneName,speciesOutName))
    if calculateP:
        pvals = clusteringSignificance(variant,normalization,speciesOneName,speciesTwoName,speciesOutName,expectation=expectation,upperLimit=graphUpperLimit)

    # Prevents divide by zero erros, and the first value isn't graphed anyway
    for key in normalization:
        normalization[key][0] = 1
    for key in graphList:
        DNDN = variant[key + ('DNDN',)].astype(float)/normalization[key + ('DN',)]
        DNDS = variant[key + ('DNDS',)].astype(float)/normalization[key + ('DN',)]
        DSDS = variant[key + ('DSDS',)].astype(float)/normalization[key + ('DS',)]

        normDNDN = DNDN/numpy.mean(DNDN[normLowerLimit:normUpperLimit])
        normDNDS = DNDS/numpy.mean(DNDS[normLowerLimit:normUpperLimit])
        normDSDS = DSDS/numpy.mean(DSDS[normLowerLimit:normUpperLimit])
        
        fig = plt.figure()
        if speciesOutName == None:
            fig.suptitle(key[0] + "Ancestral Comparison")
        else:
            fig.suptitle(key[0] + " vs " + key[1] + ", " + " Outgroup: " + key[2])
        ax = fig.add_subplot(111)

        # Smoothing here
        if smoothing:
            normDNDN = windowSmoothing(normDNDN,window)
            normDNDS = windowSmoothing(normDNDS,window)
            normDSDS = windowSmoothing(normDSDS,window)

        if calculateP:
            ax.plot(range(1,graphUpperLimit + 1),normDNDN[1:graphUpperLimit + 1],'g-',label=('DNDN p = ' + str(pvals[key + ('within',)])))
            ax.plot(range(1,graphUpperLimit + 1),normDNDS[1:graphUpperLimit + 1],'y-',label='DNDS')
            ax.plot(range(1,graphUpperLimit + 1),normDSDS[1:graphUpperLimit + 1],'b-',label='DSDS')
        else:
            ax.plot(range(1,graphUpperLimit + 1),normDNDN[1:graphUpperLimit + 1],'g-',label='DNDN')
            ax.plot(range(1,graphUpperLimit + 1),normDNDS[1:graphUpperLimit + 1],'y-',label='DNDS')
            ax.plot(range(1,graphUpperLimit + 1),normDSDS[1:graphUpperLimit + 1],'b-',label='DSDS')
        oneLine = ax.plot(range(1,graphUpperLimit + 1),numpy.ones(graphUpperLimit),'-')
        plt.setp(oneLine,'linestyle','--','color','black')
        plt.xlabel('Distance (codons)')
        plt.ylabel('Normalized Conditional Probability of Second Mutation')
        plt.legend()
        if save:
            if fileSpecificName != '':
                plt.savefig(savePath + '/' + fileSpecificName + '_' + key[0] + '-' + key[1] + '-' + key[2] + '.png',bbox_inches='tight')
            else:
                plt.savefig(savePath + '/' + key[0] + '-' + key[1] + '-' + key[2] + '.png',bbox_inches='tight')
            plt.close(fig)
    
    if graphingOptions == 4 or graphingOptions == 5 or graphingOptions == 6 or graphingOptions == 7:
        fig = plt.figure()
        fig.suptitle(speciesOneName + ", " + speciesTwoName + "; Outgroup: " + speciesOutName)
        ax = fig.add_subplot(111)
        DNDNOne = variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDN')].astype(float)/normalization[(speciesOneName,speciesTwoName,speciesOutName,'DN')]
        normDNDNOne = DNDNOne/numpy.mean(DNDNOne[normLowerLimit:normUpperLimit])
        DNDNTwo = variant[(speciesTwoName,speciesOneName,speciesOutName,'DNDN')].astype(float)/normalization[(speciesTwoName,speciesOneName,speciesOutName,'DN')]
        normDNDNTwo = DNDNTwo/numpy.mean(DNDNTwo[normLowerLimit:normUpperLimit])
        DNDNBetween = variant[(speciesOneName,speciesTwoName,speciesOutName,'DNDNcross')].astype(float)/normalization[(speciesOneName,speciesTwoName,speciesOutName,'DNcross')]
        normDNDNBetween = DNDNBetween / numpy.mean(DNDNBetween[normLowerLimit:normUpperLimit])

        # Smoothing here
        if smoothing:
            normDNDNOne = windowSmoothing(normDNDNOne,window)
            normDNDNTwo = windowSmoothing(normDNDNTwo,window)
            normDNDNBetween = windowSmoothing(normDNDNBetween,window)
        if calculateP:
            ax.plot(range(1,graphUpperLimit + 1),normDNDNOne[1:graphUpperLimit + 1],'g-',label='Within ' + speciesOneName + ' DNDN p = ' + str(pvals[(speciesOneName,speciesTwoName,speciesOutName,'between')]))
            ax.plot(range(1,graphUpperLimit + 1),normDNDNTwo[1:graphUpperLimit + 1],'r-',label='Within ' + speciesTwoName + ' DNDN p = ' + str(pvals[(speciesTwoName,speciesOneName,speciesOutName,'between')]))
            ax.plot(range(1,graphUpperLimit + 1),normDNDNBetween[1:graphUpperLimit + 1],'b-',label='Between ' + speciesOneName + ' and ' + speciesTwoName + ' DNDN')
        else:
            ax.plot(range(1,graphUpperLimit + 1),normDNDNOne[1:graphUpperLimit + 1],'g-',label='Within ' + speciesOneName + ' DNDN')
            ax.plot(range(1,graphUpperLimit + 1),normDNDNTwo[1:graphUpperLimit + 1],'r-',label='Within ' + speciesTwoName + ' DNDN')
            ax.plot(range(1,graphUpperLimit + 1),normDNDNBetween[1:graphUpperLimit + 1],'b-',label='Between ' + speciesOneName + ' and ' + speciesTwoName + ' DNDN')
        oneLine = ax.plot(range(1,graphUpperLimit + 1),numpy.ones(graphUpperLimit),'-')
        plt.setp(oneLine,'linestyle','--','color','black')
        plt.xlabel('Distance (codons)')
        plt.ylabel('Normalized Conditional Probability of Second Mutation')
        plt.legend()
        if save:
            if fileSpecificName != '':
                plt.savefig(savePath + '/' + fileSpecificName + '_' + speciesOneName + '-' + speciesTwoName + '-' + speciesOutName + '_cross.png',bbox_inches='tight')
            else:
                plt.savefig(savePath + '/' + speciesOneName + '-' + speciesTwoName + '-' + speciesOutName + '_cross.png',bbox_inches='tight')
            plt.close(fig)
    if not save:
        plt.show()
    return

def graphNotPolarizedClusters(variant,normalization,speciesOneName,speciesTwoName,
                          graphUpperLimit=400,normLowerLimit=300,normUpperLimit=350,save=False,savePath='',
                          fileSpecificName='',smoothing = False, window = 5):

    # Prevents divide by zero erros, and the first value isn't graphed anyway
    for key in normalization:
        normalization[key][0] = 1

    DNDN = variant['DNDN'].astype(float)/normalization['DN']
    DNDS = variant['DNDS'].astype(float)/(normalization['DN'] + normalization['DS'])
    DSDS = variant['DSDS'].astype(float)/normalization['DS']
    
    normDNDN = DNDN/numpy.mean(DNDN[normLowerLimit:normUpperLimit])
    normDNDS = DNDS/numpy.mean(DNDS[normLowerLimit:normUpperLimit])
    normDSDS = DSDS/numpy.mean(DSDS[normLowerLimit:normUpperLimit])
    
    fig = plt.figure()
    fig.suptitle(speciesOneName + " vs " + speciesTwoName + " Not Polarized")
    ax = fig.add_subplot(111)

    # Smoothing here
    if smoothing:
        normDNDN = windowSmoothing(normDNDN,window)
        normDNDS = windowSmoothing(normDNDS,window)
        normDSDS = windowSmoothing(normDSDS,window)

    ax.plot(range(1,graphUpperLimit + 1),normDNDN[1:graphUpperLimit + 1],'g-',label='DNDN')
    ax.plot(range(1,graphUpperLimit + 1),normDNDS[1:graphUpperLimit + 1],'y-',label='DNDS')
    ax.plot(range(1,graphUpperLimit + 1),normDSDS[1:graphUpperLimit + 1],'b-',label='DSDS')
    oneLine = ax.plot(range(1,graphUpperLimit + 1),numpy.ones(graphUpperLimit),'-')
    plt.setp(oneLine,'linestyle','--','color','black')
    plt.xlabel('Distance (codons)')
    plt.ylabel('Normalized Conditional Probability of Second Mutation')
    plt.legend()
    if save:
        if fileSpecificName != '':
            plt.savefig(savePath + '/' + fileSpecificName + '_' + speciesOneName + '-' + speciesTwoName + '_notPolarized.png',bbox_inches='tight')
        else:
            plt.savefig(savePath + '/' + speciesOneName + '-' + speciesTwoName + '_notPolarized.png',bbox_inches='tight')
        plt.close(fig)
    
    if not save:
        plt.show()
    return


def windowSmoothing(data,windowSize=5):
    if windowSize%2 != 1:
        raise ValueError('Window size must be odd')
    # Cheating below, kind of. What I really want is floor and ceiling, respectively, but this does the job.
    deltaLower = windowSize//2
    deltaUpper = windowSize//2 + 1
    return [sum(data[max([i-deltaLower,0]):min([i+deltaUpper,len(data)])])/float(len(data[max([i-deltaLower,0]):min([i+deltaUpper,len(data)])])) for i in range(0,len(data))]

def clusteringMetric(variant,normalization,clusteringLowerLimit=1,
                     clusteringUpperLimit=100,normLowerLimit=300,
                     normUpperLimit=350):
    # Restricting the length of the vector should elimintate numpy divide by zero warnings
    # if the selected range to calulate the "clustering metric" is ok
    clustering = variant[:normUpperLimit].astype(float)/normalization[:normUpperLimit]

    return sum(clustering[clusteringLowerLimit:clusteringUpperLimit] - numpy.mean(clustering[normLowerLimit:normUpperLimit]))


