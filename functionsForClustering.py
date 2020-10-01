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

# Used 
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

# Used
def clusteringSameMutation(mutation_positions,clustering_calc_length = 500):
    mutation_positions = numpy.array(mutation_positions)//3
    clusteringDistances = mutation_positions - mutation_positions[:,None]
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances > 0,
        clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    return clusteringCount

# Used
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

# Used
def clusteringDifferentMutation(mutation_positions_1,mutation_positions_2,clustering_calc_length = 500):
    mutation_positions_1 = numpy.array(mutation_positions_1)//3
    mutation_positions_2 = numpy.array(mutation_positions_2)//3
    clusteringDistances = numpy.abs(mutation_positions_1 - mutation_positions_2[:,None])
    clusteringDistances = clusteringDistances[numpy.logical_and(clusteringDistances > 0,
        clusteringDistances <= clustering_calc_length)]
    clusteringCount = numpy.histogram(clusteringDistances,bins = range(clustering_calc_length + 2))[0]
    return clusteringCount

# Used
def analyticalNormSameMut(numMutations, geneLength, clustering_calc_length = 500):
    codonRange = numpy.arange(clustering_calc_length + 1)
    numMutations = float(numMutations)
    tempHolder = (numMutations*(numMutations-1))/(geneLength-1) - (numMutations*(numMutations-1))/(geneLength*(geneLength - 1))*codonRange
    #tempHolder = -((float(numMutations)**2-numMutations)/geneLength**2)*codonRange + (float(numMutations)**2-numMutations)/geneLength
    tempHolder[tempHolder < 0] = 0
    return tempHolder

# Used
def analyticalNormDiffMut(numMutations1, numMutations2, geneLength, clustering_calc_length = 500):
    codonRange = numpy.arange(clustering_calc_length + 1)
    numMutations1 = float(numMutations1)
    numMutations2 = float(numMutations2)
    tempHolder = ((2*numMutations1*numMutations2)/(geneLength-1)) - ((2*numMutations1*numMutations2)/(geneLength*(geneLength-1)))*codonRange
    #tempHolder = -(2*float(numMutations1*numMutations2)/geneLength**2)*codonRange + 2*float(numMutations1*numMutations2)/geneLength
    #tempHolder = -(2*(float(numMutations1)*numMutations2)/geneLength**2)*codonRange + 2*(float(numMutations1)*numMutations2)/geneLength
    tempHolder[tempHolder < 0] = 0
    return tempHolder

# Used, maybe rename this
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

# Used
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

# Used
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

# Used
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

# Used
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

        # exactly same
        if (str(firstSeq[codonStart:codonEnd]) == str(secondSeq[codonStart:codonEnd]) 
            and str(firstSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd])): continue
        # same (though differing from outgroup) 
        elif str(firstSeq[codonStart:codonEnd]) == str(secondSeq[codonStart:codonEnd]): continue
        # Difference, but not polarizable 
        elif (str(firstSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd]) and
              str(secondSeq[codonStart:codonEnd]) != str(outSeq[codonStart:codonEnd])):
            try:
                if translationTable[str(firstSeq[codonStart:codonEnd])] == translationTable[str(secondSeq[codonStart:codonEnd])]:
                    nonPolarizableDS = nonPolarizableDS + 1
                else:
                    nonPolarizableDN = nonPolarizableDN + 1
            except KeyError:
                translationProblems = translationProblems + 1
        # Difference in first seq
        elif str(secondSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd]):
            try:
                if translationTable[str(firstSeq[codonStart:codonEnd])] == translationTable[str(outSeq[codonStart:codonEnd])]:
                    speciesOneDS.append(index)
                    numSpeciesOneDS = numSpeciesOneDS + 1
                else:
                    speciesOneDN.append(index)
                    numSpeciesOneDN = numSpeciesOneDN + 1
            except KeyError:
                translationProblems = translationProblems + 1
        # Difference in secondSeq
        elif str(firstSeq[codonStart:codonEnd]) == str(outSeq[codonStart:codonEnd]):
            try:
                if translationTable[str(secondSeq[codonStart:codonEnd])] == translationTable[str(outSeq[codonStart:codonEnd])]:
                    speciesTwoDS.append(index)
                    numSpeciesTwoDS = numSpeciesTwoDS + 1
                else:
                    speciesTwoDN.append(index)
                    numSpeciesTwoDN = numSpeciesTwoDN + 1
            except KeyError:
                translationProblems = translationProblems + 1
        else:
            assert False
    
    return {'Ambiguous Codon':ambiguousGenotypeOrBreak,'Translation Problems':translationProblems,
            'Non-polarizable DS':nonPolarizableDS, 'Non-polarizable DN':nonPolarizableDN,
            'Species One DS Count':numSpeciesOneDS,'Species One DN Count':numSpeciesOneDN,
            'Species Two DS Count':numSpeciesTwoDS,'Species Two DN Count':numSpeciesTwoDN,
            'Species One DS':speciesOneDS,'Species One DN':speciesOneDN,
            'Species Two DS':speciesTwoDS,'Species Two DN':speciesTwoDN,
            'Sequence Length':seqLength,'Stop Codon Count':stopCodons}

# Used
def windowSmoothing(data,windowSize=5):
    if windowSize%2 != 1:
        raise ValueError('Window size must be odd')
    # Cheating below, kind of. What I really want is floor and ceiling, respectively, but this does the job.
    deltaLower = windowSize//2
    deltaUpper = windowSize//2 + 1
    return [sum(data[max([i-deltaLower,0]):min([i+deltaUpper,len(data)])])/float(len(data[max([i-deltaLower,0]):min([i+deltaUpper,len(data)])])) for i in range(0,len(data))]




