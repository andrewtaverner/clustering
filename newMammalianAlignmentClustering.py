'''
Created on Aug 27, 2018

@author: PiTav

Some new code to cluster specific species from the 20 species alignments (mammal)
'''

from Bio import SeqIO
from clusteringClasses import *
import cPickle
import numpy

allClusteringGenes = cPickle.load(open('data/20_mammalian_alignment/genesToCluster.p'))
goClusteringDict = cPickle.load(open('data/20_mammalian_alignment/HumanLookupDict.p'))
genesToProcess = set(list(allClusteringGenes) + goClusteringDict.keys())
tempGoInfo = cPickle.load(open('data/20_mammalian_alignment/GO_dictionary_CC.p'))
listOfGO = numpy.unique(tempGoInfo.values())

def clusterVertebrates(speciesOneName, speciesTwoName, speciesOutName):
    # Some pre-processing of alignment file
    sequenceIterator = SeqIO.parse("data/20_mammalian_alignment/filteredRefGene.exonNuc.fa","fasta")
    
    relevantInd = []
    
    for seqInd,currSeq in enumerate(sequenceIterator):
        if seqInd == 20: break
        if speciesOneName in currSeq.id or speciesTwoName in currSeq.id or speciesOutName in currSeq.id:
            relevantInd.append(seqInd)
    
    if len(relevantInd) != 3:
        print("One or more of the species is not in the alignment") 
        assert len(relevantInd) == 3
    
    sequenceIterator = SeqIO.parse("data/20_mammalian_alignment/filteredRefGene.exonNuc.fa","fasta")
    
    clusteringHolder = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)
    goClusteringHolder = dict()
    goClusteringHolder['all'] = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)
    for currGO in listOfGO:
        goClusteringHolder[currGO] = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)
    
    numberOfGenes = 0
    genesProcessed = set()
    genesRemoved = set()
    extraGenes = []
    geneName = ""
    totalExonCounter = 0
    tempSeqHolder = dict()
    exonsPerSpecies = 0
    totalNumExons = 0
    totalNumGenes = len(genesToProcess)
    
    for seqInd, currSeq in enumerate(sequenceIterator):
        # We always get information about the gene from the human sequence:
        if seqInd % 20 == 0 and geneName == "":
            tempHolder = currSeq.id.split("_")
            geneName = '_'.join(tempHolder[0:2])
            if geneName not in genesToProcess:
                geneName = ""
            else:
                exonsPerSpecies = int(tempHolder[4])
                totalNumExons = exonsPerSpecies * 3
                tempSeqHolder = dict()
                totalExonCounter = 0
        if seqInd % 20 == 0 and geneName != "":
            tempHolder = currSeq.id.split("_")
            assert "_".join(tempHolder[0:2]) == geneName
        if geneName != "" and seqInd % 20 in relevantInd:
            tempSeqHolder[currSeq.id] = currSeq.seq
            totalExonCounter += 1
        if geneName != "" and totalExonCounter == totalNumExons:
            geneSeqHolder = {speciesOneName:[], speciesTwoName:[], speciesOutName:[]}
            for currSpec in (speciesOneName,speciesTwoName,speciesOutName):
                for currExon in range(exonsPerSpecies):
                    geneSeqHolder[currSpec].append(str(tempSeqHolder['%s_%s_%s_%s' % (geneName,currSpec,currExon+1,exonsPerSpecies)]))
            tempGeneHolder = protein.proteinExons(geneName,speciesOneName,speciesTwoName,speciesOutName,geneSeqHolder)
            if tempGeneHolder.filterAlignment() == -1:
                genesRemoved.add(geneName)
                geneName = ""
                continue
            
            if tempGeneHolder in genesProcessed:
                extraGenes.append(geneName)
                geneName = ""
                continue
            
            print("Clustering gene " + str(numberOfGenes + 1) + " of " + str(totalNumGenes))
            numberOfGenes += 1
            
            genesProcessed.add(geneName)
            tempGeneHolder.processGene()
            tempGeneHolder.processGeneIntrons()
            if geneName in allClusteringGenes:
                clusteringHolder.parseProteinObject(tempGeneHolder,parseIntronInfo = True)
            if geneName in goClusteringDict.keys():
                goClusteringHolder['all'].parseProteinObject(tempGeneHolder)
                for currGO in goClusteringDict[geneName]:
                    goClusteringHolder[currGO].parseProteinObject(tempGeneHolder)
            geneName = ""
    
    temp = {'DNDN Simulation':clusteringHolder.nonPolarizedDNDNIntronSimulation,
            'DSDS Simulation':clusteringHolder.nonPolarizedDSDSIntronSimulation}
    clusteringHolder.nonPolarizedDNDNIntronSimulation = None
    clusteringHolder.nonPolarizedDSDSIntronSimulation = None
    cPickle.dump({'Clustering Holder':clusteringHolder,'Genes Processed':genesProcessed,
                  'Genes Removed':genesRemoved,'Extra Genes':extraGenes}
                 ,open("results/" + "_".join([speciesOneName,speciesTwoName,speciesOutName])+"_mammalian.p",'w'))
    
    cPickle.dump(temp, open("results/" + "_".join([speciesOneName,speciesTwoName,speciesOutName])+"_mammalian_intronSimulation.p",'w'))
    saveName = 'results/%s_%s_%s_GO_clustering_mammalian.p' % (speciesOneName,speciesTwoName,speciesOutName)
    cPickle.dump(goClusteringHolder,open(saveName,'w'))

listToCluster = [
                 ('hg38','ponAbe2','nomLeu3'),
                 #('hg38','panTro4','gorGor3'),
                 #('hg38','gorGor3','ponAbe2'),
                 #('hg38','ponAbe2','papAnu2'), # This is not the outgroup that is closest to the pair -- it is the original clustering triplet for primate
                 #('hg38','macFas5','calJac3'),
                 #('hg38','calJac3','otoGar3'),
                 #('hg38','otoGar3','tupChi1'),
                 #('rheMac3','macFas5','papAnu2'),
                 #('macFas5','papAnu2','chlSab2'),
                 #('calJac3','saiBol1','hg38'),
                 #('mm10','rn6','mesAur1'),
                 #('mesAur1','criGri1','micOch1'),
                 #('mm10','mesAur1','jacJac1'),
                 #('mm10','jacJac1','speTri2'),
                 #('oryCun2','ochPri3','mm10'), # This one appears to be broken
                 #('chiLan1','octDeg1','cavPor3'),
                 #('chiLan1','cavPor3','hetGla2'),
                 #('chiLan1','hetGla2','mm10'),
                 #('oviAri3','capHir1','bosTau8'),
                 #('oviAri3','bosTau8','panHod1'),
                 #('turTru2','orcOrc1','bosTau8')
                ]

listAlreadyClustered = set()

for currComp in listToCluster:
    if currComp not in listAlreadyClustered:
        #try:
            print(" ".join(currComp) + "\n")
            clusterVertebrates(currComp[0],currComp[1],currComp[2])
        #except Exception as e:
        #    print(e)
        #    print("\n")
    listAlreadyClustered.add(currComp)
