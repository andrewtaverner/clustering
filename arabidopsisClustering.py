'''
Created on Sept 11, 2019

@author: PiTav

Some new code to use the new clustering functions that are oop-ish
'''

import pickle as cPickle
import time
import numpy
import glob
import functionsForClustering as clust
from itertools import groupby
from itertools import combinations
from Bio import SeqIO
from clusteringClasses import *

fileReadName = 'data/arabidopsis/CrAtAl.prank.sorted.fasta/Carubv*fasta*'

clusteringLength = 500

speciesOneName = 'Atha'
speciesTwoName = 'Alyr'
speciesOutName = 'Crub'

fileList = glob.glob(fileReadName)

clusteringHolder = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)

geneNameList = []
speciesOneClustering = []
speciesTwoClustering = []

for a,geneName in enumerate(fileList):
    print(a)
    
    sequenceHolder = SeqIO.to_dict(SeqIO.parse(open(geneName),'fasta'))
    
    tempHolder = set(sequenceHolder.keys())
    
    tempSpeciesTwoName = [x for x in tempHolder if x[0:2] == 'AT']
    tempHolder -= set(tempSpeciesTwoName)
    tempSpeciesOutName = [x for x in tempHolder if x[0:6] == 'Carubv']
    tempHolder -= set(tempSpeciesOutName)
    tempSpeciesOneName = tempHolder.pop()
    
    currSeqName = {speciesOneName:tempSpeciesOneName,speciesTwoName:tempSpeciesTwoName[0],speciesOutName:tempSpeciesOutName[0]}
    
    if len(currSeqName) < 3:
        continue
    
    niceGeneName = geneName.split("/")[3].split('_')[0]
    
    tempGeneHolder = protein(nice ame,currSeqName[speciesOneName],currSeqName[speciesTwoName],currSeqName[speciesOutName],sequenceHolder)
    
    if tempGeneHolder.filterAlignment() == -1:
        continue
    
    tempGeneHolder.processGene()
    
    geneNameList.append(niceGeneName)
    
    clusteringHolder.parseProteinObject(tempGeneHolder)
    
saveName = 'results_oop_code/%s_%s_%s_clustering.p' % (speciesOneName,speciesTwoName,speciesOutName)

cPickle.dump(clusteringHolder,open(saveName,'wb'))
