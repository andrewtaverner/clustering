'''
Created on Aug 27, 2018

@author: PiTav

Some new code to use the new clustering functions that are oop-ish
'''

# Remember to rerun for Dsim, Dsec, Dmel combo

import pickle as cPickle
import time
import numpy
import glob
import functionsForClustering as clust
from itertools import groupby
from itertools import combinations
from Bio import SeqIO
from clusteringClasses import *

fileReadName = 'data/drosophila/aligned_PRANK/FBpp*_aligned.fasta*'

clusteringLength = 500

# Remember when you're setting Dyak to be one of the species, you have to pick Dyak_NY73 or Dyak_Tai18E2
speciesOneName = 'Dmel'
speciesTwoName = 'Dsim'
speciesOutName = 'Dyak_Tai18E2'

fileList = glob.glob(fileReadName)

clusteringHolder = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)

geneNameList = []
speciesOneClustering = []
speciesTwoClustering = []

for a,geneName in enumerate(fileList):
    print(a)
    
    sequenceHolder = SeqIO.to_dict(SeqIO.parse(open(geneName),'fasta'))
    
    currSeqName = dict([(x,y) for y in sequenceHolder.keys() for x in [speciesOneName,speciesTwoName,speciesOutName] if x in y])
    
    if len(currSeqName) < 3:
        continue
    
    niceGeneName = geneName.split("/")[3].split('_')[0]
    
    tempGeneHolder = protein(niceGeneName,currSeqName[speciesOneName],currSeqName[speciesTwoName],currSeqName[speciesOutName],sequenceHolder)
    
    if tempGeneHolder.filterAlignment() == -1:
        continue
    
    tempGeneHolder.processGene()
    
    geneNameList.append(niceGeneName)
    
    clusteringHolder.parseProteinObject(tempGeneHolder)


saveName = 'results_oop_code/%s_%s_%s_clustering.p' % (speciesOneName,speciesTwoName,speciesOutName)

cPickle.dump(clusteringHolder,open(saveName,'wb'))
