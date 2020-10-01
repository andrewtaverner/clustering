'''
Created on Aug 27, 2018

@author: PiTav

Yeast clustering
Some new code to use the new clustering functions that are oop-ish
'''

import pickle as cPickle
import glob
from Bio import SeqIO
from clusteringClasses import *
from collections import Counter

fileReadName = 'data/yeast/coding/*codon.mfa'

clusteringLength = 500

speciesOneName = 'Scer'
speciesTwoName = 'Spar'
speciesOutName = 'Smik'

clusteringHolder = genomeWideClustering(speciesOneName,speciesTwoName,speciesOutName)

fileList = glob.glob(fileReadName)
fileList = [x for x in fileList if 'NOSGD' not in x] # remove alignments if they don't have a SGB

geneCounts = Counter([x.split('_')[1] for x in fileList])
fileList = [x for x in fileList if geneCounts[x.split('_')[1]] == 1] # remove any proteins that exist more than once

for a,geneName in enumerate(fileList):
    print(a)
    
    sequenceHolder = SeqIO.to_dict(SeqIO.parse(open(geneName),'fasta'))
    
    currSeqName = dict([(x,y) for y in sequenceHolder.keys() for x in [speciesOneName,speciesTwoName,speciesOutName] if x in y])
    
    if len(currSeqName) != 3:
        continue
    
    niceGeneName = geneName.split("/")[3].split('_')[1]
    
    tempGeneHolder = protein(niceGeneName,currSeqName[speciesOneName],currSeqName[speciesTwoName],currSeqName[speciesOutName],sequenceHolder)
    
    if tempGeneHolder.filterAlignment() == -1:
        continue
    
    tempGeneHolder.processGene()
    
    clusteringHolder.parseProteinObject(tempGeneHolder)
    
cPickle.dump(clusteringHolder,open('results_oop_code/' + speciesOneName + '_' + speciesTwoName + '_' + speciesOutName + '_clustering.p','wb'))







