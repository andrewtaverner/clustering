from clusteringClasses import *
import functionsForClustering as clust
import argparse
from Bio import SeqIO
from collections import Counter


def outputStuff(listOfSites, refSeqHolder, speciesOneName, speciesTwoName, popOneSeqHolder, popTwoSeqHolder, passedSites):
    for currSite in listOfSites:
        print('Site: {} -- {}'.format(currSite, currSite in passedSites))
        [print('Ref {}: {}'.format(x,str(Seq.Seq(''.join([x if x != "-" else "N" for x in refSeqHolder[x][currSite:(currSite+3)]])).translate()))) for x in refSeqHolder]
        print(f'Population {speciesOneName}: ', end = "")
        tempHolder = Counter([str(x[currSite:(currSite + 3)].seq.translate()) for x in popOneSeqHolder.values()])
        print(', '.join(['{}:{}'.format(x,tempHolder[x]) for x in tempHolder]))
        print(f'Population {speciesTwoName}: ', end = "")
        tempHolder = Counter([str(x[currSite:(currSite + 3)].seq.translate()) for x in popTwoSeqHolder.values()])
        print(', '.join(['{}:{}'.format(x,tempHolder[x]) for x in tempHolder]) + '\n')
        

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--speciesOneName', default = 'speciesOne')
    parser.add_argument('--speciesTwoName', default = 'speciesTwo')
    parser.add_argument('--speciesOutName', default = 'speciesOut')
    parser.add_argument('--geneName', default = 'geneOfInterest')
    parser.add_argument('--fastaFile', default = None)
    parser.add_argument('--speciesOnePopFasta', default = None)
    parser.add_argument('--speciesTwoPopFasta', default = None)
    parser.add_argument('--noHeader', action = 'store_true')
    parser.add_argument('--speciesOneVariants', default = None)
    parser.add_argument('--speciesTwoVariants', default = None)
    parser.add_argument('--sequenceLength', default = None, type = int, help = "Sequence length in nucleotides")
    parser.add_argument('--debug', action = 'store_true')
    
    args = parser.parse_args()
    
    if args.fastaFile is not None:
        tempFastaFile = SeqIO.to_dict(SeqIO.parse(args.fastaFile,'fasta'))
        refFastaFile = dict([(x, tempFastaFile[x]) for x in tempFastaFile if x in [args.speciesOneName, args.speciesTwoName, args.speciesOutName]])
        
        tempGeneHolder = protein(args.geneName, args.speciesOneName, args.speciesTwoName, args.speciesOutName, refFastaFile)
        
        if tempGeneHolder.filterAlignment() == -1:
            print("Alignmnent didn't pass filters -- Calculating anyway")
            tempGeneHolder = protein(args.geneName, args.speciesOneName, args.speciesTwoName, args.speciesOutName, refFastaFile)
        
        #import ipdb; ipdb.set_trace()
        
        if args.debug:
            tempGeneHolder.processGene()
            speciesOneList = tempGeneHolder.speciesOneDN
            speciesTwoList = tempGeneHolder.speciesTwoDN
        
        if args.speciesOnePopFasta is not None:
            tempFastaFile1 = SeqIO.to_dict(SeqIO.parse(args.speciesOnePopFasta,'fasta'))
            tempGeneHolder.addPopulationData(1, tempFastaFile1)
        else:
            tempFastaFile1 = dict()
        if args.speciesTwoPopFasta is not None:
            tempFastaFile2 = SeqIO.to_dict(SeqIO.parse(args.speciesTwoPopFasta,'fasta'))
            tempGeneHolder.addPopulationData(2, tempFastaFile2)
        else:
            tempFastaFile2 = dict()
        
        result = tempGeneHolder.geneClusteringFixedDiff()
    
    else:
        
        fakeFastaFile = {args.speciesOneName:"A"*args.sequenceLength,
                         args.speciesTwoName:"A"*args.sequenceLength,
                         args.speciesOutName:"A"*args.sequenceLength}
        
        tempGeneHolder = protein(args.geneName, args.speciesOneName, args.speciesTwoName, args.speciesOutName, fakeFastaFile)
        
        tempGeneHolder.seqLength = int(args.sequenceLength)
        
        if args.speciesOneVariants is not None:
            speciesOneVariants = sorted([int(x)//3*3 for x in args.speciesOneVariants.split(',')])
        else:
            speciesOneVariants = []
        if args.speciesTwoVariants is not None:
            speciesTwoVariants = sorted([int(x)//3*3 for x in args.speciesTwoVariants.split(',')])
        else:
            speciesTwoVariants = []
        
        tempGeneHolder.setPolarizedDivergenceList(speciesOneVariants,[],
                                                  speciesTwoVariants,[])
        
        #import ipdb; ipdb.set_trace()
        
        result = tempGeneHolder.geneClusteringFixedDiff(skipVariantFinding = True)
    
    if not args.noHeader:
        print(f"{args.speciesOneName}\t{args.speciesTwoName}")
    print("{}\t{}".format(*result))
    
    if args.debug:
        print("Number of subs found in each:\n{}\t{}\n".format(len(tempGeneHolder.speciesOneDN),len(tempGeneHolder.speciesTwoDN)))
        print('{}_positions: {}'.format(args.speciesOneName, ','.join([str(x) for x in tempGeneHolder.speciesOneDN])))
        outputStuff(speciesOneList, refFastaFile, args.speciesOneName, args.speciesTwoName ,tempFastaFile1, tempFastaFile2, tempGeneHolder.speciesOneDN)
        print('{}_positions: {}'.format(args.speciesTwoName, ','.join([str(x) for x in tempGeneHolder.speciesTwoDN])))
        outputStuff(speciesTwoList, refFastaFile, args.speciesOneName, args.speciesTwoName, tempFastaFile1, tempFastaFile2, tempGeneHolder.speciesTwoDN)

    
if __name__ == '__main__':
    main()
