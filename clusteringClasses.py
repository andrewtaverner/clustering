'''
Created on Aug 26, 2017

@author: PiTav

This is the new version of this script
in the finalizedForPub folder
''' 

import functionsForClustering as clust
import numpy
from Bio import Seq
import matplotlib.pyplot as plt
from scipy import stats

class protein(object):
    
    def __init__(self, geneName, speciesOneName, speciesTwoName, speciesOutName, seqHolder, chrPos=None):
        self.name = geneName
        self.chrPos = chrPos
        self.speciesOneName = speciesOneName
        self.speciesTwoName = speciesTwoName
        self.speciesOutName = speciesOutName
        try:
            self.firstSeq = seqHolder[speciesOneName].seq
            self.secondSeq = seqHolder[speciesTwoName].seq
            self.outSeq = seqHolder[speciesOutName].seq
        except AttributeError: 
            try: 
                self.firstSeq = Seq.Seq(seqHolder[speciesOneName])
                self.secondSeq = Seq.Seq(seqHolder[speciesTwoName])
                self.outSeq = Seq.Seq(seqHolder[speciesOutName])
            except TypeError:
                self.firstSeq = seqHolder[speciesOneName]
                self.secondSeq = seqHolder[speciesTwoName]
                self.outSeq = seqHolder[speciesOutName]
        self.seqLength = 0
        self.clusteringHolder = dict()
        self.clusteringIntronHolder = dict()
        self.polarizedClusteringHolder = dict()
        self.propertyClusteringHolder = dict()
        self.DNLists = None 
        self.DSLists = None 
        self.speciesOneDN = None
        self.speciesTwoDN = None
        self.speciesOneDS = None
        self.speciesTwoDS = None
        self.propertyClusteringList = dict() # possibly should remove these variables
        self.intronPositions = None
        self.reasonForFailure = None
        self.speciesOnePopHolder = None
        self.speciesTwoPopHolder = None
        self.speciesOutPopHolder = None
    
    # For this method, the seqHolder should contain a list of strings that are exons in the correct order!
    # Must be a list, even if there is only a single exon, otherwise, disaster!
    # This only really works if the intron positions are conserved and information is taken from species one
    @classmethod
    def proteinExons(cls, geneName, speciesOneName, speciesTwoName, speciesOutName, seqHolder, chrPos=None):
        # The last one is the end of the gene, so remove:
        intronPositions = numpy.cumsum([len(x) for x in seqHolder[speciesOneName]])[:-1] 
        # Use this to get the introns that occur nicely between codons to be decimal numbers:
        betweenCodons = (intronPositions%3 == 0)*0.5
        intronPositions//=3
        intronPositions = intronPositions.astype(float)
        intronPositions -= betweenCodons
        tempHolder = dict()
        try:
            tempHolder[speciesOneName] = ''.join([str(x.seq) for x in seqHolder[speciesOneName]])
            tempHolder[speciesTwoName] = ''.join([str(x.seq) for x in seqHolder[speciesTwoName]])
            tempHolder[speciesOutName] = ''.join([str(x.seq) for x in seqHolder[speciesOutName]])
        except AttributeError:
            tempHolder[speciesOneName] = ''.join([str(x) for x in seqHolder[speciesOneName]])
            tempHolder[speciesTwoName] = ''.join([str(x) for x in seqHolder[speciesTwoName]])
            tempHolder[speciesOutName] = ''.join([str(x) for x in seqHolder[speciesOutName]])
        
        tempHolder = cls(geneName,speciesOneName,speciesTwoName,speciesOutName,tempHolder,chrPos)
        tempHolder.intronPositions = intronPositions
        
        return tempHolder
            
    def addIntronPos(self,intronPositions):
        # This method assumes that the intron positions are in nucleotide space and indexed like python do
        betweenCodons = (intronPositions%3 == 0)*0.5
        intronPositions//=3
        intronPositions = intronPositions.astype(float)
        intronPositions -= betweenCodons
        self.intronPositions = intronPositions
        
    def setDivergenceList(self,DNList,DSList):
        self.DNLists = DNList
        self.DSLists = DSList
    
    def setPolarizedDivergenceList(self,DNListOne,DSListOne,DNListTwo,DSListTwo):
        self.speciesOneDN = DNListOne
        self.speciesTwoDN = DNListTwo
        self.speciesOneDS = DSListOne
        self.speciesTwoDS = DSListTwo
    
    # This is a wrapper function from functionsForClustering
    def filterAlignment(self, lenMultipleOfThree = True, minAlnLen = 300, 
                        eachGapMultipleOfThree = True, totalGapLengthMultipleOfThree = True,
                        gapContentThresh = .2, DNContentThresh = .2, DSContentThresh = 1, 
                        prematureStopCodon = True):
        result = clust.filterAlignment(self.firstSeq,self.secondSeq,self.outSeq,
                                       lenMultipleOfThree, minAlnLen, eachGapMultipleOfThree,
                                       totalGapLengthMultipleOfThree, gapContentThresh, 
                                       DNContentThresh, DSContentThresh, prematureStopCodon)
        if result[0] == -1:
            self.reasonForFailure = result[1]
            self.firstSeq = None
            self.secondSeq = None
            self.outSeq = None
            return -1
        else: 
            return 1
    
    def processGene(self):
        result = clust.variantsNotPolarized(self.firstSeq,self.secondSeq)
        self.setDivergenceList(result['DN List'],result['DS List'])
        self.seqLength = len(self.firstSeq)
        temp_data = clust.clusteringSameMutation(self.DNLists)
        temp_norm = clust.analyticalNormSameMut(len(self.DNLists),self.seqLength/3)
        self.clusteringHolder['DNDN'] = clusteringObject(temp_data,temp_norm)
        
        temp_data = clust.clusteringSameMutation(self.DSLists)
        temp_norm = clust.analyticalNormSameMut(len(self.DSLists),self.seqLength/3)
        self.clusteringHolder['DSDS'] = clusteringObject(temp_data,temp_norm)
        
        #temp_data = clust.clusteringDifferentMutation(self.DSLists,self.DNLists)
        #temp_norm = clust.analyticalNormDiffMut(len(self.DSLists),len(self.DNLists),self.seqLength/3)
        #self.clusteringHolder['DNDS'] = clusteringObject(temp_data,temp_norm)
        
        result = clust.variantsByParsimony(self.firstSeq,self.secondSeq,self.outSeq)
        self.setPolarizedDivergenceList(result['Species One DN'],result['Species One DS'],
                                        result['Species Two DN'],result['Species Two DS'])
        
        temp_data = clust.clusteringSameMutation(self.speciesOneDN)
        temp_norm = clust.analyticalNormSameMut(len(self.speciesOneDN),self.seqLength/3)
        self.polarizedClusteringHolder['DNDN1'] = clusteringObject(temp_data,temp_norm)
        
        temp_data = clust.clusteringSameMutation(self.speciesTwoDN)
        temp_norm = clust.analyticalNormSameMut(len(self.speciesTwoDN),self.seqLength/3)
        self.polarizedClusteringHolder['DNDN2'] = clusteringObject(temp_data,temp_norm)
        
        temp_data = clust.clusteringDifferentMutation(self.speciesOneDN,self.speciesTwoDN)
        temp_norm = clust.analyticalNormDiffMut(len(self.speciesOneDN),len(self.speciesTwoDN),self.seqLength/3)
        self.polarizedClusteringHolder['DNDNcross'] = clusteringObject(temp_data,temp_norm)
        
        chargeTable = {'A': 0, 'C': 0, 'E': -1, 'D': -1, 'G': 0, 'F': 0, 'I': 0, 'H': 1, 
                       'K': 1, 'M': 0, 'L': 0, 'N': 0, 'Q': 0, 'P': 0, 'S': 0, 'R': 1, 
                       'T': 0, 'W': 0, 'V': 0, 'Y': 0}
        polarityTable = {'A': 1, 'C': 1, 'E': 0, 'D': 0, 'G': 0, 'F': 1, 'I': 1, 'H': -1, 
                         'K': -1, 'M': 1, 'L': 1, 'N': -1, 'Q': 0, 'P': -1, 'S': 0, 
                         'R': -1, 'T': 0, 'W': 1, 'V': 1, 'Y': 1}
        sizeTable = {'A': 67, 'C': 86, 'E': 109, 'D': 91, 'G': 48, 'F': 135, 'I': 124, 
                     'H': 118, 'K': 135, 'M': 124, 'L': 124, 'N': 96, 'Q': 114, 'P': 90, 
                     'S': 73, 'R': 148, 'T': 93, 'W': 163, 'V': 105, 'Y': 141}
        
        result = clust.propertyClusteringNew(self.firstSeq,self.secondSeq,self.outSeq,
            self.speciesOneDN,self.speciesTwoDN,chargeTable)
        self.propertyClusteringHolder['Charge'] = dict([(currKey,result[currKey]) for currKey in result.keys() if ("Inc" not in currKey and "Dec" not in currKey)])
        result = clust.propertyClusteringNew(self.firstSeq,self.secondSeq,self.outSeq,
            self.speciesOneDN,self.speciesTwoDN,polarityTable)
        self.propertyClusteringHolder['Polarity'] = dict([(currKey,result[currKey]) for currKey in result.keys() if ("Inc" not in currKey and "Dec" not in currKey)])
        result = clust.propertyClusteringNew(self.firstSeq,self.secondSeq,self.outSeq,
            self.speciesOneDN,self.speciesTwoDN,sizeTable)
        self.propertyClusteringHolder['Size'] = dict([(currKey,result[currKey]) for currKey in result.keys() if ("Inc" not in currKey and "Dec" not in currKey)])
    
    def processGeneIntrons(self,calculatePolarized = False,num_simulations = 100):
        if not calculatePolarized:
            if self.DNLists == None:
                result = clust.variantsNotPolarized(self.firstSeq,self.secondSeq)
                self.setDivergenceList(result['DN List'],result['DS List'])
                self.seqLength = len(self.firstSeq)
            tempHolder = clust.clusteringSameMutationIntrons(self.DNLists,self.intronPositions,self.seqLength/3,num_simulations)
            self.clusteringIntronHolder['DNDN'] = clusteringObject(tempHolder['data'],tempHolder['norm'])
            tempHolder = clust.clusteringSameMutationIntrons(self.DSLists,self.intronPositions,self.seqLength/3,num_simulations)
            self.clusteringIntronHolder['DSDS'] = clusteringObject(tempHolder['data'],tempHolder['norm'])
        #if calculatePolarized:
        #    self.polarizedIntronClusteringHolder = dict()
        #    try:
        #        self.polarizedClusteringHolder
        #    except AttributeError:
        #        result = variantsByParsimony(self.firstSeq,self.secondSeq,self.outSeq)
        #        self.setPolarizedDivergenceList(result['Species One DN'],result['Species One DS'],
        #                                        result['Species Two DN'],result['Species Two DS'])
        #    tempHolder = clust.clusteringSameMutationIntrons(self.speciesOneDN,self.intronPositions,self.seqLength/3,num_simulations)
        #    self.polarizedIntronClusteringHolder['DNDN1'] = clusteringObject(tempHolder['data'],tempHolder['norm'])
        #    tempHolder = clust.clusteringSameMutationIntrons(self.speciesTwoDN,self.intronPositions,self.seqLength/3,num_simulations)
        #    self.polarizedIntronClusteringHolder['DNDN2'] = clusteringObject(tempHolder['data'],tempHolder['norm'])

    def geneClusteringIntrons(self,length=50):
        speciesOneClustering = sum(self.polarizedIntronClusteringHolder['DNDN1'].data[1:length+1])
        speciesOneNorm = sum(self.polarizedIntronClusteringHolder['DNDN1'].norm[1:length+1])
        speciesTwoClustering = sum(self.polarizedIntronClusteringHolder['DNDN2'].data[1:length+1])
        speciesTwoNorm = sum(self.polarizedIntronClusteringHolder['DNDN2'].norm[1:length+1])
        betweenSpeciesClustering = sum(self.polarizedIntronClusteringHolder['DNDNcross'].data[1:length+1])
        betweenSpeciesNorm = sum(self.polarizedIntronClusteringHolder['DNDNcross'].norm[1:length+1])
        
        if speciesOneNorm == 0: speciesOneNorm = 1
        if speciesTwoNorm == 0: speciesTwoNorm = 1
        if betweenSpeciesNorm == 0: betweenSpeciesNorm = 1
        
        return([speciesOneClustering/speciesOneNorm - betweenSpeciesClustering/betweenSpeciesNorm,
            speciesTwoClustering/speciesTwoNorm - betweenSpeciesClustering/betweenSpeciesNorm])

class clusteringObject(object):
    
    def __init__(self,data,norm):
        self.data = data
        self.norm = norm

class genomeWideClustering(object):
    
    def __init__(self,speciesOneName,speciesTwoName,speciesOutName,clusteringLength = 500):
        
        #self.totalGeneLength = 0
        #self.totalDNMutations = 0
        #self.totalDSMutations = 0
        #self.totalMutations = 0
        self.geneLengthList = []
        self.DNMutationList = []
        self.DSMutationList = []
        self.speciesOneDNMutationList = []
        self.speciesTwoDNMutationList = []
        self.speciesOneDSMutationList = []
        self.speciesTwoDSMutationList = []
        
        self.numberOfGenes = 0
        
        self.speciesOneName = speciesOneName
        self.speciesTwoName = speciesTwoName
        self.speciesOutName = speciesOutName
        
        self.nonPolarizedDNDNClustering = numpy.zeros(clusteringLength + 1)
        #self.nonPolarizedDNDSClustering = numpy.zeros(clusteringLength + 1)
        self.nonPolarizedDSDSClustering = numpy.zeros(clusteringLength + 1)
        self.speciesOneDNDNClustering = numpy.zeros(clusteringLength + 1)
        self.speciesTwoDNDNClustering = numpy.zeros(clusteringLength + 1)
        self.btwnDNDNClustering = numpy.zeros(clusteringLength + 1)

        self.nonPolarizedDNDNClusteringNorm = numpy.zeros(clusteringLength + 1)
        #self.nonPolarizedDNDSClusteringNorm = numpy.zeros(clusteringLength + 1)
        self.nonPolarizedDSDSClusteringNorm = numpy.zeros(clusteringLength + 1)
        self.speciesOneDNDNClusteringNorm = numpy.zeros(clusteringLength + 1)
        self.speciesTwoDNDNClusteringNorm = numpy.zeros(clusteringLength + 1)
        self.btwnDNDNClusteringNorm = numpy.zeros(clusteringLength + 1)

        self.speciesOneChargeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesOneChargeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.speciesOneSizeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesOneSizeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.speciesOnePolarityClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesOnePolarityClusteringRein = numpy.zeros(clusteringLength + 1)
        
        self.speciesTwoChargeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesTwoChargeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.speciesTwoSizeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesTwoSizeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.speciesTwoPolarityClusteringComp = numpy.zeros(clusteringLength + 1)
        self.speciesTwoPolarityClusteringRein = numpy.zeros(clusteringLength + 1)
                
        self.btwnChargeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.btwnChargeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.btwnSizeClusteringComp = numpy.zeros(clusteringLength + 1)
        self.btwnSizeClusteringRein = numpy.zeros(clusteringLength + 1)
        self.btwnPolarityClusteringComp = numpy.zeros(clusteringLength + 1)
        self.btwnPolarityClusteringRein = numpy.zeros(clusteringLength + 1)
        
        self.nonPolarizedDNDNIntronClustering = numpy.zeros(clusteringLength + 1)
        #self.nonPolarizedDNDSIntronClustering = numpy.zeros(clusteringLength + 1)
        self.nonPolarizedDSDSIntronClustering = numpy.zeros(clusteringLength + 1)
        self.nonPolarizedDNDNIntronClusteringNorm = numpy.zeros(clusteringLength + 1)
        #self.nonPolarizedDNDSIntronClusteringNorm = numpy.zeros(clusteringLength + 1)
        self.nonPolarizedDSDSIntronClusteringNorm = numpy.zeros(clusteringLength + 1)
        
        self.nonPolarizedDNDNIntronSimulation = []
        self.nonPolarizedDSDSIntronSimulation = []
        
    def parseProteinObject(self,currProtein,parseIntronInfo = False):
        self.geneLengthList.append(currProtein.seqLength)
        self.DNMutationList.append(len(currProtein.DNLists))
        self.DSMutationList.append(len(currProtein.DSLists))
        self.speciesOneDNMutationList.append(len(currProtein.speciesOneDN))
        self.speciesOneDSMutationList.append(len(currProtein.speciesOneDS))
        self.speciesTwoDNMutationList.append(len(currProtein.speciesTwoDN))
        self.speciesTwoDSMutationList.append(len(currProtein.speciesTwoDS))
        
        self.numberOfGenes += 1
        
        self.nonPolarizedDNDNClustering += currProtein.clusteringHolder['DNDN'].data
        self.nonPolarizedDNDNClusteringNorm += currProtein.clusteringHolder['DNDN'].norm
        self.nonPolarizedDSDSClustering += currProtein.clusteringHolder['DSDS'].data
        self.nonPolarizedDSDSClusteringNorm += currProtein.clusteringHolder['DSDS'].norm
        #self.nonPolarizedDNDSClustering += currProtein.clusteringHolder['DNDS'].data
        #self.nonPolarizedDNDSClusteringNorm += currProtein.clusteringHolder['DNDS'].norm
        
        self.speciesOneDNDNClustering += currProtein.polarizedClusteringHolder['DNDN1'].data
        self.speciesOneDNDNClusteringNorm += currProtein.polarizedClusteringHolder['DNDN1'].norm
        self.speciesTwoDNDNClustering += currProtein.polarizedClusteringHolder['DNDN2'].data
        self.speciesTwoDNDNClusteringNorm += currProtein.polarizedClusteringHolder['DNDN2'].norm
        self.btwnDNDNClustering += currProtein.polarizedClusteringHolder['DNDNcross'].data
        self.btwnDNDNClusteringNorm += currProtein.polarizedClusteringHolder['DNDNcross'].norm
        
        self.speciesOneChargeClusteringComp += currProtein.propertyClusteringHolder['Charge']['Species One Comp']
        self.speciesOneChargeClusteringRein += currProtein.propertyClusteringHolder['Charge']['Species One Rein']
        self.speciesOneSizeClusteringComp += currProtein.propertyClusteringHolder['Size']['Species One Comp']
        self.speciesOneSizeClusteringRein += currProtein.propertyClusteringHolder['Size']['Species One Rein']
        self.speciesOnePolarityClusteringComp += currProtein.propertyClusteringHolder['Polarity']['Species One Comp']
        self.speciesOnePolarityClusteringRein += currProtein.propertyClusteringHolder['Polarity']['Species One Rein']
        
        self.speciesTwoChargeClusteringComp += currProtein.propertyClusteringHolder['Charge']['Species Two Comp']
        self.speciesTwoChargeClusteringRein += currProtein.propertyClusteringHolder['Charge']['Species Two Rein']
        self.speciesTwoSizeClusteringComp += currProtein.propertyClusteringHolder['Size']['Species Two Comp']
        self.speciesTwoSizeClusteringRein += currProtein.propertyClusteringHolder['Size']['Species Two Rein']
        self.speciesTwoPolarityClusteringComp += currProtein.propertyClusteringHolder['Polarity']['Species Two Comp']
        self.speciesTwoPolarityClusteringRein += currProtein.propertyClusteringHolder['Polarity']['Species Two Rein']
        
        self.btwnChargeClusteringComp += currProtein.propertyClusteringHolder['Charge']['Cross Species Comp']
        self.btwnChargeClusteringRein += currProtein.propertyClusteringHolder['Charge']['Cross Species Rein']
        self.btwnSizeClusteringComp += currProtein.propertyClusteringHolder['Size']['Cross Species Comp']
        self.btwnSizeClusteringRein += currProtein.propertyClusteringHolder['Size']['Cross Species Rein']
        self.btwnPolarityClusteringComp += currProtein.propertyClusteringHolder['Polarity']['Cross Species Comp']
        self.btwnPolarityClusteringRein += currProtein.propertyClusteringHolder['Polarity']['Cross Species Rein']
        
        if parseIntronInfo:
            self.nonPolarizedDNDNIntronClustering += currProtein.clusteringIntronHolder['DNDN'].data
            self.nonPolarizedDNDNIntronSimulation.append(currProtein.clusteringIntronHolder['DNDN'].norm)
            self.nonPolarizedDNDNIntronClusteringNorm += numpy.sum(currProtein.clusteringIntronHolder['DNDN'].norm,axis = 0)/currProtein.clusteringIntronHolder['DNDN'].norm.shape[0]
            self.nonPolarizedDSDSIntronClustering += currProtein.clusteringIntronHolder['DSDS'].data
            self.nonPolarizedDSDSIntronClusteringNorm += numpy.sum(currProtein.clusteringIntronHolder['DSDS'].norm,axis = 0)/currProtein.clusteringIntronHolder['DSDS'].norm.shape[0]
            self.nonPolarizedDSDSIntronSimulation.append(currProtein.clusteringIntronHolder['DSDS'].norm)
            #self.nonPolarizedDNDSIntronClustering += currProtein.clusteringIntronHolder['DNDS'].data
            #self.nonPolarizedDNDSIntronClusteringNorm += currProtein.clusteringIntronHolder['DNDS'].norm
            
    def calcNonPolarized(self, smooth = True , normalizeAsym = False):
        DNDN = self.nonPolarizedDNDNClustering/self.nonPolarizedDNDNClusteringNorm
        #DNDS = self.nonPolarizedDNDSClustering/self.nonPolarizedDNDSClusteringNorm
        DSDS = self.nonPolarizedDSDSClustering/self.nonPolarizedDSDSClusteringNorm
        if normalizeAsym:
            DNDN = DNDN/numpy.mean(DNDN[80:121])
            #DNDS = DNDS/numpy.mean(DNDS[80:121])
            DSDS = DSDS/numpy.mean(DSDS[80:121])
        if smooth:
            DNDN = [0] + clust.windowSmoothing(DNDN[1:])
            #DNDS = [0] + clust.windowSmoothing(DNDS[1:])
            DSDS = [0] + clust.windowSmoothing(DSDS[1:])
        return {'DNDN':DNDN,'DSDS':DSDS}
    
    def plotNonPolarized(self,plotTitle = None, showLegend = True, normalizeAsym = False, minMaxX = (0,100), plotIntron = False, saveFigure = None, ax = None, plotDNDS = True):
        if not plotIntron:
            DNDN = self.nonPolarizedDNDNClustering/self.nonPolarizedDNDNClusteringNorm
            #DNDS = self.nonPolarizedDNDSClustering/self.nonPolarizedDNDSClusteringNorm
            DSDS = self.nonPolarizedDSDSClustering/self.nonPolarizedDSDSClusteringNorm
        if plotIntron:
            DNDN = self.nonPolarizedDNDNIntronClustering/self.nonPolarizedDNDNIntronClusteringNorm
            #DNDS = self.nonPolarizedDNDSIntronClustering/self.nonPolarizedDNDSIntronClusteringNorm
            DSDS = self.nonPolarizedDSDSIntronClustering/self.nonPolarizedDSDSIntronClusteringNorm
        if normalizeAsym:
            DNDN = DNDN/numpy.mean(DNDN[80:121])
            #DNDS = DNDS/numpy.mean(DNDS[80:121])
            DSDS = DSDS/numpy.mean(DSDS[80:121])
        if plotDNDS:
            minY = min([min(clust.windowSmoothing(DNDN[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(DSDS[1:])[minMaxX[0]:minMaxX[1]])])
            maxY = max([max(clust.windowSmoothing(DNDN[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(DSDS[1:])[minMaxX[0]:minMaxX[1]])])
        else:
            minY = min([min(clust.windowSmoothing(DNDN[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(DSDS[1:])[minMaxX[0]:minMaxX[1]])])
            maxY = max([max(clust.windowSmoothing(DNDN[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(DSDS[1:])[minMaxX[0]:minMaxX[1]])])
        if ax == None:
            fig = plt.figure(figsize = (6,6))
            ax = fig.add_subplot()
        ax.plot(range(1,501),clust.windowSmoothing(DNDN[1:]),color='green',label='DNDN', linewidth = 1.5)
        ax.plot(range(1,501),clust.windowSmoothing(DSDS[1:]),color='blue',label='DSDS', linewidth = 1.5)
        ax.scatter(range(1,501),DNDN[1:],color = 'green', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
        ax.scatter(range(1,501),DSDS[1:],color = 'blue', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
        #if plotDNDS:
        #    ax.plot(range(1,501),clust.windowSmoothing(DNDS[1:]),color='orange',label='DNDS', linewidth = 1.5)
        #    ax.scatter(range(1,501),DNDS[1:],color = 'orange', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
        ax.set(xlim = (minMaxX[0],minMaxX[1]),ylim = (minY*.9,maxY*1.1))
        if showLegend:
            ax.legend()
        if plotTitle != None:
            ax.title.set_text(plotTitle)
        if saveFigure != None:
            plt.savefig(saveFigure)
            
    def calcPolarized(self, smooth = True, normalizeAsym = False):
        DNDN1 = self.speciesOneDNDNClustering/self.speciesOneDNDNClusteringNorm
        DNDN2 = self.speciesTwoDNDNClustering/self.speciesTwoDNDNClusteringNorm
        DNDNcross = self.btwnDNDNClustering/self.btwnDNDNClusteringNorm
        if normalizeAsym:
            DNDN1 = DNDN1/numpy.mean(DNDN1[80:121])
            DNDN2 = DNDN2/numpy.mean(DNDN2[80:121])
            DNDNcross = DNDNcross/numpy.mean(DNDNcross[80:121])
        if smooth:
            DNDN1 = [0] + clust.windowSmoothing(DNDN1[1:])
            DNDN2 = [0] + clust.windowSmoothing(DNDN2[1:])
            DNDNcross = [0] + clust.windowSmoothing(DNDNcross[1:])
        return {'DNDN1':DNDN1,'DNDN2':DNDN2,'DNDNcross':DNDNcross}
    
    def plotPolarized(self,plotTitle = None, showLegend = True, normalizeAsym = False, customNameOne = None, customNameTwo = None, customNameOut = None, minMaxX = (0,100), saveFigure = None, ax = None):
        if customNameOne == None:
            customNameOne = 'DNDN ' + self.speciesOneName
        if customNameTwo == None:
            customNameTwo = 'DNDN ' + self.speciesTwoName
        if customNameOut == None:
            customNameOut = 'DNDN cross'
        DNDN1 = self.speciesOneDNDNClustering/self.speciesOneDNDNClusteringNorm
        DNDN2 = self.speciesTwoDNDNClustering/self.speciesTwoDNDNClusteringNorm
        DNDNcross = self.btwnDNDNClustering/self.btwnDNDNClusteringNorm
        if normalizeAsym:
            DNDN1 = DNDN1/numpy.mean(DNDN1[80:121])
            DNDN2 = DNDN2/numpy.mean(DNDN2[80:121])
            DNDNcross = DNDNcross/numpy.mean(DNDNcross[80:121])
        minY = min([min(clust.windowSmoothing(DNDN1[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(DNDN2[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(DNDNcross[1:])[minMaxX[0]:minMaxX[1]])])
        maxY = max([max(clust.windowSmoothing(DNDN1[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(DNDN2[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(DNDNcross[1:])[minMaxX[0]:minMaxX[1]])])
        if ax == None:
            fig = plt.figure(figsize = (6,6))
            ax = fig.add_subplot()
        ax.plot(range(1,501),clust.windowSmoothing(DNDN1[1:]),color = 'green', label = customNameOne, linewidth = 1.5)
        ax.scatter(range(1,501),DNDN1[1:],color = 'green', alpha = 0.5, marker = '+',s = 18,linewidth = 1)
        ax.plot(range(1,501),clust.windowSmoothing(DNDN2[1:]),color = 'orange', label = customNameTwo, linewidth = 1.5)
        ax.scatter(range(1,501),DNDN2[1:],color = 'orange', alpha = 0.5, marker = '+',s=18,linewidth = 1)
        ax.plot(range(1,501),clust.windowSmoothing(DNDNcross[1:]),color = 'blue', label = customNameOut, linewidth = 1.5)
        ax.scatter(range(1,501),DNDNcross[1:],color = 'blue', alpha = 0.5, marker = '+',s=18,linewidth = 1)
        ax.set(xlim = (minMaxX[0],minMaxX[1]), ylim = (minY*0.9,maxY*1.1))
        if showLegend:
            ax.legend()
        if plotTitle != None:
            ax.title.set_text(plotTitle)
        if saveFigure != None:
            plt.savefig(saveFigure)
    
    def plotProperty(self,propertyType,compOrRein,plotTitle = None, showLegend = True,customNameOne = None, customNameTwo = None, customNameOut = None, minMaxY = None, minMaxX = (0,100), saveFigure = None, ax = None):
        if customNameOne == None:
            customNameOne = self.speciesOneName + " lineage"
        if customNameTwo == None:
            customNameTwo = self.speciesTwoName + " lineage"
        if customNameOut == None:
            customNameOut = 'Between lineage'
        if propertyType == "Charge":
            comp1 = self.speciesOneChargeClusteringComp
            comp2 = self.speciesTwoChargeClusteringComp
            compCross = self.btwnChargeClusteringComp
            rein1 = self.speciesOneChargeClusteringRein
            rein2 = self.speciesTwoChargeClusteringRein
            reinCross = self.btwnChargeClusteringRein
        elif propertyType == "Size":
            comp1 = self.speciesOneSizeClusteringComp
            comp2 = self.speciesTwoSizeClusteringComp
            compCross = self.btwnSizeClusteringComp
            rein1 = self.speciesOneSizeClusteringRein
            rein2 = self.speciesTwoSizeClusteringRein
            reinCross = self.btwnSizeClusteringRein
        elif propertyType == "Polarity":
            comp1 = self.speciesOnePolarityClusteringComp
            comp2 = self.speciesTwoPolarityClusteringComp
            compCross = self.btwnPolarityClusteringComp
            rein1 = self.speciesOnePolarityClusteringRein
            rein2 = self.speciesTwoPolarityClusteringRein
            reinCross = self.btwnPolarityClusteringRein
        else:
            print("Unknown type selected")
            return
        comp1 = comp1/self.speciesOneDNDNClustering
        comp2 = comp2/self.speciesTwoDNDNClustering
        compCross = compCross/self.btwnDNDNClustering
        rein1 = rein1/self.speciesOneDNDNClustering
        rein2 = rein2/self.speciesTwoDNDNClustering
        reinCross = reinCross/self.btwnDNDNClustering
        if minMaxY == None:
            minY = min([min(clust.windowSmoothing(comp1[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(comp2[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(compCross[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(rein1[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(rein2[1:])[minMaxX[0]:minMaxX[1]]),min(clust.windowSmoothing(reinCross[1:])[minMaxX[0]:minMaxX[1]])])
            maxY = max([max(clust.windowSmoothing(comp1[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(comp2[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(compCross[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(rein1[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(rein2[1:])[minMaxX[0]:minMaxX[1]]),max(clust.windowSmoothing(reinCross[1:])[minMaxX[0]:minMaxX[1]])])
            yRange = maxY - minY
            minY = minY - yRange * 0.4
            maxY = maxY + yRange * 0.4
        else:
            minY = minMaxY[0]
            maxY = minMaxY[1]
        if ax == None:
            fig = plt.figure(figsize = (6,6))
            ax = fig.add_subplot()
        if compOrRein == "Comp":
            ax.plot(range(1,501),clust.windowSmoothing(comp1[1:]),color = 'green', label = customNameOne, linewidth = 1.5)
            ax.scatter(range(1,501),comp1[1:],color = 'green', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
            ax.plot(range(1,501),clust.windowSmoothing(comp2[1:]),color = 'orange', label = customNameTwo, linewidth = 1.5)
            ax.scatter(range(1,501),comp2[1:],color = 'orange', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
            ax.plot(range(1,501),clust.windowSmoothing(compCross[1:]),color = 'blue', label = customNameOut, linewidth = 1.5)
            ax.scatter(range(1,501),compCross[1:],color = 'blue', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
        elif compOrRein == "Rein":
            ax.plot(range(1,501),clust.windowSmoothing(rein1[1:]),color = 'green', label = customNameOne, linewidth = 1.5)
            ax.scatter(range(1,501),rein1[1:],color = 'green', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
            ax.plot(range(1,501),clust.windowSmoothing(rein2[1:]),color = 'orange', label = customNameTwo, linewidth = 1.5)
            ax.scatter(range(1,501),rein2[1:],color = 'orange', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
            ax.plot(range(1,501),clust.windowSmoothing(reinCross[1:]),color = 'blue', label = customNameOut, linewidth = 1.5)
            ax.scatter(range(1,501),reinCross[1:],color = 'blue', alpha = 0.5, marker = '+', s = 18, linewidth = 1)
        ax.set(xlim = (minMaxX[0],minMaxX[1]), ylim = (minY,maxY))
        if showLegend:
            ax.legend()
        if plotTitle != None:
            ax.title.set_text(plotTitle)
        if saveFigure != None:
            plt.savefig(saveFigure)
    
    def helperGetPropertySpecified(self, whichSpecies = None, propertyType = None, compOrRein = None):
        # whichSpecies == 3 indicates between
        if whichSpecies == 1:
            if propertyType == 'Charge' and compOrRein == 'Comp':
                clustObs = self.speciesOneChargeClusteringComp
            elif propertyType == 'Size' and compOrRein == 'Comp':
                clustObs = self.speciesOneSizeClusteringComp
            elif propertyType == 'Polarity' and compOrRein == 'Comp':
                clustObs = self.speciesOnePolarityClusteringComp
            elif propertyType == 'Charge' and compOrRein == 'Rein':
                clustObs = self.speciesOneChargeClusteringRein
            elif propertyType == 'Size' and compOrRein == 'Rein':
                clustObs = self.speciesOneSizeClusteringRein
            elif propertyType == 'Polarity' and compOrRein == 'Rein':
                clustObs = self.speciesOnePolarityClusteringRein
        elif whichSpecies == 2:
            if propertyType == 'Charge' and compOrRein == 'Comp':
                clustObs = self.speciesTwoChargeClusteringComp
            elif propertyType == 'Size' and compOrRein == 'Comp':
                clustObs = self.speciesTwoSizeClusteringComp
            elif propertyType == 'Polarity' and compOrRein == 'Comp':
                clustObs = self.speciesTwoPolarityClusteringComp
            elif propertyType == 'Charge' and compOrRein == 'Rein':
                clustObs = self.speciesTwoChargeClusteringRein
            elif propertyType == 'Size' and compOrRein == 'Rein':
                clustObs = self.speciesTwoSizeClusteringRein
            elif propertyType == 'Polarity' and compOrRein == 'Rein':
                clustObs = self.speciesTwoPolarityClusteringRein
        elif whichSpecies == 3:
            if propertyType == 'Charge' and compOrRein == 'Comp':
                clustObs = self.btwnChargeClusteringComp
            elif propertyType == 'Size' and compOrRein == 'Comp':
                clustObs = self.btwnSizeClusteringComp
            elif propertyType == 'Polarity' and compOrRein == 'Comp':
                clustObs = self.btwnPolarityClusteringComp
            elif propertyType == 'Charge' and compOrRein == 'Rein':
                clustObs = self.btwnChargeClusteringRein
            elif propertyType == 'Size' and compOrRein == 'Rein':
                clustObs = self.btwnSizeClusteringRein
            elif propertyType == 'Polarity' and compOrRein == 'Rein':
                clustObs = self.btwnPolarityClusteringRein
        return clustObs
    
    def calculateSignificance(self,clustType = 'DNDN', normType = None, propertyType = None, compOrRein = None, significanceLength = 30):
        # Norm type is only used for property clustering
        # supposed to be either DNDN1, DNDN2, or DNDNbtwn, depending on what kind of property
        if clustType == 'DNDN':
            clustObs = self.nonPolarizedDNDNClustering
            clustNorm = self.nonPolarizedDNDNClusteringNorm
        #elif clustType == 'DNDS':
        #    clustObs = self.nonPolarizedDNDSClustering
        #    clustNorm = self.nonPolarizedDNDSClusteringNorm
        elif clustType == 'DSDS':
            clustObs = self.nonPolarizedDSDSClustering
            clustNorm = self.nonPolarizedDSDSClusteringNorm
        elif clustType == "DNDNvsDSDS":
            clustObs = self.nonPolarizedDNDNClustering
            clustNorm = self.nonPolarizedDSDSClustering
        elif clustType == 'DNDN1vsbtwn':
            clustObs = self.speciesOneDNDNClustering
            clustNorm = self.btwnDNDNClustering
        elif clustType == 'DNDN2vsbtwn':
            clustObs = self.speciesTwoDNDNClustering
            clustNorm = self.btwnDNDNClustering
        elif clustType == 'Property':
            if normType == 'DNDN1':
                clustNorm = self.speciesOneDNDNClustering
                clustObs = self.helperGetPropertySpecified(1,propertyType,compOrRein)
            elif normType == 'DNDN2':
                clustNorm = self.speciesTwoDNDNClustering
                clustObs = self.helperGetPropertySpecified(2,propertyType,compOrRein)
            elif normType == 'DNDNbtwn':
                clustNorm = self.btwnDNDNClustering
                clustObs = self.helperGetPropertySpecified(3,propertyType,compOrRein)
            elif normType == "DNDN1vsbtwn":
                clustObs = self.helperGetPropertySpecified(1,propertyType,compOrRein)
                clustNorm = self.helperGetPropertySpecified(3,propertyType,compOrRein)
            elif normType == "DNDN2vsbtwn":
                clustObs = self.helperGetPropertySpecified(2,propertyType,compOrRein)
                clustNorm = self.helperGetPropertySpecified(3,propertyType,compOrRein)
        
        observedPairs = numpy.sum(clustObs[1:significanceLength+1])
        expectedPairs = numpy.sum(clustNorm[1:significanceLength+1])
        normalizedExpectation = clustNorm[1:significanceLength+1]/expectedPairs*observedPairs
        
        return stats.chisquare(clustObs[1:significanceLength+1],normalizedExpectation)[1]
    
    def calculateSignificanceProperty(self,whichSpecies,propertyType = "Charge",compOrRein = "Comp",significanceLength = 20):
        withinYes = sum(self.helperGetPropertySpecified(whichSpecies,propertyType,compOrRein)[1:significanceLength+1])
        if whichSpecies == 1:
            withinNo = sum(self.speciesOneDNDNClustering[1:significanceLength+1]) - withinYes
        elif whichSpecies == 2:
            withinNo = sum(self.speciesTwoDNDNClustering[1:significanceLength+1]) - withinYes
        betweenYes = sum(self.helperGetPropertySpecified(3,propertyType,compOrRein)[1:significanceLength+1])
        betweenNo = sum(self.btwnDNDNClustering[1:significanceLength+1]) - betweenYes
        conTable = numpy.array([[withinYes,withinNo],[betweenYes,betweenNo]])
        print(conTable)
        return stats.chi2_contingency(conTable,lambda_="log-likelihood")
    
    def clusteringExcess(self,NPorP,length=30):
        if NPorP == "Nonpolarized":
            DNDN = self.nonPolarizedDNDNClustering/self.nonPolarizedDNDNClusteringNorm
            #DNDS = self.nonPolarizedDNDSClustering/self.nonPolarizedDNDSClusteringNorm
            DSDS = self.nonPolarizedDSDSClustering/self.nonPolarizedDSDSClusteringNorm
            #DNDN = DNDN/numpy.mean(DNDN[80:121])
            #DNDS = DNDS/numpy.mean(DNDS[80:121])
            #DSDS = DSDS/numpy.mean(DSDS[80:121])
            #DNDN_excess = sum(DNDN[1:length+1]-1)
            #DNDS_excess = sum(DNDS[1:length+1]-1)
            #DSDS_excess = sum(DSDS[1:length-1]-1)
            DNDN_excess = sum(DNDN[1:length+1]-numpy.mean(DNDN[80:121]))
            #DNDS_excess = sum(DNDS[1:length+1]-numpy.mean(DNDS[80:121]))
            DSDS_excess = sum(DSDS[1:length-1]-numpy.mean(DSDS[80:121]))
            #return {"DNDN":DNDN_excess,"DNDS":DNDS_excess,"DSDS":DSDS_excess}
            return {"DNDN":DNDN_excess,"DSDS":DSDS_excess}
        elif NPorP == "Polarized":
            DNDN1 = self.speciesOneDNDNClustering/self.speciesOneDNDNClusteringNorm
            DNDN2 = self.speciesTwoDNDNClustering/self.speciesTwoDNDNClusteringNorm
            DNDNcross = self.btwnDNDNClustering/self.btwnDNDNClusteringNorm
            #DNDN1 = DNDN1/numpy.mean(DNDN1[80:121])
            #DNDN2 = DNDN2/numpy.mean(DNDN2[80:121])
            #DNDNcross = DNDNcross/numpy.mean(DNDNcross[80:121])
            #DNDN1_excess = sum(DNDN1[1:length+1]-1)
            #DNDN2_excess = sum(DNDN2[1:length+1]-1)
            #DNDNcross_excess = sum(DNDNcross[1:length+1]-1)
            DNDN1_excess = sum(DNDN1[1:length+1]-numpy.mean(DNDN1[80:121]))
            DNDN2_excess = sum(DNDN2[1:length+1]-numpy.mean(DNDN2[80:121]))
            DNDNcross_excess = sum(DNDNcross[1:length+1]-numpy.mean(DNDNcross[80:121]))
            return {"DNDN1":DNDN1_excess,"DNDN2":DNDN2_excess,"DNDNcross":DNDNcross_excess}
        else:
            print("Unknown option selected")
            return
    
    def propertyClusteringExcess(self,propertyType,compOrRein,length=10):
        if propertyType == "Charge":
            comp1 = self.speciesOneChargeClusteringComp
            comp2 = self.speciesTwoChargeClusteringComp
            compCross = self.btwnChargeClusteringComp
            rein1 = self.speciesOneChargeClusteringRein
            rein2 = self.speciesTwoChargeClusteringRein
            reinCross = self.btwnChargeClusteringRein
        elif propertyType == "Size":
            comp1 = self.speciesOneSizeClusteringComp
            comp2 = self.speciesTwoSizeClusteringComp
            compCross = self.btwnSizeClusteringComp
            rein1 = self.speciesOneSizeClusteringRein
            rein2 = self.speciesTwoSizeClusteringRein
            reinCross = self.btwnSizeClusteringRein
        elif propertyType == "Polarity":
            comp1 = self.speciesOnePolarityClusteringComp
            comp2 = self.speciesTwoPolarityClusteringComp
            compCross = self.btwnPolarityClusteringComp
            rein1 = self.speciesOnePolarityClusteringRein
            rein2 = self.speciesTwoPolarityClusteringRein
            reinCross = self.btwnPolarityClusteringRein
        else:
            print("Unknown type selected")
            return
        if compOrRein == "Comp":
            comp1 = comp1/self.speciesOneDNDNClustering
            comp2 = comp2/self.speciesTwoDNDNClustering
            compCross = compCross/self.btwnDNDNClustering
            comp1 = comp1/numpy.mean(comp1[80:121])
            comp2 = comp2/numpy.mean(comp2[80:121])
            compCross = compCross/numpy.mean(compCross[80:121])
            comp1_excess = sum(comp1[1:length+1] - numpy.mean(comp1[80:121]))
            comp2_excess = sum(comp2[1:length+1] - numpy.mean(comp2[80:121]))
            compCross_excess = sum(compCross[1:length+1] - numpy.mean(compCross[80:121]))
            return {"comp1":comp1_excess,"comp2":comp2_excess,"compCross":compCross_excess}
        elif compOrRein == "Rein":
            rein1 = rein1/self.speciesOneDNDNClustering
            rein2 = rein2/self.speciesTwoDNDNClustering
            reinCross = reinCross/self.btwnDNDNClustering
            rein1 = rein1/numpy.mean(rein1[80:121])
            rein2 = rein2/numpy.mean(rein2[80:121])
            reinCross = reinCross/numpy.mean(reinCross[80:121])
            rein1_excess = sum(rein1[1:length+1] - numpy.mean(rein1[80:121]))
            rein2_excess = sum(rein2[1:length+1] - numpy.mean(rein2[80:121]))
            reinCross_excess = sum(reinCross[1:length+1] - numpy.mean(reinCross[80:121]))
            return {"rein1":rein1_excess,"rein2":rein2_excess,"reinCross":reinCross_excess}
        else:
            print("Unknown type selected")
            return
    

