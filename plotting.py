'''
Created on Dec 16, 2019

@author: PiTav

This is the script to make the final figures, only for python3

'''

from clusteringClasses import *
import pickle
import matplotlib.pyplot as plt
from matplotlib import gridspec

dataHuman = pickle.load(open('results_oop_code/vertebrate_clustering/hg38_ponAbe2_nomLeu3_mammalian.p','rb'),encoding = 'latin1')
dataYeast = pickle.load(open('results_oop_code/Scer_Smik_Skud_clustering.p','rb'),encoding = 'latin1')
dataDrosophila = pickle.load(open('results_oop_code/Dmel_Dsim_Dyak_Tai18E2_clustering.p','rb'),encoding = 'latin1')
dataArabidopsis = pickle.load(open('results_oop_code/Atha_Alyr_Crub_clustering.p','rb'),encoding = 'latin1')

plt.rcParams['legend.loc'] = 'upper right'

fig = plt.figure(figsize = (7,7), constrained_layout=True)
ax_big = fig.add_subplot()
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_big.set_xlabel('Distance in Codons',fontsize = 12)
ax_big.set_ylabel('Clustering',fontsize = 12)
gs = gridspec.GridSpec(2,2,figure = fig)
ax = fig.add_subplot(gs[0,0])
dataYeast.plotNonPolarized(plotTitle = 'Saccharomyces',showLegend = False, normalizeAsym = True, ax = ax, plotDNDS = False)
ax = fig.add_subplot(gs[0,1])
dataDrosophila.plotNonPolarized(plotTitle = 'Drosophila', normalizeAsym = True, ax = ax, plotDNDS = False)
ax = fig.add_subplot(gs[1,0])
dataArabidopsis.plotNonPolarized(plotTitle = 'Arabidopsis',showLegend = False, normalizeAsym = True, ax = ax, plotDNDS = False)
ax.set_xlabel(" ") # Cheat to make the axis tick numbers move out of the way
ax.set_ylabel(" ")
ax = fig.add_subplot(gs[1,1])
dataHuman['Clustering Holder'].plotNonPolarized(plotTitle = 'Primate',showLegend = False, normalizeAsym = True, ax = ax, plotDNDS = False)
fig.suptitle('Nonpolarized Clustering')
plt.savefig('results_oop_code/finalizedFigures/nonPolarizedClusteringFigure1.png')
plt.show()

dataYeast.calculateSignificance(significanceLength = 20)
dataDrosophila.calculateSignificance(significanceLength = 20)
dataArabidopsis.calculateSignificance(significanceLength = 20)
dataHuman['Clustering Holder'].calculateSignificance(significanceLength = 20)

dataYeast.calculateSignificance(clustType = "DNDNvsDSDS",significanceLength = 20)
dataDrosophila.calculateSignificance(clustType = "DNDNvsDSDS",significanceLength = 20)
dataArabidopsis.calculateSignificance(clustType = "DNDNvsDSDS",significanceLength = 20)
dataHuman['Clustering Holder'].calculateSignificance(clustType = "DNDNvsDSDS",significanceLength = 20)

## This is the code for intron plotting for primate data:
## To be run locally:
#dataHumanCI = pickle.load(open('results/hg38_ponAbe2_nomLeu3_intronSimulation_CI.p','rb'),encoding = 'latin1')
#
##ax.fill_between(range(501),data['DNDNlowerCI'],data['DNDNupperCI'],color = 'green',alpha = 0.5)
#dataHuman['Clustering Holder'].plotNonPolarized(plotIntron = True,plotDNDS = False)
#plt.fill_between(range(501),dataHumanCI['DSDSlowerCI'],dataHumanCI['DSDSupperCI'],color = 'blue',alpha = 0.25)
#plt.fill_between(range(501),dataHumanCI['DNDNlowerCI'],dataHumanCI['DNDNupperCI'],color = 'green',alpha = 0.25)
#ax = fig.add_subplot()


fig = plt.figure(figsize = (7,7), constrained_layout=True)
ax_big = fig.add_subplot()
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_big.set_xlabel('Distance in Codons',fontsize = 12)
ax_big.set_ylabel('Clustering',fontsize = 12)
gs = gridspec.GridSpec(2,2,figure = fig)
ax = fig.add_subplot(gs[0,0])
dataYeast.plotPolarized(plotTitle = 'Saccharomyces',
                        customNameOne = 'within S. cerevisiae**', 
                        customNameTwo = 'within S. mikatae*', 
                        customNameOut = 'between species', 
                        normalizeAsym = True, ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,1])
dataDrosophila.plotPolarized(plotTitle = 'Drosophila', 
                             customNameOne = 'within D. melanogaster**', 
                             customNameTwo = 'within D. simulans*', 
                             customNameOut = 'between species', 
                             normalizeAsym = True, ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,0])
dataArabidopsis.plotPolarized(plotTitle = 'Arabidopsis',
                              customNameOne = 'within A. thaliana**', 
                              customNameTwo = 'within A. lyrata**', 
                              customNameOut = 'between species',
                              normalizeAsym = True, ax = ax, minMaxX = (0,40))
ax.set_xlabel(" ") # Cheat to make the axis tick numbers move out of the way
ax.set_ylabel(" ")
ax = fig.add_subplot(gs[1,1])
dataHuman['Clustering Holder'].plotPolarized(plotTitle = 'Primate',
                                             customNameOne = 'within H. sapiens*', 
                                             customNameTwo = 'within P. abelii**', 
                                             customNameOut = 'between species',
                                             normalizeAsym = True, ax = ax, minMaxX = (0,40))
fig.suptitle('Polarized Clustering')
plt.savefig('results_oop_code/finalizedFigures/polarizedClusteringFigure2.png')
plt.show()

dataYeast.calculateSignificance(clustType = "DNDN1vsbtwn",significanceLength = 20)
dataDrosophila.calculateSignificance(clustType = "DNDN1vsbtwn",significanceLength = 20)
dataArabidopsis.calculateSignificance(clustType = "DNDN1vsbtwn",significanceLength = 20)
dataHuman['Clustering Holder'].calculateSignificance(clustType = "DNDN1vsbtwn",significanceLength = 20)

dataYeast.calculateSignificance(clustType = "DNDN2vsbtwn",significanceLength = 20)
dataDrosophila.calculateSignificance(clustType = "DNDN2vsbtwn",significanceLength = 20)
dataArabidopsis.calculateSignificance(clustType = "DNDN2vsbtwn",significanceLength = 20)
dataHuman['Clustering Holder'].calculateSignificance(clustType = "DNDN2vsbtwn",significanceLength = 20)


fig = plt.figure(figsize = (14,7), constrained_layout=True)
gs = gridspec.GridSpec(2,4,figure = fig)
ax_big = fig.add_subplot()
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_big.set_xlabel('Distance in Codons',fontsize = 12)
ax = fig.add_subplot(gs[0,0])
dataYeast.plotProperty('Charge','Comp',plotTitle = 'Saccharomyces',
                       customNameOne = 'within S. cerevisiae**', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Compensatory Clustering Fraction',fontsize = 12)
ax = fig.add_subplot(gs[0,1])
dataDrosophila.plotProperty('Charge','Comp',plotTitle = 'Drosophila', 
                            customNameOne = 'within D. melanogaster*', 
                            customNameTwo = 'within D. simulans**', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,2])
dataArabidopsis.plotProperty('Charge','Comp',plotTitle = 'Arabidopsis',
                             customNameOne = 'within A. thaliana**', 
                             customNameTwo = 'within A. lyrata**', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,3])
dataHuman['Clustering Holder'].plotProperty('Charge','Comp',plotTitle = 'Primate',
                                            customNameOne = 'within H. sapiens', 
                                            customNameTwo = 'within P. abelii*', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))

ax = fig.add_subplot(gs[1,0])
dataYeast.plotProperty('Charge','Rein', showLegend = True,
                       customNameOne = 'within S. cerevisiae**', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Reinforcing Clustering Fraction',fontsize = 12)
ax.set_xlabel(' ')
ax = fig.add_subplot(gs[1,1])
dataDrosophila.plotProperty('Charge','Rein', showLegend = True,
                            customNameOne = 'within D. melanogaster**', 
                            customNameTwo = 'within D. simulans*', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,2])
dataArabidopsis.plotProperty('Charge','Rein', showLegend = True,
                             customNameOne = 'within A. thaliana*', 
                             customNameTwo = 'within A. lyrata**', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,3])
dataHuman['Clustering Holder'].plotProperty('Charge','Rein', showLegend = True,
                                            customNameOne = 'within H. sapiens', 
                                            customNameTwo = 'within P. abelii', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))
fig.suptitle('Charge Clustering')
plt.savefig('results_oop_code/finalizedFigures/chargeClusteringFigure3.png')
plt.show()

dataYeast.calculateSignificanceProperty(1,"Charge","Comp",10)
dataDrosophila.calculateSignificanceProperty(1,"Charge","Comp",10)
dataArabidopsis.calculateSignificanceProperty(1,"Charge","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Charge","Comp",10)

dataYeast.calculateSignificanceProperty(2,"Charge","Comp",10)
dataDrosophila.calculateSignificanceProperty(2,"Charge","Comp",10)
dataArabidopsis.calculateSignificanceProperty(2,"Charge","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Charge","Comp",10)

dataYeast.calculateSignificanceProperty(1,"Charge","Rein",10)
dataDrosophila.calculateSignificanceProperty(1,"Charge","Rein",10)
dataArabidopsis.calculateSignificanceProperty(1,"Charge","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Charge","Rein",10)

dataYeast.calculateSignificanceProperty(2,"Charge","Rein",10)
dataDrosophila.calculateSignificanceProperty(2,"Charge","Rein",10)
dataArabidopsis.calculateSignificanceProperty(2,"Charge","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Charge","Rein",10)


fig = plt.figure(figsize = (14,7), constrained_layout=True)
ax_big = fig.add_subplot()
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_big.set_xlabel('Distance in Codons',fontsize = 12)
gs = gridspec.GridSpec(2,4,figure = fig)
ax = fig.add_subplot(gs[0,0])
dataYeast.plotProperty('Size','Comp',plotTitle = 'Saccharomyces',
                       customNameOne = 'within S. cerevisiae**', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Compensatory Clustering',fontsize = 12)
ax = fig.add_subplot(gs[0,1])
dataDrosophila.plotProperty('Size','Comp',plotTitle = 'Drosophila', 
                            customNameOne = 'within D. melanogaster', 
                            customNameTwo = 'within D. simulans', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,2])
dataArabidopsis.plotProperty('Size','Comp',plotTitle = 'Arabidopsis',
                             customNameOne = 'within A. thaliana*', 
                             customNameTwo = 'within A. lyrata*', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,3])
dataHuman['Clustering Holder'].plotProperty('Size','Comp',plotTitle = 'Primate',
                                            customNameOne = 'within H. sapiens', 
                                            customNameTwo = 'within P. abelii', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))

ax = fig.add_subplot(gs[1,0])
dataYeast.plotProperty('Size','Rein', showLegend = True,
                       customNameOne = 'within S. cerevisiae*', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Reinforcing Clustering',fontsize = 12)
ax.set_xlabel(' ')
ax = fig.add_subplot(gs[1,1])
dataDrosophila.plotProperty('Size','Rein', showLegend = True,
                            customNameOne = 'within D. melanogaster', 
                            customNameTwo = 'within D. simulans*', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,2])
dataArabidopsis.plotProperty('Size','Rein', showLegend = True,
                             customNameOne = 'within A. thaliana**', 
                             customNameTwo = 'within A. lyrata**', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,3])
dataHuman['Clustering Holder'].plotProperty('Size','Rein', showLegend = True,
                                            customNameOne = 'within H. sapiens', 
                                            customNameTwo = 'within P. abelii', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))
fig.suptitle('Size Clustering')
plt.savefig('results_oop_code/finalizedFigures/sizeClustering.png')
plt.show()


dataYeast.calculateSignificanceProperty(1,"Size","Comp",10)
dataDrosophila.calculateSignificanceProperty(1,"Size","Comp",10)
dataArabidopsis.calculateSignificanceProperty(1,"Size","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Size","Comp",10)

dataYeast.calculateSignificanceProperty(2,"Size","Comp",10)
dataDrosophila.calculateSignificanceProperty(2,"Size","Comp",10)
dataArabidopsis.calculateSignificanceProperty(2,"Size","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Size","Comp",10)

dataYeast.calculateSignificanceProperty(1,"Size","Rein",10)
dataDrosophila.calculateSignificanceProperty(1,"Size","Rein",10)
dataArabidopsis.calculateSignificanceProperty(1,"Size","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Size","Rein",10)

dataYeast.calculateSignificanceProperty(2,"Size","Rein",10)
dataDrosophila.calculateSignificanceProperty(2,"Size","Rein",10)
dataArabidopsis.calculateSignificanceProperty(2,"Size","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Size","Rein",10)


fig = plt.figure(figsize = (14,7), constrained_layout=True)
ax_big = fig.add_subplot()
ax_big.spines['top'].set_color('none')
ax_big.spines['bottom'].set_color('none')
ax_big.spines['left'].set_color('none')
ax_big.spines['right'].set_color('none')
ax_big.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
ax_big.set_xlabel('Distance in Codons',fontsize = 12)
gs = gridspec.GridSpec(2,4,figure = fig)
ax = fig.add_subplot(gs[0,0])
dataYeast.plotProperty('Polarity','Comp',plotTitle = 'Saccharomyces',
                       customNameOne = 'within S. cerevisiae*', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Compensatory Clustering',fontsize = 12)
ax = fig.add_subplot(gs[0,1])
dataDrosophila.plotProperty('Polarity','Comp',plotTitle = 'Drosophila', 
                            customNameOne = 'within D. melanogaster', 
                            customNameTwo = 'within D. simulans*', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,2])
dataArabidopsis.plotProperty('Polarity','Comp',plotTitle = 'Arabidopsis',
                             customNameOne = 'within A. thaliana*', 
                             customNameTwo = 'within A. lyrata*', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[0,3])
dataHuman['Clustering Holder'].plotProperty('Polarity','Comp',plotTitle = 'Primate',
                                            customNameOne = 'within H. sapiens*', 
                                            customNameTwo = 'within P. abelii', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))

ax = fig.add_subplot(gs[1,0])
dataYeast.plotProperty('Polarity','Rein', showLegend = True,
                       customNameOne = 'within S. cerevisiae*', 
                       customNameTwo = 'within S. mikatae*', 
                       customNameOut = 'between species', 
                       ax = ax, minMaxX = (0,40))
ax.set_ylabel('Reinforcing Clustering',fontsize = 12)
ax.set_xlabel(' ')
ax = fig.add_subplot(gs[1,1])
dataDrosophila.plotProperty('Polarity','Rein', showLegend = True,
                            customNameOne = 'within D. melanogaster*', 
                            customNameTwo = 'within D. simulans*', 
                            customNameOut = 'between species', 
                            ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,2])
dataArabidopsis.plotProperty('Polarity','Rein', showLegend = True,
                             customNameOne = 'within A. thaliana*', 
                             customNameTwo = 'within A. lyrata**', 
                             customNameOut = 'between species',
                             ax = ax, minMaxX = (0,40))
ax = fig.add_subplot(gs[1,3])
dataHuman['Clustering Holder'].plotProperty('Polarity','Rein', showLegend = True,
                                            customNameOne = 'within H. sapiens', 
                                            customNameTwo = 'within P. abelii', 
                                            customNameOut = 'between species',
                                            ax = ax, minMaxX = (0,40))
fig.suptitle('Polarity Clustering')
plt.savefig('results_oop_code/finalizedFigures/polarityClustering.png')
plt.show()

dataYeast.calculateSignificanceProperty(1,"Polarity","Comp",10)
dataDrosophila.calculateSignificanceProperty(1,"Polarity","Comp",10)
dataArabidopsis.calculateSignificanceProperty(1,"Polarity","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Polarity","Comp",10)

dataYeast.calculateSignificanceProperty(2,"Polarity","Comp",10)
dataDrosophila.calculateSignificanceProperty(2,"Polarity","Comp",10)
dataArabidopsis.calculateSignificanceProperty(2,"Polarity","Comp",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Polarity","Comp",10)

dataYeast.calculateSignificanceProperty(1,"Polarity","Rein",10)
dataDrosophila.calculateSignificanceProperty(1,"Polarity","Rein",10)
dataArabidopsis.calculateSignificanceProperty(1,"Polarity","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(1,"Polarity","Rein",10)

dataYeast.calculateSignificanceProperty(2,"Polarity","Rein",10)
dataDrosophila.calculateSignificanceProperty(2,"Polarity","Rein",10)
dataArabidopsis.calculateSignificanceProperty(2,"Polarity","Rein",10)
dataHuman['Clustering Holder'].calculateSignificanceProperty(2,"Polarity","Rein",10)





