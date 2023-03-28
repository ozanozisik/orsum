#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
@author: Ozan, Morgane T.

orsum enrichment analysis results filtering tool.
This code gets the list of gene set/pathway IDs resulting from an
enrichment analysis and lets the significant subterms represent their
less significant subterms.
"""

from termCombinationLib import *
from plotFunctions import orsum_plot
from argparse import ArgumentParser, SUPPRESS
import os


VERSION='1.5.0'


def argumentParserFunction():
	"""
	"""
	##Parse arguments
	parser = ArgumentParser(description = 'orsum summarizes enrichment results', add_help = False)
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	# Add back help
	optional.add_argument('-h', '--help', action = 'help', default = SUPPRESS, help = 'Show this help message and exit')
	# Add version
	optional.add_argument('-v', '--version', action = 'version', version=VERSION)
	# required arguments
	required.add_argument('--gmt', required = True, help = 'Path of the GMT file.')
	required.add_argument('--files', required = True, nargs = '+', help = 'Paths of the enrichment result files.')
	# optional arguments
	optional.add_argument('--fileAliases', nargs = '+', default=None, help = 'Aliases for input enrichment result files to be used in orsum results')
	optional.add_argument('--outputFolder', default = ".", help = 'Path for the output result files. If it is not specified, results are written to the current directory.')
	optional.add_argument('--maxRepSize', type = int, default = int(1E6), help = 'The maximum size of a representative term. Terms larger than this will not be discarded but also will not be used to represent other terms. By default, it is larger than any annotation term (1E6), which means that it has no effect.')
	optional.add_argument('--maxTermSize', type = int, default = int(1E6), help = 'The maximum size of the terms to be processed. Larger terms will be discarded. By default, it is larger than any annotation term (1E6), which means that it has no effect.')
	optional.add_argument('--minTermSize', type = int, default = 10, help = 'The minimum size of the terms to be processed. Smaller terms will be discarded. By default, minTermSize = 10')
	optional.add_argument('--numberOfTermsToPlot', type = int, default = 50, help = 'The number of representative terms to be presented in barplot and heatmap. By default (and maximum), numberOfTermsToPlot = 50')
	return(parser)

if __name__ == "__main__":
	# Command-line interface
	parser = argumentParserFunction()
	args = parser.parse_args()
	argsDict = vars(args)

	# Parameters
	gmtPath=argsDict['gmt']
	tbsFiles=argsDict['files']#to-be-summarized files
	fileAliases=argsDict['fileAliases']
	outputFolder=argsDict['outputFolder']
	maxRepresentativeTermSize=argsDict['maxRepSize']
	maxTermSize=argsDict['maxTermSize']
	minTermSize=argsDict['minTermSize']
	numberOfTermsToPlot=argsDict['numberOfTermsToPlot']


	if fileAliases is None:
		fileAliases=[]
		for fname in tbsFiles:
			fileAliases.append(os.path.basename(fname))
	else:
		if len(fileAliases)!=len(tbsFiles):
			print('Number of input files and aliases do not match.\n')
			exit()


	if outputFolder[-1]!=os.sep:
		outputFolder=outputFolder+os.sep

	if not os.path.isdir(outputFolder):
		os.mkdir(outputFolder)


	logFile=open(outputFolder+'log.txt', 'w')
	logFile.write('orsum '+VERSION+'\n\n')
	for k,v in argsDict.items():
		logFile.write("{}:\t{}\n".format(k,v))
	logFile.write('\n')

	#Read the GMT file.
	#The GMT file contains one gene set/pathway at each line, with no header.
	#A line starts with gene set code, then gene set name, then gene IDs, each
	#separated by tab. As long as the gene IDs are consistent throughout the GMT
	#file it doesn't matter which ID is used, only the overlaps between
	#different gene sets are checked.

	geneSetsDict, gsIDToGsNameDict=readGmtFile(gmtPath)



	tbsGsIDsList=[]
	for tf in tbsFiles:
		print()
		print('Processing', tf)
		logFile.write('\nProcessing {}\n'.format(tf))
		tbsGsIDs=readTBSListFile(tf)

		tbsGsIDsRUT=removeUnknownTermsFromTBS(tbsGsIDs, geneSetsDict)
		difRUT=len(tbsGsIDs)-len(tbsGsIDsRUT)
		if(difRUT>1):
			print('{} terms are not in GMT, they are removed.'.format(difRUT))
			logFile.write('{} terms are not in GMT, they are removed.\n'.format(difRUT))
		elif(difRUT==1):
			print('{} term is not in GMT, it is removed.'.format(difRUT))
			logFile.write('{} term is not in GMT, it is removed.\n'.format(difRUT))

		tbsGsIDsRTS=removeTermsSmallerThanMinTermSize(tbsGsIDsRUT, geneSetsDict, minTermSize)
		difRTS=len(tbsGsIDsRUT)-len(tbsGsIDsRTS)
		if(difRTS>1):
			print('{} terms are smaller than minTermSize={}, they are removed.'.format(difRTS, minTermSize))
			logFile.write('{} terms are smaller than minTermSize={}, they are removed.\n'.format(difRTS, minTermSize))
		elif(difRTS==1):
			print('{} term is smaller than minTermSize={}, it is removed.'.format(difRTS, minTermSize))
			logFile.write('{} term is smaller than minTermSize={}, it is removed.\n'.format(difRTS, minTermSize))
		
		tbsGsIDsRTL=removeTermsLargerThanMaxTermSize(tbsGsIDsRTS, geneSetsDict, maxTermSize)
		difRTL=len(tbsGsIDsRTS)-len(tbsGsIDsRTL)
		if(difRTL>1):
			print('{} terms are larger than maxTermSize={}, they are removed.'.format(difRTL, maxTermSize))
			logFile.write('{} terms are larger than maxTermSize={}, they are removed.\n'.format(difRTL, maxTermSize))
		elif(difRTL==1):
			print('{} term is larger than maxTermSize={}, it is removed.'.format(difRTL, maxTermSize))
			logFile.write('{} term is larger than maxTermSize={}, it is removed.\n'.format(difRTL, maxTermSize))
		
		tbsGsIDsFinal=tbsGsIDsRTL.copy()

		if len(tbsGsIDsFinal)==0:
			print('There is no term left to be summarized from this input file. A possible reason is that IDs in the input file do not match the IDs in the GMT file. Another possible reason is setting minTermSize parameter too high. Please check your command, the GMT file and input files.')
			logFile.write('There is no term left to be summarized from this input file. A possible reason is that IDs in the input file do not match the IDs in the GMT file. Another possible reason is setting minTermSize parameter too high. Please check your command, the GMT file and input files.\n')
		else:
			tbsGsIDsList.append(tbsGsIDsFinal)

	#termSummary is a list, each element is a list that contains
	#term ID, the list of represented terms, rank
	termSummary=initializeTermSummary(tbsGsIDsList)


	logFile.write('\n')


	#Rules: (function name, description).
	#This rule is run by default if there are multiple enrichment results
	multipleListsUnifyRule=(recurringTermsUnified,'Same terms in multiple lists are unified')

	supertermRepresentsLessSignificantSubtermRule=(supertermRepresentsLessSignificantSubterm, 'Superterms represent their less significant (worse ranked) subterms. This includes equal terms, i.e. the terms that annotate exactly the same set of genes.')


	print('\n\n')
	
	if(len(tbsGsIDsList)==0):
		print('There is no valid file to be summarized.')
		logFile.write('There is no valid file to be summarized.\n')
		exit()
	elif(len(tbsGsIDsList)==1):
		print('Initial term number: {}\n'.format(len(termSummary)))
		logFile.write('Initial term number: {}\n\n'.format(len(termSummary)))
	else:
		print('Initial term number (recurring terms in different lists are not merged yet, each one is counted): {}\n'.format(len(termSummary)))
		logFile.write('Initial term number (recurring terms in different lists are not merged yet, each one is counted): {}\n\n'.format(len(termSummary)))
		#Apply multipleListsUnifyRule to unify the same terms from multiple lists
		print(multipleListsUnifyRule[1])
		logFile.write(multipleListsUnifyRule[1]+'\n')
		termSummary=applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, multipleListsUnifyRule[0])
		print('Representing term number: {}\n'.format(len(termSummary)))
		logFile.write('Representing term number: {}\n\n'.format(len(termSummary)))

	#Apply rule
	print(supertermRepresentsLessSignificantSubtermRule[1])
	logFile.write(supertermRepresentsLessSignificantSubtermRule[1]+'\n')
	termSummary=applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, supertermRepresentsLessSignificantSubtermRule[0])
	print('Representing term number:{}\n'.format(len(termSummary)))
	logFile.write('Representing term number:{}\n\n'.format(len(termSummary)))

	fileName=outputFolder+'filteredResult'
	

	writeTermSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, fileAliases, fileName+'-Detailed.tsv', fileName+'-Summary.tsv')
	writeHTMLSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, fileAliases, fileName+'.html')
	writeRepresentativeToRepresentedIDsFile(termSummary, fileName+'IDMapping.tsv')
	orsum_plot(fileName+'-Summary.tsv', outputFolder, numberOfTermsToPlot)

	writeTermSummaryFileClustered(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, fileAliases, fileName+'-SummaryClustered.tsv', numberOfTermsToPlot)
	orsum_plot(fileName+'-SummaryClustered.tsv', outputFolder, numberOfTermsToPlot, heatmapName = 'HeatmapClustered')
	os.remove(fileName+'-SummaryClustered.tsv')

	logFile.close()
