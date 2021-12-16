#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on Wed Mar 24 13:55:58 2021

@author: Ozan, Morgane T.

orsum enrichment analysis results filtering tool
This code gets the list of ranked gene set/pathway codes resulting from an
enrichment analysis and summarizes them at different levels applying more
relaxed rules at each step.
"""

from termCombinationLib import *
from plotFunctions import orsum_plot
from argparse import ArgumentParser, SUPPRESS
import os


def argumentParserFunction():
	"""
	"""
	##Parse arguments
	parser = ArgumentParser(description = 'This tool summarizes enrichment results', add_help = False)
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	# Add back help
	optional.add_argument('-h', '--help', action = 'help', default = SUPPRESS, help = 'show this help message and exit')
	# required arguments
	required.add_argument('--gmt', required = True, help = 'Path of the GMT file.')
	required.add_argument('--files', required = True, nargs = '+', help = 'Paths of the enrichment result files.')
	# optional arguments
	optional.add_argument('--outputFolder', default = ".", help = 'Path of the output result files. If it is not specified, results are written to the current directory.')
	optional.add_argument('--maxRepSize', type = int, default = 2000, help = 'The maximum size of a representative term. Terms bigger than this will not be discarded but also will not be used to represent other terms. By default, maxRepSize = 2000')
	optional.add_argument('--minTermSize', type = int, default = 10, help = 'The minimum size of the terms to be processed. Smaller terms will be discarded. By default, minTermSize = 10')
	optional.add_argument('--rules', type = int, nargs = '*', help = 'List of ordered numbers of the rules to apply while summarizing. By default, all the rules from 1 to 5 are run.')
	optional.add_argument('--numberOfTermsToPlot', type = int, default = 50, help = 'The number of representative terms to be presented in barplot and heatmap. By default, numberOfTermsToPlot = 50')
	optional.add_argument('--outputAll', action = 'store_true', help = 'When this option is used, a summary file is created after applying each rule, otherwise only final summary is created')
	return(parser)

if __name__ == "__main__":
	# Command-line interface
	parser = argumentParserFunction()
	args = parser.parse_args()
	argsDict = vars(args)

	# Parameters
	gmtPath=argsDict['gmt']
	tbsFiles=argsDict['files']
	outputAll=argsDict['outputAll']
	maxRepresentativeTermSize=argsDict['maxRepSize']
	minTermSize=argsDict['minTermSize']
	rulesToUseNo=argsDict['rules']
	numberOfTermsToPlot=argsDict['numberOfTermsToPlot']
	outputFolder=argsDict['outputFolder']



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
		tbsGsIDs=readTBSListFile(tf)

		tbsGsIDsRUT=removeUnknownTermsFromTBS(tbsGsIDs, geneSetsDict)
		difRUT=len(tbsGsIDs)-len(tbsGsIDsRUT)
		if(difRUT>0):
			print(difRUT, 'terms are not in GMT, they are removed.')

		tbsGsIDsRTS=removeTermsSmallerThanMinTermSize(tbsGsIDsRUT, geneSetsDict, minTermSize)
		difRTS=len(tbsGsIDsRUT)-len(tbsGsIDsRTS)
		if(difRTS>0):
			print(difRTS, 'terms are smaller than minTermSize, they are removed.')

		if len(tbsGsIDs)==0:
			print('The IDs in an input file do not match the IDs in the GMT file. Please check your command, the GMT file and input files, and be sure that they are consistent.')
			exit()
		tbsGsIDsList.append(tbsGsIDs)


	termSummary=initializeTermSummary(tbsGsIDsList)


	#Rules are listed here as tuples (function name, explanation, rule no).
	allRules=[
		(supertermRepresentsLessSignificantSubterm, 'Super <- Sub (w worse rank) || Superterms represent their less significant (worse ranked) subterms. This includes equal terms, i.e. the terms that annotate exactly the same set of genes.', 1),
		(subtermRepresentsLessSignificantSimilarSuperterm, 'Sub <- Super (w worse rank, subterm is 75% of superterm) || Subterms represent their less significant  superterms if subterm gene set constitutes 75% of the superterm gene set.', 2),
		(subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower, 'Sub <- Super (w worse rank less representative power) || Subterms represent their less significant and with less representative power superterms. Representative power is the number of terms they represent.', 3),
		(commonSupertermInListRepresentsSubtermsWithLessRepresentativePower, 'Common Super <- (Sub Sub) (w less or equal representative power) || Pairs of subterms are represented by their common superterms that have more representative power.', 4),
		(supertermRepresentsSubtermLargerThanMaxRep, 'Large Super <- Large Sub || Superterms larger than maxRepSize represent their less significant subterms which are also larger than maxRepSize', 5)
	]


	#This rule is run by default if there are multiple enrichment results
	multipleListsUnifyRule=(recurringTermsUnified,'Same terms in multiple lists are unified')

	#rulesToApply contains the rules to be applied, based on arguments (by default
	#all rules are applied)
	if rulesToUseNo==None:
		rulesToApply=allRules
	else:
		try:
			rulesToApply=[allRules[int(r)-1] for r in rulesToUseNo]
		except IndexError:
			print('For the rules parameter, numbers between 1-5 (inclusive) must be given')
			exit()
		except ValueError:
			print('For the rules parameter, numbers between 1-5 (inclusive) must be given')
			exit()



	#Rules are applied one by one.
	#After application of each rule, previously applied rules are applied again
	#starting from the first.


	if outputFolder[-1]!=os.sep:
		outputFolder=outputFolder+os.sep

	if not os.path.isdir(outputFolder):
		os.mkdir(outputFolder)

	print('\n\n')
	if(len(tbsFiles)==1):
		print('Initial term number:',len(termSummary),'\n')
	else:
		print('Initial term number (recurring terms in different lists are not merged yet, each one is counted):',len(termSummary),'\n')

		#Apply multipleListsUnifyRule to unify the same terms from multiple lists
		print(multipleListsUnifyRule[1])
		termSummary=applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, multipleListsUnifyRule[0])
		print('Representing term number:',len(termSummary),'\n')


	processStep=1
	for i in range(len(rulesToApply)):
		#Apply rule
		print(rulesToApply[i][1])
		termSummary=applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, rulesToApply[i][0])

		#Apply previously applied rules again
		for j in range(0,i):
			termSummary=applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, rulesToApply[j][0])

		print('Representing term number:',len(termSummary),'\n')

		fileName=outputFolder+'termSummary'+'{:02d}'.format(i+1)+'-Rule'+'{:02d}'.format(rulesToApply[i][2])

		if(outputAll or i==len(rulesToApply)-1):
			#writeTermSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Detailed.tsv', outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Summary.tsv')
			writeTermSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, fileName+'-Detailed.tsv', fileName+'-Summary.tsv')
			writeHTMLSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, fileName+'.html')

		if(i==len(rulesToApply)-1):
			#orsum_plot(outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Summary.tsv', outputFolder, 50)
			orsum_plot(fileName+'-Summary.tsv', outputFolder, numberOfTermsToPlot)

		processStep=processStep+1

