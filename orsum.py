#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on Wed Mar 24 13:55:58 2021

@author: Ozan

Gene Set Summary
This code gets the list of ranked gene set/pathway codes resulting from an enrichment analysis and summarizes them at different levels applying more relaxed rules at each step.

"""

from termCombinationLib import *
from plotFunctions import *
from argparse import ArgumentParser, SUPPRESS
import os


parser = ArgumentParser(description='This tool summarizes enrichment results', add_help=False)
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# Add back help
optional.add_argument(
    '-h',
    '--help',
    action='help',
    default=SUPPRESS,
    help='show this help message and exit'
)
# required arguments
required.add_argument('--gmt', required=True, help='Path for the GMT file.')
required.add_argument('--hierarchyFile', required=True, help='Path for the hierarchy file. It is created if the file does not exist.')
required.add_argument('--files', required=True, nargs='+', help='Paths for the enrichment result files.')
# optional arguments
optional.add_argument('--outputFolder', default=".", help='Path for the output result files. If it is not specified, results are written to the current directory.')
optional.add_argument('--createHF', action='store_true', help='Creates the hierarchy file when this is used, otherwise tries to read, if absent creates it.')
optional.add_argument('--rules', type=int, nargs='*', help='List of ordered numbers of the rules to apply while summarizing. By default, all the rules from 1 to 9 are run.')
optional.add_argument('--maxRepSize', type=int, default=2000, help='The maximum size of a representative term. Terms bigger than this will not be discarded but also will not be used to represent other terms. By default, maxRepSize = 2000')
optional.add_argument('--outputAll', action='store_true', help='When this option is used, a summary file is created after applying each rule, otherwise only final summary is created')
args = parser.parse_args()
argsDict=vars(args)
#print(argsDict)



gmtPath=argsDict['gmt']
termHierarchyFile=argsDict['hierarchyFile']
createHF=argsDict['createHF']
tbsFiles=argsDict['files']
outputAll=argsDict['outputAll']

maxRepresentativeTermSize=int(argsDict['maxRepSize'])

#We read the GMT file and create term hiearchy for terms which have subset superset relationship. We don't use any prior information for hierarchy.
#The GMT file contains one gene set/pathway at each line, with no header. A line starts with gene set code, then gene set name, then gene IDs, each separated by tab. As long as the gene IDs are consistent throughout the GMT file it doesn't matter which ID is used, we only check overlaps between different gene sets.<br />
#Reactome and GO:BP GMT files are provided in the demo. You can obtain new versions from gProfiler data sources tab.<br />
#As creation of hierarchy file takes time for GO, you can create it once and use it.<br />
#Hierarchy files are provided in the demo.
geneSetsDict, gsIDToGsNameDict=readGmtFile(gmtPath)


#Based on the option provided by the user creates or reads existing hierarchy
#file. If IOError occurs the file is created.

if(createHF):#createHF is boolean
	#Create term hiearchy, save it to a file, and use it in the rest of the run
	hierarchyDict=createTermHierarchy(geneSetsDict, termHierarchyFile)
else:
	if os.path.exists(termHierarchyFile):
		#Read hierarchy from an existing hierarchy file
		hierarchyDict=readHierarchyFile(termHierarchyFile)
	else:
		print('Could not find hierarchy file, creating...')
		hierarchyDict=createTermHierarchy(geneSetsDict, termHierarchyFile)


#We read enrichment results to be summarized. These files must only contain ranked gene set IDs which are consistent with the GMT file.


if len(set(hierarchyDict.keys()).intersection(set(geneSetsDict.keys())))==0:
	if createHF:
		print('The IDs in the hierarchy file do not match the IDs in the GMT file. As --createHF parameter is used something must have gone wrong during hiearchy file creation.' )
		exit()
	else:
		print('The IDs in the hierarchy file do not match the IDs in the GMT file. Please check your command, the GMT file and the hierarchy file, and be sure that they are consistent.')
		exit()


tbsGsIDsList=[]
originalTermsSet=set()#Set of all the terms given in user supplied lists
for tf in tbsFiles:
	tbsGsIDs=readTBSListFile(tf)
	tbsGsIDs=removeUnknownGeneSetsFromTBS(tbsGsIDs, geneSetsDict)
	if len(tbsGsIDs)==0:
		print('The IDs in an input file do not match the IDs in the GMT file. Please check your command, the GMT file and input files, and be sure that they are consistent.')
		exit()
	tbsGsIDsList.append(tbsGsIDs)
	originalTermsSet.update(tbsGsIDs)

termSummary=initializeTermSummary(tbsGsIDsList)



#rulesToApply is a list of tuples that contains function names for the rules and their explanations. Rules are applied one by one according to this list. After application of each rule, previously applied rules are applied again starting from the first.

rulesToUseNo=argsDict['rules']

allRules=[
	(supertermRepresentsSubterm, 'Super <- Sub (w less Rank) || Superterms represent their less significant / lower ranked subterms', 1),
	(commonParentInListRepresentsTerms, 'Parent in list || Terms with a common parent are represented by the parent if the parent is in the original list', 2),
	(commonGrandparentInListRepresentsTerms, 'Grandparent in list || Terms with a common grandparent are represented by the grandparent if the grandparent is in the original list', 3),
	(commonParentGrandparentInListRepresentsTerms, 'Parent/Grandparent || Terms in which one\'s parent is other\'s grandparent are represented by this ancestor term if it is in the original list', 4),
	(subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower, 'Sub <- Super (w less rank less representative power) || Subterms represent their superterms with less significance / lower rank and less representative power. Representative power is the number of terms they represent.', 5),
	(subtermRepresentsSlightlyLowerRankedSuperterm, 'Sub <- Super (less Rank) || Subterms represent their slightly lower ranked superterms. Rank tolerance is set to 1.', 6),
	(commonParentRepresentsTerms, 'Parent || Terms with a common parent are represented by the parent (even if the parent is not in the original list)', 7),
	(commonGrandparentRepresentsTerms, 'Grandparent || Terms with a common grandparent are represented by the grandparent (even if the grandparent is not in the original list)', 8),
	(commonParentGrandparentRepresentsTerms, 'Parent/Grandparent || Terms in which one\'s parent is other\'s grandparent are represented by this ancestor term (even if it is not in the original list)', 9)
]


if rulesToUseNo==None:
	rulesToApply=allRules
else:
	try:
		rulesToApply=[allRules[int(r)-1] for r in rulesToUseNo]
	except IndexError:
		print('For the rules parameter, numbers between 1-9 (inclusive) must be given')
		exit()
	except ValueError:
		print('For the rules parameter, numbers between 1-9 (inclusive) must be given')
		exit()


multipleListsUnifyRule=(recurringTermsUnified,'Same terms in multiple lists are unified')


outputFolder=argsDict['outputFolder']
if outputFolder[-1]!=os.sep:
	outputFolder=outputFolder+os.sep

if not os.path.isdir(outputFolder):
	os.mkdir(outputFolder)

print('\n\n')
if(len(tbsFiles)==1):
	print('Initial term number:',len(termSummary),'\n')
else:
	print('Initial term number (recurring terms in different lists are not merged yet, each one is counted):',len(termSummary),'\n')

	#Applying default rule to unify same terms in multiple lists
	print(multipleListsUnifyRule[1])
	termSummary=applyRule(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, multipleListsUnifyRule[0])
	print('Representing term number:',len(termSummary),'\n')

processStep=1
#strAppliedRuleNos=''
for i in range(len(rulesToApply)):
	print(rulesToApply[i][1])
	#strAppliedRuleNos=strAppliedRuleNos+str(rulesToApply[i][2])
	termSummary=applyRule(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, rulesToApply[i][0])

	#Applying previously applied rules again
	for j in range(0,i):
		termSummary=applyRule(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, rulesToApply[j][0])

	print('Representing term number:',len(termSummary),'\n')

	if(outputAll or i==len(rulesToApply)-1):
		writeTermSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Detailed.tsv', outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Summary.tsv')

	if(i==len(rulesToApply)-1):
		orsum_plot(outputFolder+'termSummary'+str(rulesToApply[i][2])+'-Summary.tsv', outputFolder, 50)

	processStep=processStep+1
