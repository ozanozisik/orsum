#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Ozan

Functions for filtering
"""

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import pandas as pd

##############################################################################



def applyRule(termSummary, termIdToGenesDict, maxRepresentativeTermSize, process):
	"""
	This function applies the specified rule ("process") on each pair of terms.

	:param list termSummary: Representative term list to be summarized with the application of rules. It is a list, each element is a list that contains term ID, the list of represented terms, rank
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
	:param function process: The rule to be applied.
	:return: **termSummary** (*list*) – Representative term list after applying the rule
	"""

	#Starting with the top terms, term pairs are picked from the termSummary
	#and the rule is applied to them.
	for idNo in range(len(termSummary)-1):
		termId=termSummary[idNo][0]
		if termId!=-1:#Check if the term is still a representative term
			for idNo2 in range(idNo+1, len(termSummary)):
				termId2=termSummary[idNo2][0]
				if  termId2!=-1:#Check if the term is still a representative term
					termSummary=process(termSummary, termIdToGenesDict, maxRepresentativeTermSize, idNo, idNo2, termId, termId2)

	#Remove terms that are represented by other terms
	termSummary=[e for e in termSummary if e[0]!=-1]
	#Sort termSummary by rank (first term has the best/smallest rank)
	termSummary.sort(key=lambda x: x[2])

	return termSummary


##############################################################################
##############################################################################
##############################################################################

'''
General description on how rules work:

Rules take information on two terms from the termSummary.
These two terms are evaluated whether they satisfy the conditions of the rule.
If term A will represent term B, the terms represented by term B
(this includes itself) are appended to the list of terms represented by
term A. Term B's termId which is stored in termSummary[idNoB][0] is set to -1
to mark that it is not a representative term any more. termId information is
not lost because it is already copied under term A (termSummary[idNoA][1]).
'''


def recurringTermsUnified(termSummary, termIdToGenesDict, maxRepresentativeTermSize, idNo, idNo2, termId, termId2):
	'''
	Recurring terms coming from multiple lists are unified.
	This rule is run by default if there are multiple enrichment results.

	:param list termSummary: Representative term list to be summarized
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
	:param int idNo: Index of the first term handled in this function in termSummary
	:param int idNo: Index of the second term handled in this function in termSummary
	:param str termId: Term ID of the first term handled in this function
	:param str termId2: Term ID of the second term handled in this function
	:return: **termSummary** (*list*) – Representative term list after applying the rule
	'''

	if(termId==termId2):
		#Terms represented by the second term are copied under the first term
		for termRepresentedByCoveredTerm in termSummary[idNo2][1]:
			if termRepresentedByCoveredTerm not in termSummary[idNo][1]:
				termSummary[idNo][1].append(termRepresentedByCoveredTerm)
		termSummary[idNo2][0]=-1
		#Minimum of the ranks is assigned as the rank for this representative term
		termSummary[idNo][2]=min(termSummary[idNo][2], termSummary[idNo2][2])
	return termSummary



def supertermRepresentsLessSignificantSubterm(termSummary, termIdToGenesDict, maxRepresentativeTermSize, idNo, idNo2, termId, termId2):
	'''
	Superterms represent their subterms that are less significant.
	The rule also works for the terms with the same list of genes.

	:param list termSummary: Representative term list to be summarized
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
	:param int idNo: Index of the first term handled in this function in termSummary
	:param int idNo: Index of the second term handled in this function in termSummary
	:param str termId: Term ID of the first term handled in this function
	:param str termId2: Term ID of the second term handled in this function
	:return: **termSummary** (*list*) – Representative term list after applying the rule
	'''

	geneSet1=termIdToGenesDict[termId]#Supposed to be superset and representative
	geneSet2=termIdToGenesDict[termId2]#Supposed to be subset and be represented
	if(len(geneSet1)<=maxRepresentativeTermSize):
		if geneSet1.issuperset(geneSet2):
			#Terms represented by the second term are copied under the first term
			for termRepresentedByCoveredTerm in termSummary[idNo2][1]:
				if termRepresentedByCoveredTerm not in termSummary[idNo][1]:
					termSummary[idNo][1].append(termRepresentedByCoveredTerm)
			termSummary[idNo2][0]=-1
	return termSummary



##############################################################################
##############################################################################
##############################################################################



def readGmtFile(gmtPath):
	"""
	Read GMT file. In GMT file each line consists of term ID, term name and genes,
	all tab separated, no header.

	:param str gmtPath: Path of the GMT file
	:return: **termIdToGenesDict** (*dict*) – Dictionary mapping term IDs to set of genes.
	:return: **termIdToTermNameDict** (*dict*) – Dictionary mapping term IDs to term names.
	"""

	termIdToGenesDict=dict()#term ID to set of genes mapping
	termIdToTermNameDict=dict()#term ID to term name mapping
	try:
		f=open(gmtPath, 'r')
		lines=f.readlines()
		for line in lines:
			tokens=line.strip().split('\t')
			termId=tokens[0]
			termName=tokens[1]
			genes=set(tokens[2:])
			termIdToGenesDict[termId]=genes
			termIdToTermNameDict[termId]=termName
		f.close()
	except IOError:
		print("I/O error while reading gmt file.")
	return termIdToGenesDict, termIdToTermNameDict


##############################################################################
##############################################################################
##############################################################################


def readInputEnrichmentResultFile(inputEnrichmentResultFile):
	"""
	Read the input enrichment result file.
	Each line consists of a term ID, sorted from the most significant to the least

	:param str inputEnrichmentResultFile: Path of the enrichment result files
	:return: **termIdsList** (*list*) – List of term IDs to be summarized
	"""
	termIdsList=[]
	try:
		f=open(inputEnrichmentResultFile, 'r')
		lines=f.readlines()
		for line in lines:
			termIdsList.append(line.strip().strip('\"'))
		f.close()
	except IOError:
		print("I/O error while reading a file to be summarized.")
		exit()

	return termIdsList


def removeUnknownTerms(termIdsList, termIdToGenesDict):
	"""
	Remove unknown terms

	:param list termIdsList: Term IDs list
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:return: **termIdsList** (*list*) – Term IDs list after removal of unknown terms
	"""

	termIdsList=termIdsList.copy()

	idsToRemove=set()
	for termId in termIdsList:
		if not(termId in termIdToGenesDict.keys()):
			idsToRemove.add(termId)
	for idToRemove in idsToRemove:
		print(idToRemove, 'is not in gmt file')
		termIdsList.remove(idToRemove)
	return termIdsList


def removeTermsSmallerThanMinTermSize(termIdsList, termIdToGenesDict, minTermSize):
	"""
	Remove terms smaller than minTermSize

	:param list termIdsList: Term IDs list
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:param int minTermSize: The minimum size of the terms to be processed. Smaller terms are discarded.
	:return: **termIdsList** (*list*) – Term IDs list after removal of small terms
	"""
	termIdsList=termIdsList.copy()

	idsToRemove=set()
	for termId in termIdsList:
		if len(termIdToGenesDict[termId])<minTermSize:
			idsToRemove.add(termId)
	for idToRemove in idsToRemove:
			termIdsList.remove(idToRemove)
	return termIdsList


def removeTermsLargerThanMaxTermSize(termIdsList, termIdToGenesDict, maxTermSize):
	"""
	Remove terms larger than maxTermSize

	:param list termIdsList: Term IDs list
	:param dict termIdToGenesDict: Dictionary mapping term IDs to set of genes.
	:param int maxTermSize: The maximum size of the terms to be processed. Larger terms are discarded.
	:return: **termIdsList** (*list*) – Term IDs list after removal of large terms
	"""
	termIdsList=termIdsList.copy()

	idsToRemove=set()
	for termId in termIdsList:
		if len(termIdToGenesDict[termId])>maxTermSize:
			idsToRemove.add(termId)
	for idToRemove in idsToRemove:
			termIdsList.remove(idToRemove)
	return termIdsList

##############################################################################
##############################################################################
##############################################################################


def initializeTermSummary(termIdsListList):
	"""
	Initialize the term summary.
	Term summary is a list where each element is also a list describing
	the representative term.
	Each element in term summary contains
	term ID, list of represented terms, rank.
	Initially all terms are representative terms and represent themselves.

	:param list termIdsListList: List of terms obtained from enrichment results files
	:return: **termSummary** (*list*) – Initial representative term list
	"""

	termSummary=[]
	for termIdsListNo in range(len(termIdsListList)):
		termIdsList=termIdsListList[termIdsListNo]
		for idNo in range(len(termIdsList)):
			termId=termIdsList[idNo]
			rank=idNo+1
			termSummary.append([termId, [termId], rank])
	termSummary.sort(key=lambda x: x[2])#Ascending sort based on rank/score
	return termSummary


##############################################################################
##############################################################################
##############################################################################



def writeTermSummaryFile(termSummary, termIdToGenesDict, termIdToTermNameDict, termIdsListList, fileAliases, termSummaryFile, termSummaryFile2):
	'''
	Writes the results.
	'''
	try:
		#Detailed
		f=open(termSummaryFile, 'w')

		f.write('Representing term id\tRepresenting term name\tRepresenting term rank')
		for termIdsListNo in range(len(termIdsListList)):
			f.write('\t'+fileAliases[termIdsListNo]+ ' term id')
			f.write('\t'+fileAliases[termIdsListNo]+ ' term name')
			f.write('\t'+fileAliases[termIdsListNo]+ ' term rank')
		f.write('\n')

		for ts in termSummary:#For each representation
			f.write(ts[0]+'\t'+termIdToTermNameDict[ts[0]]+'\t'+str(ts[2])+'\n')

			mtr=np.empty([len(ts[1]), len(termIdsListList)*3],dtype=(np.unicode_, 10000))
			row=0
			for representedTerm in ts[1]:#For each represented term
				for termIdsListNo in range(len(termIdsListList)):#For each input enrichment result
					if representedTerm in termIdsListList[termIdsListNo]:
						mtr[row, termIdsListNo*3]=representedTerm
						mtr[row, termIdsListNo*3+1]=termIdToTermNameDict[representedTerm]
						mtr[row, termIdsListNo*3+2]=termIdsListList[termIdsListNo].index(representedTerm)+1
					else:
						mtr[row, termIdsListNo*3]=''
						mtr[row, termIdsListNo*3+1]=''
						mtr[row, termIdsListNo*3+2]=''
				row=row+1

			for r in range(row):
				f.write('\t\t')
				for c in range(mtr.shape[1]):
					f.write('\t'+mtr[r,c])
				f.write('\n')
		f.close()


		#Summary of summary
		f=open(termSummaryFile2, 'w')

		f.write('Representing term id\tRepresenting term name\tRepresenting term size\tRepresenting term rank\tRepresented term number')
		for termIdsListNo in range(len(termIdsListList)):
			f.write('\t'+fileAliases[termIdsListNo]+ ' term rank')
		f.write('\n')

		#For each representative term and for each input enrichment result, 
		#find the best rank from the terms represented by that representative term
		for ts in termSummary:
			f.write(ts[0]+'\t'+termIdToTermNameDict[ts[0]]+'\t'+str(len(termIdToGenesDict[ts[0]]))+'\t'+str(ts[2])+'\t'+str(len(ts[1])))

			for termIdsListNo in range(len(termIdsListList)):
				found=None
				for representedTerm in ts[1]:#for each represented term
					try:
						rank=termIdsListList[termIdsListNo].index(representedTerm)
						if found==None or rank<found:
							found=rank
					except ValueError:
						pass
				if found!=None:
					found=found+1 #rank starts from 1
				f.write('\t'+str(found))
			f.write('\n')
		f.close()


	except IOError:
		print("I/O error while writing term summary file.")


def writeHTMLSummaryFile(termSummary, termIdToGenesDict, termIdToTermNameDict, termIdsListList, fileAliases, termSummaryFile):
	'''
	Writes the results to an HTML file
	'''
	try:
		#Detailed
		f=open(termSummaryFile, 'w')
		f.write('<!DOCTYPE html>\n')
		f.write('<html>\n')
		f.write('<head>\n')
		f.write('<title>orsum result</title>\n')
		f.write('</head>\n')
		f.write('<body>\n')

		for ts in termSummary:
			#f.write(getTextForTSElement(ts, termIdToTermNameDict))
			f.write(getTextForTSElementMultiEnrichment(ts, termIdToGenesDict, termIdToTermNameDict, termIdsListList, fileAliases))

		f.write('</body>\n')
		f.write('</html>\n')

		f.close()



	except IOError:
		print("I/O error while writing term summary file.")




def writeTermSummaryFileClustered(termSummary, termIdToGenesDict, termIdToTermNameDict, termIdsListList, fileAliases, termSummaryFile, nbTerm):
	'''
	Writes the top results as clustered. The purpose of this file is to be
	consumed by the plot function to create a clustered heatmap.
	It can then be deleted (which is done by orsum.py).
	'''
	try:
				
		ranksPerInputFile=dict()
		for termIdsListNo in range(len(termIdsListList)):
			ranksPerInputFile[fileAliases[termIdsListNo]]=[]

		for ts in termSummary:
			for termIdsListNo in range(len(termIdsListList)):
				found=None
				for representedTerm in ts[1]:
					try:
						rank=termIdsListList[termIdsListNo].index(representedTerm)
						if found==None or rank<found:
							found=rank
					except ValueError:
						pass
				if found!=None:
					found=found+1 #rank starts from 1
				ranksPerInputFile[fileAliases[termIdsListNo]].append(found)
		
		dfRanksPerInputFile=pd.DataFrame.from_dict(ranksPerInputFile)
		
		
		#dfScaledRanksPerInputFile= (dfRanksPerInputFile - dfRanksPerInputFile.min()) / (dfRanksPerInputFile.max() - dfRanksPerInputFile.min())
		
		
		quartiles = dfRanksPerInputFile.quantile(q=[0.25, 0.50, 0.75])
				
		
		dfLinkageInput=dfRanksPerInputFile.copy(deep=True)
		
		for col in dfRanksPerInputFile.columns:
			dfLinkageInput.loc[(dfRanksPerInputFile[col]<=quartiles.loc[0.25,col]), col]=1
			dfLinkageInput.loc[((dfRanksPerInputFile[col]>quartiles.loc[0.25,col]) & (dfRanksPerInputFile[col]<=quartiles.loc[0.50,col])), col]=2
			dfLinkageInput.loc[((dfRanksPerInputFile[col]>quartiles.loc[0.50,col]) & (dfRanksPerInputFile[col]<=quartiles.loc[0.75,col])), col]=3
			dfLinkageInput.loc[(dfRanksPerInputFile[col]>quartiles.loc[0.75,col]), col]=4
			dfLinkageInput.loc[(dfRanksPerInputFile[col].isna()), col]=5
		
		
		dfLinkageInput=dfLinkageInput.iloc[0:nbTerm]
		
		
		
		Z = linkage(dfLinkageInput, 'average')
		dn = dendrogram(Z)
		
		
		#Summary of summary
		f=open(termSummaryFile, 'w')

		f.write('Representing term id\tRepresenting term name\tRepresenting term size\tRepresenting term rank\tRepresented term number')
		for termIdsListNo in range(len(termIdsListList)):
			f.write('\t'+fileAliases[termIdsListNo]+ ' term rank')
		f.write('\n')
		
		
		for leaf in dn['leaves']:
			ts=termSummary[leaf]
			f.write(ts[0]+'\t'+termIdToTermNameDict[ts[0]]+'\t'+str(len(termIdToGenesDict[ts[0]]))+'\t'+str(ts[2])+'\t'+str(len(ts[1])))

			for termIdsListNo in range(len(termIdsListList)):
				found=ranksPerInputFile[fileAliases[termIdsListNo]][leaf]
				if found==None:
					found='None'
				f.write('\t'+str(found))
			f.write('\n')
		
		for tsNo in range(nbTerm, len(termSummary)):
			ts=termSummary[tsNo]
			f.write(ts[0]+'\t'+termIdToTermNameDict[ts[0]]+'\t'+str(len(termIdToGenesDict[ts[0]]))+'\t'+str(ts[2])+'\t'+str(len(ts[1])))

			for termIdsListNo in range(len(termIdsListList)):
				found=ranksPerInputFile[fileAliases[termIdsListNo]][tsNo]
				if found==None:
					found='None'
				f.write('\t'+str(found))
			f.write('\n')
		
		f.close()

	except IOError:
		print("I/O error while writing term summary file.")




def writeRepresentativeToRepresentedIDsFile(termSummary, outputFile):
	'''
	Writes ID mapping between representative terms and represented terms
	'''
	f=open(outputFile, 'w')
	f.write('Representative\tRepresented')
	for ts in termSummary:
		representativeID=ts[0]
		for representedID in ts[1]:
			f.write('\n')
			f.write(representativeID+'\t'+representedID)
	f.close()



def getTextForTSElementMultiEnrichment(ts, termIdToGenesDict, termIdToTermNameDict, termIdsListList, fileAliases):
	txt='\n'
	txt=txt+'<details>'+'\n'
	#txt=txt+'<summary>'+ts[0]+' '+termIdToTermNameDict[ts[0]]+' '+str(ts[2])+'</summary>'+'\n'
	txt=txt+'<summary>'+ts[0]+' '+termIdToTermNameDict[ts[0]]+'</summary>'+'\n'

	for termIdsListNo in range(len(termIdsListList)):
		txt=txt+'\t'+'<p style="margin-left:40px">'+'\n'
		if(len(fileAliases)>1):
			txt=txt+'\t'+fileAliases[termIdsListNo]+'<br>'+'\n'
		for representedTerm in ts[1]:
			if representedTerm in termIdsListList[termIdsListNo]:
				txt=txt+'\t'+representedTerm+' '+termIdToTermNameDict[representedTerm]+' (rank: '+str(termIdsListList[termIdsListNo].index(representedTerm)+1)+', term size: '+ str(len(termIdToGenesDict[representedTerm])) +')<br>'+'\n'
		txt=txt+'\t'+'<br>'+'\n'
		txt=txt+'\t'+'</p>'+'\n'
	txt=txt+'</details>'+'\n'
	return txt


