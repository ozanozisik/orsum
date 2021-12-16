#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Nov 30 13:57:07 2020

@author: Ozan
"""

import numpy as np
import os

##############################################################################



def applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, process):
	"""
	This function applies the specified rule ("process") on each pair of terms.

	:param list termSummary: Term list to be summarized with the application of each rule.
	:param dict geneSetsDict: Dictionary mapping gene set (term) IDs to set of genes.
	:param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent terms smaller than this size.
	:param function process: The rule to be applied.
	"""

	#Starting with the top terms, term pairs are picked from the termSummary
	#and the rule is applied to them.
	for idNo in range(len(termSummary)-1):
		gsID=termSummary[idNo][0]
		if gsID!=-1:#Check if the term is still a representative term
			for idNo2 in range(idNo+1, len(termSummary)):
				gsID2=termSummary[idNo2][0]
				if  gsID2!=-1:#Check if the term is still a representative term
					termSummary=process(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2)

	#Remove terms that are represented by other terms
	termSummary=[e for e in termSummary if e[0]!=-1]
	#Sort termSummary by rank (first term has the best/lowest rank)
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
term A. Term B's gsID which is stored in termSummary[idNoB][0] is set to -1
to mark that it is not a representative term any more. gsID information is
not lost because it is already copied under term A (termSummary[idNoA][1]).
'''

#In combination methods sending gsIDs as parameters is important
#because when a method assigns -1 in place of gsID, this affects following
#iterations in the calling loop when gsID is obtained by termSummary[idNo][0]



def recurringTermsUnified(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Recurring terms coming from multiple lists are unified.
	This rule is run by default if there are multiple enrichment results.
	This rule is hidden from the user.
	'''
	if(gsID==gsID2):
		#Terms represented by the second term are copied under the first term
		for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
			if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
				termSummary[idNo][1].append(reprTermsOfCoveredTerm)
		termSummary[idNo2][0]=-1
		termSummary[idNo][2]=min(termSummary[idNo][2], termSummary[idNo2][2])
	return termSummary



def supertermRepresentsLessSignificantSubterm(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Superterms with more significance represent their subterms with
	less significance.
	The rule also works for the terms with the same list of genes.
	'''
	geneSet1=geneSetsDict[gsID]#Supposed to be superset and representative
	geneSet2=geneSetsDict[gsID2]#Supposed to be subset and be represented
	if(len(geneSet1)<=maxRepresentativeTermSize):
		if geneSet1.issuperset(geneSet2):
			#Terms represented by the second term are copied under the first term
			for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
				if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
					termSummary[idNo][1].append(reprTermsOfCoveredTerm)
			termSummary[idNo2][0]=-1
	return termSummary



def subtermRepresentsLessSignificantSimilarSuperterm(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Subterms with more significance represent their superterms whose gene list
	is at least 75% composed of the subterm
	'''

	geneSet1=geneSetsDict[gsID]#Supposed to be subset and representative
	geneSet2=geneSetsDict[gsID2]#Supposed to be superset and be represented

	#Representative term size is checked as usual.
	#Represented term size is also checked as it is larger.
	if(len(geneSet1)<=maxRepresentativeTermSize and len(geneSet2)<=maxRepresentativeTermSize):
		if geneSet2.issuperset(geneSet1) and len(geneSet1)/len(geneSet2)>0.75:
			#Terms represented by the second term are copied under the first term
			for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
				if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
					termSummary[idNo][1].append(reprTermsOfCoveredTerm)
			termSummary[idNo2][0]=-1
	return termSummary



def subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Subterms with more significance represent their superterms with less
	significance and less representative power.
	Representative power is the number of terms they represent.
	'''
	geneSet1=geneSetsDict[gsID]#Supposed to be subset and representative
	geneSet2=geneSetsDict[gsID2]#Supposed to be superset and be represented

	#Representative term size is checked as usual.
	#Represented term size is also checked as it is larger.
	if(len(geneSet1)<=maxRepresentativeTermSize and len(geneSet2)<=maxRepresentativeTermSize):
		if geneSet1.issubset(geneSet2):
			#Representative power comparison
			if(len(termSummary[idNo][1])>len(termSummary[idNo2][1])):
				#Terms represented by the second term are copied under the first term
				for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
					if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
						termSummary[idNo][1].append(reprTermsOfCoveredTerm)
				termSummary[idNo2][0]=-1
	return termSummary



def commonSupertermInListRepresentsSubtermsWithLessRepresentativePower(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms that have a common superterm that has more or equal representative
	power than them are represented by this superterm
	'''
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]

	if len(geneSet1)<=maxRepresentativeTermSize and len(geneSet2)<=maxRepresentativeTermSize:

		#All representative terms in termSummary are checked excluding
		#the considered two terms
		for ts in termSummary:
			if ts[0]!=-1 and ts[0]!=gsID and ts[0]!=gsID2:
				geneSetOther=geneSetsDict[ts[0]]
				if(len(geneSetOther)<=maxRepresentativeTermSize):
					if geneSetOther.issuperset(geneSet1) and geneSetOther.issuperset(geneSet2):
						#Representative power comparison
						if(len(ts[1])>len(termSummary[idNo][1]) and len(ts[1])>len(termSummary[idNo2][1])):

							#The term to represent the considered two terms might have better or worse rank.
							#Best (minimum) rank among these three terms will be its new rank.
							#In the rules that more significant covers less significant this code
							#is not needed because the term to be representative has already better rank.
							ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])


							#Terms represented by the considered terms are copied under the representative term

							for reprTermsOfCoveredTerm in termSummary[idNo][1]:
								if reprTermsOfCoveredTerm not in ts[1]:
									ts[1].append(reprTermsOfCoveredTerm)
							termSummary[idNo][0]=-1


							for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
								if reprTermsOfCoveredTerm not in ts[1]:
									ts[1].append(reprTermsOfCoveredTerm)
							termSummary[idNo2][0]=-1

							break

	return termSummary




def supertermRepresentsSubtermLargerThanMaxRep(termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	This rule is to collect large terms into groups, it only works for terms
	that are larger than maxRepresentativeTermSize.
	Superterms larger than maxRepresentativeTermSize represent their less
	significant subterms which are also larger than maxRepresentativeTermSize
	'''
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]
	if len(geneSet1)>maxRepresentativeTermSize and len(geneSet2)>maxRepresentativeTermSize:
		if geneSet1.issuperset(geneSet2):
			for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
				if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
					termSummary[idNo][1].append(reprTermsOfCoveredTerm)
			termSummary[idNo2][0]=-1
	return termSummary




##############################################################################
##############################################################################
##############################################################################



def readGmtFile(gmtPath):
	"""
	Read GMT file. In GMT file each line consists of gene set ID, gene set name and genes,
	all tab separated, no header.

	:param str gmtPath: Path of the GMT file
	"""

	geneSetsDict=dict()#gsID to set of genes mapping
	gsIDToGsNameDict=dict()#gsID to gene set name mapping
	try:
		f=open(gmtPath, 'r')
		lines=f.readlines()
		for line in lines:
			tokens=line.strip().split('\t')
			gsID=tokens[0]
			gsName=tokens[1]
			genes=set(tokens[2:])
			geneSetsDict[gsID]=genes
			gsIDToGsNameDict[gsID]=gsName
		f.close()
	except IOError:
		print("I/O error while reading gmt file.")
	return geneSetsDict, gsIDToGsNameDict



##############################################################################
##############################################################################
##############################################################################

def countTotalRepresentedTermNumber(termSummary):
	tn=0
	for ts in termSummary:
		tn=tn+len(ts[1])
	return tn

##############################################################################
##############################################################################
##############################################################################


def readTBSListFile(tbsFile):
	"""
	Read the input enrichment result file ("to be summarized" file).
	Each line consists of a gene set ID, sorted from most important to least

	:param str tbsFile: Path of the enrichment result files
	"""
	tbsGsIDs=[]
	try:
		f=open(tbsFile, 'r')
		lines=f.readlines()
		for line in lines:
			tbsGsIDs.append(line.strip())
		f.close()
	except IOError:
		print("I/O error while reading a file to be summarized.")

	return tbsGsIDs


def removeUnknownTermsFromTBS(tbsGsIDs, geneSetsDict):
	"""
	Remove unknown gene sets

	:param str tbsFile: Path of the enrichment result files
	:param dict geneSetsDict: Dictionary mapping gene set (term) IDs to set of genes.
	"""

	tbsGsIDs=tbsGsIDs.copy()

	idsToRemove=set()
	for tbsGsID in tbsGsIDs:
		if not(tbsGsID in geneSetsDict.keys()):
			idsToRemove.add(tbsGsID)
	for idToRemove in idsToRemove:
			print(idToRemove, 'is not in gmt file')
			tbsGsIDs.remove(idToRemove)
	return tbsGsIDs


def removeTermsSmallerThanMinTermSize(tbsGsIDs, geneSetsDict, minTermSize):
	"""
	Remove terms smaller than minTermSize

	:param str tbsFile: Path of the enrichment result files
	:param dict geneSetsDict: Dictionary mapping gene set (term) IDs to set of genes.
	:param int minTermSize: The minimum size of the terms to be processed. Smaller terms are discarded.
	"""
	tbsGsIDs=tbsGsIDs.copy()

	idsToRemove=set()
	for tbsGsID in tbsGsIDs:
		if len(geneSetsDict[tbsGsID])<minTermSize:
			idsToRemove.add(tbsGsID)
	for idToRemove in idsToRemove:
			tbsGsIDs.remove(idToRemove)
	return tbsGsIDs

##############################################################################
##############################################################################
##############################################################################


def initializeTermSummary(tbsGsIDsList):
	"""
	Initialize the term summary.
	All terms represent themselves.
	Term summaries are stored in a list.
	Each member contains gene set (term) ID, list of represented terms,
	rank (currently rank of the best term).
	rank starts from 1

	:param list tbsGsIDsList: List of terms obtained from enrichment results files
	"""
	termSummary=[]
	for tbsGsIDsNo in range(len(tbsGsIDsList)):
		tbsGsIDs=tbsGsIDsList[tbsGsIDsNo]
		for idNo in range(len(tbsGsIDs)):
			gsID=tbsGsIDs[idNo]
			rank=idNo+1
			termSummary.append([gsID, [gsID], rank])
	termSummary.sort(key=lambda x: x[2])#Ascending sort based on rank/score
	return termSummary



def writeTermSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, termSummaryFile, termSummaryFile2):
	'''
	Writes the results.
	'''
	try:
		#Detailed
		f=open(termSummaryFile, 'w')

		f.write('Representing term id\tRepresenting term name\tRepresenting term rank')
		for tbsGsIDsNo in range(len(tbsGsIDsList)):
			f.write('\t'+os.path.basename(tbsFiles[tbsGsIDsNo])+ ' term id')
			f.write('\t'+os.path.basename(tbsFiles[tbsGsIDsNo])+ ' term name')
			f.write('\t'+os.path.basename(tbsFiles[tbsGsIDsNo])+ ' term rank')
		f.write('\n')

		for ts in termSummary:
			f.write(ts[0]+'\t'+gsIDToGsNameDict[ts[0]]+'\t'+str(ts[2])+'\n')

			mtr=np.empty([len(ts[1]), len(tbsGsIDsList)*3],dtype=(np.unicode_, 10000))
			row=0
			for reprT in ts[1]:
				for tbsGsIDsNo in range(len(tbsGsIDsList)):
					if reprT in tbsGsIDsList[tbsGsIDsNo]:
						mtr[row, tbsGsIDsNo*3]=reprT
						mtr[row, tbsGsIDsNo*3+1]=gsIDToGsNameDict[reprT]
						mtr[row, tbsGsIDsNo*3+2]=tbsGsIDsList[tbsGsIDsNo].index(reprT)+1
					else:
						mtr[row, tbsGsIDsNo*3]=''
						mtr[row, tbsGsIDsNo*3+1]=''
						mtr[row, tbsGsIDsNo*3+2]=''
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
		for tbsGsIDsNo in range(len(tbsGsIDsList)):
			f.write('\t'+os.path.basename(tbsFiles[tbsGsIDsNo])+ ' term rank')
		f.write('\n')

		for ts in termSummary:
			f.write(ts[0]+'\t'+gsIDToGsNameDict[ts[0]]+'\t'+str(len(geneSetsDict[ts[0]]))+'\t'+str(ts[2])+'\t'+str(len(ts[1])))

			for tbsGsIDsNo in range(len(tbsGsIDsList)):
				found=None
				for reprT in ts[1]:
					try:
						rank=tbsGsIDsList[tbsGsIDsNo].index(reprT)
						if found==None or rank<found:
							found=rank
					except ValueError:
						pass
				if found!=None:
					found=found+1
				f.write('\t'+str(found))#rank starts from 1
			f.write('\n')
		f.close()


	except IOError:
		print("I/O error while writing term summary file.")




def writeHTMLSummaryFile(termSummary, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, termSummaryFile):
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
			#f.write(getTextForTSElement(ts, gsIDToGsNameDict))
			f.write(getTextForTSElementMultiEnrichment(ts, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles))

		f.write('</body>\n')
		f.write('</html>\n')

		f.close()



	except IOError:
		print("I/O error while writing term summary file.")




def getTextForTSElementMultiEnrichment(ts, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles):
	txt='\n'
	txt=txt+'<details>'+'\n'
	#txt=txt+'<summary>'+ts[0]+' '+gsIDToGsNameDict[ts[0]]+' '+str(ts[2])+'</summary>'+'\n'
	txt=txt+'<summary>'+ts[0]+' '+gsIDToGsNameDict[ts[0]]+'</summary>'+'\n'

	for tbsGsIDsNo in range(len(tbsGsIDsList)):
		txt=txt+'\t'+'<p style="margin-left:40px">'+'\n'
		if(len(tbsFiles)>1):
			txt=txt+'\t'+os.path.basename(tbsFiles[tbsGsIDsNo])+'<br>'+'\n'
		for reprT in ts[1]:
			if reprT in tbsGsIDsList[tbsGsIDsNo]:
				txt=txt+'\t'+reprT+' '+gsIDToGsNameDict[reprT]+' (rank: '+str(tbsGsIDsList[tbsGsIDsNo].index(reprT)+1)+', term size: '+ str(len(geneSetsDict[reprT])) +')<br>'+'\n'
		txt=txt+'\t'+'<br>'+'\n'
		txt=txt+'\t'+'</p>'+'\n'
	txt=txt+'</details>'+'\n'
	return txt


