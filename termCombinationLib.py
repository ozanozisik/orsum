"""
Created on Mon Nov 30 13:57:07 2020

@author: Ozan
"""

import numpy as np


##############################################################################



def applyRule(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, process):
	'''
	This function is called to apply a rule on term pairs. The function
	of the rule to be applied is given as parameter (process).
	'''

	for idNo in range(len(termSummary)-1):
		gsID=termSummary[idNo][0]
		if gsID!=-1:
			for idNo2 in range(idNo+1, len(termSummary)):
				gsID2=termSummary[idNo2][0]
				if  gsID2!=-1:
					termSummary=process(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2)
	termSummary=[e for e in termSummary if e[0]!=-1]
	termSummary.sort(key=lambda x: x[2])#Sort based on rank/score
	return termSummary


##############################################################################
##############################################################################
##############################################################################

#In combination methods sending gsIDs as parameters is important
#because when a method assigns -1 in place of gsID, this affects following
#iterations in the calling loop when gsID is obtained by termSummary[idNo][0]



def recurringTermsUnified(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms part of more than one list are unified.
	'''
	if((gsID==gsID2) or (geneSetsDict[gsID]==geneSetsDict[gsID2])):
		#This process merges recurring terms that come from multiple
		#enrichment results given as input, or gene sets that contain
		#the same genes.
		for reprTermsOfEaten in termSummary[idNo2][1]:
			if reprTermsOfEaten not in termSummary[idNo][1]:
				termSummary[idNo][1].append(reprTermsOfEaten)
		termSummary[idNo2][0]=-1
		termSummary[idNo][2]=min(termSummary[idNo][2], termSummary[idNo2][2])
	return termSummary



def supertermRepresentsSubterm(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Superterms with more significance represent their subterms with
	less significance.
	'''
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]
	if(len(geneSet1)<=maxRepresentativeTermSize):
		if geneSet1.issuperset(geneSet2):
			for reprTermsOfEaten in termSummary[idNo2][1]:
				if reprTermsOfEaten not in termSummary[idNo][1]:
					termSummary[idNo][1].append(reprTermsOfEaten)
			termSummary[idNo2][0]=-1
	return termSummary



def commonParentInListRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms with a common parent are represented by the parent
	if the parent is in the original list
	'''
	representingTerm=None
	try:
		if(hierarchyDict[gsID]==hierarchyDict[gsID2]):
			if(hierarchyDict[gsID] in originalTermsSet):
				representingTerm=hierarchyDict[gsID]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break
			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])

				for reprTermsOfEaten in termSummary[idNo][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo][0]=-1

				for reprTermsOfEaten in termSummary[idNo2][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo2][0]=-1

	return termSummary



def commonGrandparentInListRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms with a common grandparent are represented by the grandparent
	if the grandparent is in the original list
	'''
	representingTerm=None
	try:
		if(hierarchyDict[hierarchyDict[gsID]]==hierarchyDict[hierarchyDict[gsID2]]):
			if(hierarchyDict[hierarchyDict[gsID]] in originalTermsSet):
				representingTerm=hierarchyDict[hierarchyDict[gsID]]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break
			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])

				for reprTermsOfEaten in termSummary[idNo][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo][0]=-1
				for reprTermsOfEaten in termSummary[idNo2][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo2][0]=-1

	return termSummary



def commonParentGrandparentInListRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms in which one's parent is other's grandparent are represented
	by this ancestor term if it is in the original list
	'''
	representingTerm=None
	try:
		if(hierarchyDict[hierarchyDict[gsID]]==hierarchyDict[gsID2]):
			if(hierarchyDict[gsID2] in originalTermsSet):
				representingTerm=hierarchyDict[gsID2]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem
	try:
		if(hierarchyDict[gsID]==hierarchyDict[hierarchyDict[gsID2]]):
			if(hierarchyDict[gsID] in originalTermsSet):
				representingTerm=hierarchyDict[gsID]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break
			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])

				for reprTermsOfEaten in termSummary[idNo][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo][0]=-1
				for reprTermsOfEaten in termSummary[idNo2][1]:
					if reprTermsOfEaten not in ts[1]:
						ts[1].append(reprTermsOfEaten)
				termSummary[idNo2][0]=-1

	return termSummary



def subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Subterms represent their superterms with less significance / lower rank
	and less representative power.
	Representative power is the number of terms they represent.
	'''
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]
	if(len(geneSet1)<=maxRepresentativeTermSize):
		if geneSet1.issubset(geneSet2):
			if(len(termSummary[idNo][1])>len(termSummary[idNo2][1])):
				for reprTermsOfEaten in termSummary[idNo2][1]:
					if reprTermsOfEaten not in termSummary[idNo][1]:
						termSummary[idNo][1].append(reprTermsOfEaten)
				termSummary[idNo2][0]=-1
	return termSummary



def subtermRepresentsSlightlyLowerRankedSuperterm(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Subterms represent their slightly lower ranked superterms. Rank tolerance is set to 1.
	'''
	rankDiffTolerance=1
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]
	if(len(geneSet1)<=maxRepresentativeTermSize):
		if geneSet2.issuperset(geneSet1):
			if (termSummary[idNo2][2]-termSummary[idNo][2])<=rankDiffTolerance:
				for reprTermsOfEaten in termSummary[idNo2][1]:
					if reprTermsOfEaten not in termSummary[idNo][1]:
						termSummary[idNo][1].append(reprTermsOfEaten)
				termSummary[idNo2][0]=-1
	return termSummary



def commonParentRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms with a common parent are represented by the parent
	(even if the parent is not in the original list)
	'''
	representingTerm=None
	try:
		if(hierarchyDict[gsID]==hierarchyDict[gsID2]):
			representingTerm=hierarchyDict[gsID]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break

			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])
			else:
				termSummary.append([representingTerm, [], min(termSummary[idNo][2], termSummary[idNo2][2])])
				ts=termSummary[-1]

			for reprTermsOfEaten in termSummary[idNo][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo][0]=-1

			for reprTermsOfEaten in termSummary[idNo2][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo2][0]=-1

	return termSummary



def commonGrandparentRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms with a common grandparent are represented by the grandparent
	(even if the grandparent is not in the original list)
	'''
	representingTerm=None
	try:
		if(hierarchyDict[hierarchyDict[gsID]]==hierarchyDict[hierarchyDict[gsID2]]):
			representingTerm=hierarchyDict[hierarchyDict[gsID]]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break

			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])
			else:
				termSummary.append([representingTerm, [], min(termSummary[idNo][2], termSummary[idNo2][2])])
				ts=termSummary[-1]

			for reprTermsOfEaten in termSummary[idNo][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo][0]=-1

			for reprTermsOfEaten in termSummary[idNo2][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo2][0]=-1

	return termSummary



def commonParentGrandparentRepresentsTerms(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Terms in which one's parent is other's grandparent are represented
	by this ancestor term (even if it is not in the original list)
	'''
	representingTerm=None
	try:
		if(hierarchyDict[hierarchyDict[gsID]]==hierarchyDict[gsID2]):
			representingTerm=hierarchyDict[gsID2]
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem
	try:
		if(hierarchyDict[gsID]==hierarchyDict[hierarchyDict[gsID2]]):
			representingTerm=hierarchyDict[gsID]
			#print('Superterm of one is grand superterm of other')
	except KeyError:
		pass#A top term in hierarchy is encountered, no problem

	if representingTerm is not None:
		if(len(geneSetsDict[representingTerm])<=maxRepresentativeTermSize):
			tsFound=False
			for ts in termSummary:
				if ts[0]==representingTerm:
					tsFound=True
					break

			if tsFound:
				ts[2]=min(termSummary[idNo][2], termSummary[idNo2][2], ts[2])
			else:
				termSummary.append([representingTerm, [], min(termSummary[idNo][2], termSummary[idNo2][2])])
				ts=termSummary[-1]

			for reprTermsOfEaten in termSummary[idNo][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo][0]=-1

			for reprTermsOfEaten in termSummary[idNo2][1]:
				if reprTermsOfEaten not in ts[1]:
					ts[1].append(reprTermsOfEaten)
			termSummary[idNo2][0]=-1

	return termSummary






def Template(termSummary, geneSetsDict, hierarchyDict, originalTermsSet, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2):
	'''
	Template for new rules
	'''
	geneSet1=geneSetsDict[gsID]
	geneSet2=geneSetsDict[gsID2]


	return termSummary




##############################################################################
##############################################################################
##############################################################################



def readGmtFile(gmtPath):
	'''
	Reads GMT file.
	Each line consists of gene set ID, gene set name and genes,
	all tab separated, no header.
	'''
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


def createTermHierarchy(geneSetsDict, termHierarchyFile=None):
	'''
	Creates hierarchy between terms based on subset/superset relation.
	Hiearchy can be saved to a file.
	'''
	hierarchyDict=dict() #subset to superset
	for gsID, genes1 in geneSetsDict.items():
		geneSet1=geneSetsDict[gsID]
		for gsID2, genes2 in geneSetsDict.items():
			if gsID!=gsID2:
				geneSet2=geneSetsDict[gsID2]

				#When there is a superset subset relation we check whether
				#a parent was previously assigned in the hierarchy.
				#If yes, we check the size of the existing superset (say
				#superE) and new superset (superN). If superN is smaller,
				#it is assigned as the parent. The removed superset superE
				#is not lost because during its processing step, it is already
				#assigned as parent to its subterms, possibly including
				#superN
				if(geneSet1!=geneSet2):
					if geneSet1.issuperset(geneSet2):
						curSuperset=hierarchyDict.get(gsID2)
						if curSuperset is None:
							hierarchyDict[gsID2]=gsID
						else:
							if len(geneSetsDict[gsID])<len(geneSetsDict[curSuperset]):
								hierarchyDict[gsID2]=gsID
					if geneSet2.issuperset(geneSet1):
						curSuperset=hierarchyDict.get(gsID)
						if curSuperset is None:
							hierarchyDict[gsID]=gsID2
						else:
							if len(geneSetsDict[gsID2])<len(geneSetsDict[curSuperset]):
								hierarchyDict[gsID]=gsID2

	if termHierarchyFile is not None:
		try:
			f=open(termHierarchyFile, 'w')
			for k,v in hierarchyDict.items():
				f.write(k+'\t'+v+'\n')
			f.close()
		except IOError:
		    print("I/O error while writing hierarchy file.")

	return hierarchyDict


def readHierarchyFile(termHierarchyFile):
	'''
	Reads hierarchy file.
	Each line consists of two term IDs separated by tab.
	'''
	hierarchyDict=dict() #subset to superset
	try:
		f=open(termHierarchyFile, 'r')
		lines=f.readlines()
		for line in lines:
			tokens=line.strip().split('\t')
			hierarchyDict[tokens[0]]=tokens[1]
		f.close()
	except IOError:
	    print("I/O error while reading hierarchy file.")
	return hierarchyDict


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
	'''
	Reads the ranked gene set list ("to be summarized" file).
	Each line consists of a gene set ID, it should be sorted from
	most important to least
	'''
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

def removeUnknownGeneSetsFromTBS(tbsGsIDs, geneSetsDict):
	'''
	Removes unknown gene sets
	'''
	idsToRemove=set()
	for tbsGsID in tbsGsIDs:
		if not(tbsGsID in geneSetsDict.keys()):
			idsToRemove.add(tbsGsID)
	for idToRemove in idsToRemove:
			print(idToRemove, 'is not in gmt file')
			tbsGsIDs.remove(idToRemove)
	return tbsGsIDs


##############################################################################
##############################################################################
##############################################################################


def initializeTermSummary(tbsGsIDsList):
	'''
	This function initializes the term summary.
	All terms represent themselves.
	Term summaries are stored in a list.
	Each member contains representing term ID, list of represented terms,
	rank (currently rank of the best term)
	rank starts from 1
	'''
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
			f.write('\t'+tbsFiles[tbsGsIDsNo]+ ' term id')
			f.write('\t'+tbsFiles[tbsGsIDsNo]+ ' term name')
			f.write('\t'+tbsFiles[tbsGsIDsNo]+ ' term rank')
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
			f.write('\t'+tbsFiles[tbsGsIDsNo]+ ' term rank')
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
