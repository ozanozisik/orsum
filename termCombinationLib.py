#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: Ozan

Functions for filtering
"""

import os

import numpy as np

##############################################################################


def applyRule(termSummary, geneSetsDict, maxRepresentativeTermSize, process):
    """
    This function applies the specified rule ("process") on each pair of terms.

    :param list termSummary: Representative term list to be summarized with the application of rules
    :param dict geneSetsDict: Dictionary mapping term IDs to set of genes.
    :param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
    :param function process: The rule to be applied.
    :return: **termSummary** (*list*) – Representative term list after applying the rule
    """

    # Starting with the top terms, term pairs are picked from the termSummary
    # and the rule is applied to them.
    for idNo in range(len(termSummary) - 1):
        gsID = termSummary[idNo][0]
        if gsID != -1:  # Check if the term is still a representative term
            for idNo2 in range(idNo + 1, len(termSummary)):
                gsID2 = termSummary[idNo2][0]
                if gsID2 != -1:  # Check if the term is still a representative term
                    termSummary = process(
                        termSummary,
                        geneSetsDict,
                        maxRepresentativeTermSize,
                        idNo,
                        idNo2,
                        gsID,
                        gsID2,
                    )

    # Remove terms that are represented by other terms
    termSummary = [e for e in termSummary if e[0] != -1]
    # Sort termSummary by rank (first term has the best/smallest rank)
    termSummary.sort(key=lambda x: x[2])

    return termSummary


##############################################################################
##############################################################################
##############################################################################

"""
General description on how rules work:

Rules take information on two terms from the termSummary.
These two terms are evaluated whether they satisfy the conditions of the rule.
If term A will represent term B, the terms represented by term B
(this includes itself) are appended to the list of terms represented by
term A. Term B's gsID which is stored in termSummary[idNoB][0] is set to -1
to mark that it is not a representative term any more. gsID information is
not lost because it is already copied under term A (termSummary[idNoA][1]).
"""


def recurringTermsUnified(
    termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2
):
    """
    Recurring terms coming from multiple lists are unified.
    This rule is run by default if there are multiple enrichment results.

    :param list termSummary: Representative term list to be summarized
    :param dict geneSetsDict: Dictionary mapping term IDs to set of genes.
    :param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
    :param int idNo: Index of the first term handled in this function in termSummary
    :param int idNo: Index of the second term handled in this function in termSummary
    :param str gsID: Term ID of the first term handled in this function
    :param str gsID2: Term ID of the second term handled in this function
    :return: **termSummary** (*list*) – Representative term list after applying the rule
    """

    if gsID == gsID2:
        # Terms represented by the second term are copied under the first term
        for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
            if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
                termSummary[idNo][1].append(reprTermsOfCoveredTerm)
        termSummary[idNo2][0] = -1
        # Minimum of the ranks is assigned as the rank for this representative term
        termSummary[idNo][2] = min(termSummary[idNo][2], termSummary[idNo2][2])
    return termSummary


def supertermRepresentsLessSignificantSubterm(
    termSummary, geneSetsDict, maxRepresentativeTermSize, idNo, idNo2, gsID, gsID2
):
    """
    Superterms represent their subterms that are less significant.
    The rule also works for the terms with the same list of genes.

    :param list termSummary: Representative term list to be summarized
    :param dict geneSetsDict: Dictionary mapping term IDs to set of genes.
    :param int maxRepresentativeTermSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms.
    :param int idNo: Index of the first term handled in this function in termSummary
    :param int idNo: Index of the second term handled in this function in termSummary
    :param str gsID: Term ID of the first term handled in this function
    :param str gsID2: Term ID of the second term handled in this function
    :return: **termSummary** (*list*) – Representative term list after applying the rule
    """

    geneSet1 = geneSetsDict[gsID]  # Supposed to be superset and representative
    geneSet2 = geneSetsDict[gsID2]  # Supposed to be subset and be represented
    if len(geneSet1) <= maxRepresentativeTermSize:
        if geneSet1.issuperset(geneSet2):
            # Terms represented by the second term are copied under the first term
            for reprTermsOfCoveredTerm in termSummary[idNo2][1]:
                if reprTermsOfCoveredTerm not in termSummary[idNo][1]:
                    termSummary[idNo][1].append(reprTermsOfCoveredTerm)
            termSummary[idNo2][0] = -1
    return termSummary


##############################################################################
##############################################################################
##############################################################################


def readGmtFile(gmtPath):
    """
    Read GMT file. In GMT file each line consists of term ID, term name and genes,
    all tab separated, no header.

    :param str gmtPath: Path of the GMT file
    :return: **geneSetsDict** (*dict*) – Dictionary mapping term IDs to set of genes.
    :return: **gsIDToGsNameDict** (*dict*) – Dictionary mapping term IDs to term names.
    """

    geneSetsDict = dict()  # term ID to set of genes mapping
    gsIDToGsNameDict = dict()  # term ID to term name mapping
    try:
        f = open(gmtPath, "r")
        lines = f.readlines()
        for line in lines:
            tokens = line.strip().split("\t")
            gsID = tokens[0]
            gsName = tokens[1]
            genes = set(tokens[2:])
            geneSetsDict[gsID] = genes
            gsIDToGsNameDict[gsID] = gsName
        f.close()
    except IOError:
        print("I/O error while reading gmt file.")
    return geneSetsDict, gsIDToGsNameDict


##############################################################################
##############################################################################
##############################################################################


def readTBSListFile(tbsFile):
    """
    Read the input enrichment result file ("to be summarized" file).
    Each line consists of a term ID, sorted from the most significant to the least

    :param str tbsFile: Path of the enrichment result files
    :return: **tbsGsIDs** (*list*) – List of term IDs to be summarized
    """
    tbsGsIDs = []
    try:
        f = open(tbsFile, "r")
        lines = f.readlines()
        for line in lines:
            tbsGsIDs.append(line.strip().strip('"'))
        f.close()
    except IOError:
        print("I/O error while reading a file to be summarized.")
        exit()

    return tbsGsIDs


def removeUnknownTermsFromTBS(tbsGsIDs, geneSetsDict):
    """
    Remove unknown terms

    :param list tbsGsIDs: Term IDs list
    :param dict geneSetsDict: Dictionary mapping term IDs to set of genes.
    :return: **tbsGsIDs** (*list*) – Term IDs list after removal of unknown terms
    """

    tbsGsIDs = tbsGsIDs.copy()

    idsToRemove = set()
    for tbsGsID in tbsGsIDs:
        if not (tbsGsID in geneSetsDict.keys()):
            idsToRemove.add(tbsGsID)
    for idToRemove in idsToRemove:
        print(idToRemove, "is not in gmt file")
        tbsGsIDs.remove(idToRemove)
    return tbsGsIDs


def removeTermsSmallerThanMinTermSize(tbsGsIDs, geneSetsDict, minTermSize):
    """
    Remove terms smaller than minTermSize

    :param list tbsGsIDs: Term IDs list
    :param dict geneSetsDict: Dictionary mapping term IDs to set of genes.
    :param int minTermSize: The minimum size of the terms to be processed. Smaller terms are discarded.
    :return: **tbsGsIDs** (*list*) – Term IDs list after removal of small terms
    """
    tbsGsIDs = tbsGsIDs.copy()

    idsToRemove = set()
    for tbsGsID in tbsGsIDs:
        if len(geneSetsDict[tbsGsID]) < minTermSize:
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
    Term summary is a list where each element is also a list describing
    the representative term.
    Each element in term summary contains
    term ID, list of represented terms, rank.
    Initially all terms are representative terms and represent themselves.

    :param list tbsGsIDsList: List of terms obtained from enrichment results files
    :return: **termSummary** (*list*) – Initial representative term list
    """

    termSummary = []
    for tbsGsIDsNo in range(len(tbsGsIDsList)):
        tbsGsIDs = tbsGsIDsList[tbsGsIDsNo]
        for idNo in range(len(tbsGsIDs)):
            gsID = tbsGsIDs[idNo]
            rank = idNo + 1
            termSummary.append([gsID, [gsID], rank])
    termSummary.sort(key=lambda x: x[2])  # Ascending sort based on rank/score
    return termSummary


##############################################################################
##############################################################################
##############################################################################


def writeTermSummaryFile(
    termSummary,
    geneSetsDict,
    gsIDToGsNameDict,
    tbsGsIDsList,
    tbsFiles,
    fileAliases,
    termSummaryFile,
    termSummaryFile2,
):
    """
    Writes the results.
    """
    try:
        # Detailed
        f = open(termSummaryFile, "w")

        f.write("Representing term id\tRepresenting term name\tRepresenting term rank")
        for tbsGsIDsNo in range(len(tbsGsIDsList)):
            f.write("\t" + fileAliases[tbsGsIDsNo] + " term id")
            f.write("\t" + fileAliases[tbsGsIDsNo] + " term name")
            f.write("\t" + fileAliases[tbsGsIDsNo] + " term rank")
        f.write("\n")

        for ts in termSummary:
            f.write(ts[0] + "\t" + gsIDToGsNameDict[ts[0]] + "\t" + str(ts[2]) + "\n")

            mtr = np.empty(
                [len(ts[1]), len(tbsGsIDsList) * 3], dtype=(np.unicode_, 10000)
            )
            row = 0
            for reprT in ts[1]:
                for tbsGsIDsNo in range(len(tbsGsIDsList)):
                    if reprT in tbsGsIDsList[tbsGsIDsNo]:
                        mtr[row, tbsGsIDsNo * 3] = reprT
                        mtr[row, tbsGsIDsNo * 3 + 1] = gsIDToGsNameDict[reprT]
                        mtr[row, tbsGsIDsNo * 3 + 2] = (
                            tbsGsIDsList[tbsGsIDsNo].index(reprT) + 1
                        )
                    else:
                        mtr[row, tbsGsIDsNo * 3] = ""
                        mtr[row, tbsGsIDsNo * 3 + 1] = ""
                        mtr[row, tbsGsIDsNo * 3 + 2] = ""
                row = row + 1

            for r in range(row):
                f.write("\t\t")
                for c in range(mtr.shape[1]):
                    f.write("\t" + mtr[r, c])
                f.write("\n")
        f.close()

        # Summary of summary
        f = open(termSummaryFile2, "w")

        f.write(
            "Representing term id\tRepresenting term name\tRepresenting term size\tRepresenting term rank\tRepresented term number"
        )
        for tbsGsIDsNo in range(len(tbsGsIDsList)):
            f.write("\t" + fileAliases[tbsGsIDsNo] + " term rank")
        f.write("\n")

        for ts in termSummary:
            f.write(
                ts[0]
                + "\t"
                + gsIDToGsNameDict[ts[0]]
                + "\t"
                + str(len(geneSetsDict[ts[0]]))
                + "\t"
                + str(ts[2])
                + "\t"
                + str(len(ts[1]))
            )

            for tbsGsIDsNo in range(len(tbsGsIDsList)):
                found = None
                for reprT in ts[1]:
                    try:
                        rank = tbsGsIDsList[tbsGsIDsNo].index(reprT)
                        if found == None or rank < found:
                            found = rank
                    except ValueError:
                        pass
                if found != None:
                    found = found + 1
                f.write("\t" + str(found))  # rank starts from 1
            f.write("\n")
        f.close()

    except IOError:
        print("I/O error while writing term summary file.")


def writeHTMLSummaryFile(
    termSummary,
    geneSetsDict,
    gsIDToGsNameDict,
    tbsGsIDsList,
    tbsFiles,
    fileAliases,
    termSummaryFile,
):
    """
    Writes the results to an HTML file
    """
    try:
        # Detailed
        f = open(termSummaryFile, "w")
        f.write("<!DOCTYPE html>\n")
        f.write("<html>\n")
        f.write("<head>\n")
        f.write("<title>orsum result</title>\n")
        f.write("</head>\n")
        f.write("<body>\n")

        for ts in termSummary:
            # f.write(getTextForTSElement(ts, gsIDToGsNameDict))
            f.write(
                getTextForTSElementMultiEnrichment(
                    ts,
                    geneSetsDict,
                    gsIDToGsNameDict,
                    tbsGsIDsList,
                    tbsFiles,
                    fileAliases,
                )
            )

        f.write("</body>\n")
        f.write("</html>\n")

        f.close()

    except IOError:
        print("I/O error while writing term summary file.")


def writeRepresentativeToRepresentedIDsFile(termSummary, outputFile):
    """
    Writes ID mapping between representative terms and represented terms
    """
    f = open(outputFile, "w")
    f.write("Representative\tRepresented")
    for ts in termSummary:
        representativeID = ts[0]
        for representedID in ts[1]:
            f.write("\n")
            f.write(representativeID + "\t" + representedID)
    f.close()


def getTextForTSElementMultiEnrichment(
    ts, geneSetsDict, gsIDToGsNameDict, tbsGsIDsList, tbsFiles, fileAliases
):
    txt = "\n"
    txt = txt + "<details>" + "\n"
    # txt=txt+'<summary>'+ts[0]+' '+gsIDToGsNameDict[ts[0]]+' '+str(ts[2])+'</summary>'+'\n'
    txt = (
        txt + "<summary>" + ts[0] + " " + gsIDToGsNameDict[ts[0]] + "</summary>" + "\n"
    )

    for tbsGsIDsNo in range(len(tbsGsIDsList)):
        txt = txt + "\t" + '<p style="margin-left:40px">' + "\n"
        if len(fileAliases) > 1:
            txt = txt + "\t" + fileAliases[tbsGsIDsNo] + "<br>" + "\n"
        for reprT in ts[1]:
            if reprT in tbsGsIDsList[tbsGsIDsNo]:
                txt = (
                    txt
                    + "\t"
                    + reprT
                    + " "
                    + gsIDToGsNameDict[reprT]
                    + " (rank: "
                    + str(tbsGsIDsList[tbsGsIDsNo].index(reprT) + 1)
                    + ", term size: "
                    + str(len(geneSetsDict[reprT]))
                    + ")<br>"
                    + "\n"
                )
        txt = txt + "\t" + "<br>" + "\n"
        txt = txt + "\t" + "</p>" + "\n"
    txt = txt + "</details>" + "\n"
    return txt
