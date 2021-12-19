from termCombinationLib import initializeTermSummary, applyRule, recurringTermsUnified, supertermRepresentsLessSignificantSubterm, subtermRepresentsLessSignificantSimilarSuperterm, subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower, commonSupertermInListRepresentsSubtermsWithLessRepresentativePower, supertermRepresentsSubtermLargerThanMaxRep

def test_initializeTermSummary_singleInput():
	tbsGsIDsList=[['term1', 'term2', 'term3']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	assert termSummary==[
		['term1', ['term1'], 1],
		['term2', ['term2'], 2],
		['term3', ['term3'], 3]
		]

def test_initializeTermSummary_multiInput():
	tbsGsIDsList=[['term11', 'termCommon', 'term13', 'term14'], ['term21', 'term22', 'term23', 'termCommon']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	assert termSummary==[
		['term11', ['term11'], 1],
		['term21', ['term21'], 1],
		['termCommon', ['termCommon'], 2],
		['term22', ['term22'], 2],
		['term13', ['term13'], 3],
		['term23', ['term23'], 3],
		['term14', ['term14'], 4],
		['termCommon', ['termCommon'], 4]
		]

def test_rule_recurringTermsUnified():
	tbsGsIDsList=[['term11', 'termCommon', 'term13', 'term14'], ['term21', 'term22', 'term23', 'termCommon']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict=dict()#Not important, not used in this rule
	termSummary=applyRule(termSummary, geneSetsDict, 2000, recurringTermsUnified)
	assert termSummary==[
		['term11', ['term11'], 1],
		['term21', ['term21'], 1],
		['termCommon', ['termCommon'], 2],
		['term22', ['term22'], 2],
		['term13', ['term13'], 3],
		['term23', ['term23'], 3],
		['term14', ['term14'], 4],
		]

def test_rule_supertermRepresentsLessSignificantSubterm():
	tbsGsIDsList=[['term1', 'term2', 'term3', 'term4', 'term5', 'term6']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict={
		'term1':{'A','B','C'},
		'term2':{'A','B','C','D','E','F'},
		'term3':{'A','B','C','D'},
		'term4':{'A','B','G','H'},
		'term5':{'A','B'},
		'term6':{'G','H'}
		}
	termSummary=applyRule(termSummary, geneSetsDict, 2000, supertermRepresentsLessSignificantSubterm)
	assert termSummary==[
		['term1', ['term1', 'term5'], 1],
		['term2', ['term2', 'term3'], 2],
		['term4', ['term4', 'term6'], 4],
		]


def test_rule_subtermRepresentsLessSignificantSimilarSuperterm():
	tbsGsIDsList=[['term1', 'term2', 'term3', 'term4', 'term5', 'term6']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict={
		'term1':{'A','B','C'},
		'term2':{'A','B','C','D','E','F'},
		'term3':{'A','B','C','D'},
		'term4':{'A','B','C','D','E','F','G'},
		'term5':{'A','B'},
		'term6':{'G','H'}
		}
	termSummary=applyRule(termSummary, geneSetsDict, 2000, subtermRepresentsLessSignificantSimilarSuperterm)
	assert termSummary==[
		['term1', ['term1', 'term3'], 1],
		['term2', ['term2', 'term4'], 2],
		['term5', ['term5'], 5],
		['term6', ['term6'], 6]
		]

def test_rule_subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower():
	tbsGsIDsList=[['term1', 'term2', 'term3', 'term4', 'term5']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict={
		'term1':{'A','B','C'},
		'term2':{'G','H'},
		'term3':{'A','B','C','D'},
		'term4':{'A','B'},
		'term5':{'A','B','C','D','E','F'}
		}
	termSummary=applyRule(termSummary, geneSetsDict, 2000, subtermRepresentsLessSignificantSimilarSuperterm)
	assert termSummary==[
		['term1', ['term1', 'term3'], 1],
		['term2', ['term2'], 2],
		['term4', ['term4'], 4],
		['term5', ['term5'], 5]
		]
	termSummary=applyRule(termSummary, geneSetsDict, 2000, subtermRepresentsSupertermWithLessSignificanceAndLessRepresentativePower)
	assert termSummary==[
		['term1', ['term1', 'term3', 'term5'], 1],
		['term2', ['term2'], 2],
		['term4', ['term4'], 4]
		]

def test_rule_commonSupertermInListRepresentsSubtermsWithLessRepresentativePower():
	tbsGsIDsList=[['term1', 'term2', 'term3', 'term4']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict={
		'term1':{'A','B','C'},
		'term2':{'G','H'},
		'term3':{'A','B','C','D','G','H'},
		'term4':{'A','B','D'}
		}
	termSummary=applyRule(termSummary, geneSetsDict, 2000, supertermRepresentsLessSignificantSubterm)
	assert termSummary==[
		['term1', ['term1'], 1],
		['term2', ['term2'], 2],
		['term3', ['term3','term4'], 3],
		]

	termSummary=applyRule(termSummary, geneSetsDict, 2000, commonSupertermInListRepresentsSubtermsWithLessRepresentativePower)
	assert termSummary==[
		['term3', ['term3','term4', 'term1', 'term2'], 1]
		]



def test_rule_supertermRepresentsSubtermLargerThanMaxRep():
	tbsGsIDsList=[['term1', 'term2', 'term3', 'term4', 'term5', 'term6']]
	termSummary=initializeTermSummary(tbsGsIDsList)
	geneSetsDict={
		'term1':{'A','B','C'},
		'term2':{'A','B','C','D','E','F','G','H','I','J'},
		'term3':{'A','B','C','D'},
		'term4':{'A','B','C','D','G','H','I','J'},
		'term5':{'A','B'},
		'term6':{'G','H'}
		}
	termSummary=applyRule(termSummary, geneSetsDict, 6, supertermRepresentsSubtermLargerThanMaxRep)
	assert termSummary==[
		['term1', ['term1'], 1],
		['term2', ['term2', 'term4'], 2],
		['term3', ['term3'], 3],
		['term5', ['term5'], 5],
		['term6', ['term6'], 6],
		]
