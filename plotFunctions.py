#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Created on April 2021.

@author: Morgane T.

Create plots for orsum results
"""

# LIBRARIES
import os
import numpy as np
import seaborn as sns
import pandas as ps
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# FUNCTIONS

def orsum_readResultFile(inputFile):
	"""
	Read and parse result file created by orsum.py.

	Read results file from orsum and parse it.
	Create a dataframe with the main results.
	Create an array with the results for each analysis.

	:param str inputFile: Path of the result file created by orsum.py

	:returns:
		- **df** (*pandas.DataFrame*) – Data frame with 4 columns (sizes, labels, ranks, colors)
		- **allRanks_array** (*numpy.ndarray*) – Array with the rank for each analysis contains in the results file
		- **resultsId** (*list*) – List of analysis names
		- **palette_cmap** (*matplotlib.ListedColormap*) – Color map
	"""
	# PARAMETERS
	sizes = []
	labels = []
	ranks = []
	allRanks = []
	resultsId = []

	# READ INPUT FILE FROM ORSUM
	try:
		with open(inputFile, 'r') as resultsFileHandler:
			headerLine = resultsFileHandler.readline()
			lheaderLine = headerLine.rstrip('\n').split('\t')
			for r in range(5, len(lheaderLine)):
				if lheaderLine[r].endswith(' term rank'):
					lheaderLine[r] = lheaderLine[r][:-10]
				resultsId.append(lheaderLine[r])
			for line in resultsFileHandler:
				lLine = line.rstrip('\n').split('\t')
				indRanks = []
				for r in range(5, len(lLine)):
					if(lLine[r] == 'None'):
						indRanks.append(np.nan)
					else:
						indRanks.append(int(lLine[r]))
				sizes.append(int(lLine[4]))
				labels.append(lLine[1])
				ranks.append(int(lLine[3]))
				indRanks_array = np.array(indRanks)
				allRanks.append(indRanks_array)
	except IOError:
		print("I/O error while reading result file.")

	# SELECT COLOR
	rankMax = max(ranks)
	palette = sns.cubehelix_palette(light = .8, n_colors = rankMax, as_cmap = False, reverse = True)
	palette_cmap = sns.cubehelix_palette(light = .8, n_colors = rankMax, as_cmap = True, reverse = True)
	ranksColor = []
	for rank in ranks:
		ranksColor.append(palette[rank - 1])

	# CREATE DATA FRAME
	df = ps.DataFrame({'sizes': sizes, 'labels': labels, 'ranks': ranks, 'colors': ranksColor})

	# CREATE ARRAY WITH ALL RESULTS
	allRanks_array = np.array(allRanks)

	return(df, allRanks_array, resultsId, palette_cmap)


def orsum_barplot(df, nbTerm, sizeMax, sizeMin, plotName, ticks):
	"""
	Create ans save barplot of the main results from orsum.py.

	:param pandas.DataFrame df: Data frame with 4 columns (sizes, labels, ranks, colors). Created by orsum_readResultFile() function.
	:param int nbTerm: Number of top results you want to display (MAX = 50).
	:param int sizeMax: Size max of representing term.
	:param int sizeMin: Size min of representing term.
	:param str plotName: Path and name of the plot created by this function.
	:param list ticks: List of integer. Correspond of the values for the colorbar ticks.
	"""
	# Select the nbTerm top of results
	df_filt = df[0:nbTerm]
	# Create theme and plot
	sns.set_theme(style = 'whitegrid')
	sns.set_context('paper', font_scale = 1.1)
	plt.figure(figsize = (10, 10))
	ax = sns.barplot(data = df_filt, x = 'sizes', y = 'labels', edgecolor = 'white', color = df_filt['colors'][1])
	# Labels
	plt.title('Top {} representative terms'.format(nbTerm), fontsize = 20)
	plt.xlabel('Number of represented terms', fontsize = 15)
	plt.ylabel('')
	plt.xlim(0, max(df_filt['sizes']))
	# Save and close plot
	plt.savefig(plotName, bbox_inches = 'tight', dpi = 300)

def orsum_heatmap(allRanks_array, df, nbTerm, plotName, conditionName, palette_cmap, ticks):
	"""
	Create and save heatmap of the results of each analysis from orsum.py.
	Field with ranks from each analysis.

	:param numpy.ndarray allRanks_array: Array with the rank for each analysis contains in the results file. Created by orsum_readResultFile() function.
	:param pandas.DataFrame df: Data frame with 4 columns (sizes, labels, ranks, colors). Created by orsum_readResultFile() function.
	:param int nbTerm: Number of top results you want to display (MAX = 50).
	:param str plotName: Path and name of the plot created by this function.
	:param list conditionName: List of analysis names.
	:param matplotlib.ListedColormap palette_cmap: Color map.
	:param list ticks: List of integer. Correspond of the values for the colorbar ticks.
	"""
	# Select the nbTerm top of results
	array = allRanks_array[0:nbTerm]
	df_filt = df[0:nbTerm]
	# Create theme and plot
	sns.set(font_scale = 0.5)
	plt.figure(figsize = (10, 6))
	plt.title('Representative term ranks', fontsize = 10)
	ax = sns.heatmap(array, cmap = palette_cmap, square = True, linewidths = 0.5, yticklabels = df_filt['labels'], cbar_kws = {'shrink': 0.5, 'ticks': ticks}, xticklabels = conditionName, vmin = df['ranks'].min(), vmax = df['ranks'].max())
	# Colorbar
	ax.collections[0].colorbar.set_label("Ranks")
	ax.collections[0].colorbar.ax.set_ylim(df['ranks'].max(), 0)
	# Save and close plot
	plt.savefig(plotName, bbox_inches = 'tight', dpi = 300)

def calculateQuartileFromRanks(allRanks_array):
	"""
	Calculate the quartiles values for each condition
	Create a new array with the number of quartile associated for each rank

	:param numpy.ndarray allRanks_array: Array with the rank for each analysis contains in the results file
	:return: **allQ** (*numpy.ndarray*) – Array with the quartile for each analysis
	"""
	# Calcul quartiles for each condition
	allRanks_df = ps.DataFrame(allRanks_array)
	Q_array = allRanks_df.quantile(q=[0.25, 0.50, 0.75])
	# Create new array with quartiles
	allQ = []
	for line in allRanks_array:
		quantiles = []
		for ind in range(len(line)):
			rank = line[ind]
			if (np.isnan(rank)):
				quantiles.append(np.nan)
			else:
				if rank <= Q_array[ind][0.25]:
					quantiles.append(1)
				if Q_array[ind][0.25] < rank <= Q_array[ind][0.50]:
					quantiles.append(2)
				if Q_array[ind][0.50] < rank <= Q_array[ind][0.75]:
					quantiles.append(3)
				if rank > Q_array[ind][0.75]:
					quantiles.append(4)
		quantiles_array = np.array(quantiles)
		allQ.append(quantiles_array)
	return(allQ)


def orsum_heatmap_quartile_quantitative(quartiles_array, df, nbTerm, plotName, conditionName):
	"""
	Create and save heatmap of the results of each analysis from orsum.py.
	Display quartile calculated from each condition
	Discrete color map

	:param numpy.ndarray quartiles_array: Array with the quartile for each analysis contains in the results file. Created by orsum_readResultFile() function.
	:param pandas.DataFrame df: Data frame with 4 columns (sizes, labels, ranks, colors). Created by orsum_readResultFile() function.
	:param int nbTerm: Number of top results you want to display (MAX = 50).
	:param str plotName: Path and name of the plot created by this function.
	:param list conditionName: List of analysis names.
	"""
	# Color for heatmap
	nbMax = 4
	quartilePalette_cmap = sns.cubehelix_palette(light=.8, n_colors=nbMax, as_cmap=True, reverse=True)
	bounds = [1, 2, 3, 4, 5]
	norm = colors.BoundaryNorm(bounds, quartilePalette_cmap.N)

	# Select the nbTerm top of results
	array = quartiles_array[0:nbTerm]
	df_filt = df[0:nbTerm]
	# Create theme and plot
	sns.set(font_scale = 0.5)
	plt.figure(figsize = (10, 6))
	plt.title('Representative term ranks', fontsize = 10)
	ax = sns.heatmap(array, cmap = quartilePalette_cmap, norm = norm, square=True, linewidths=0.5, yticklabels=df_filt['labels'], xticklabels = conditionName)
	# Colorbar
	ax.collections[0].colorbar.ax.set_ylim(5, 1)
	colorbar = ax.collections[0].colorbar
	colorbar.set_ticks([1.5, 2.5, 3.5, 4.5])
	colorbar.set_ticklabels(["Q1", "Q2", "Q3", "Q4"])
	ax.collections[0].colorbar.set_label("Quartiles")
	# Save and close plot
	plt.savefig(plotName, bbox_inches = 'tight', dpi = 300)
	plt.close()

def orsum_linePlot(df, plotName):
	"""
	Create and save scatterplot of the size of represented terms from orsum.py.

	:param pandas.DataFrame df: Data frame with 4 columns (sizes, labels, ranks, colors). Created by orsum_readResultFile() function.
	:param str plotName: Path and name of the plot created by this function.
	"""
	# Create theme and plot
	sns.set_theme(style = 'whitegrid')
	sns.set_context('paper', font_scale = 1.1)
	plt.figure(figsize = (10, 10))
	plt.scatter(df['ranks'], df['sizes'], c = 'purple')
	# Labels
	#plt.title('Size of each representative term', fontsize = 20)
	plt.title('Representative power vs rank of the representative term', fontsize = 20)
	plt.xlabel('Rank', fontsize = 15)
	#plt.ylabel('Number of terms inside the representative term', fontsize = 15)
	plt.ylabel('Number of represented terms', fontsize = 15)
	# Save and close plot
	plt.savefig(plotName, bbox_inches = 'tight', dpi = 300)

def createBoundaries4Colorbar(df, step):
	"""
	Create list of boundaries for colorbar.

	:param pandas.DataFrame df: Data frame with 4 columns (sizes, labels, ranks, colors). Created by orsum_readResultFile() function.
	:param int step: Difference between each number in the sequence.

	:returs: **boundariesCB** (*list*) – list of values.
	"""
	boundariesCB = list(range(0, df['ranks'].max(), step))
	boundariesCB[0] = df['ranks'].min()
	if(df['ranks'].max() - boundariesCB[-1] < step/2):
		boundariesCB[-1] = df['ranks'].max()
	else:
		boundariesCB.append(df['ranks'].max())
	return(boundariesCB)


def orsum_plot(inputFile, outputDir, threshold, heatmapName = 'Heatmap'):
	"""
	Main function.

	Initialisation of the parameters.
	Read and parse results file from orsum.py.
	Calcul bounderies for colorbar.
	Create and save severals plot.

	:param str inputFile: Path name of the orsum results file (termSummaryXX-Summary.tsv).
	:param str outputDir: Folder path name to write the results.
	:param int threshold: Number of top results you want to display (MAX = 50)
	"""
	# PARAMETERS
	barplotName = '{}{}{}'.format(outputDir, os.sep, 'Barplot')
	#heatmapName = '{}{}{}'.format(outputDir, os.sep, 'Heatmap')
	heatmapQuartQuantitativeName = '{}{}{}'.format(outputDir, os.sep, heatmapName)
	lineplotName = '{}{}{}'.format(outputDir, os.sep, 'SizesDistribution')
	# CHECK
	if(threshold > 50):
		print("Number of terms to be plotted was greater than 50, it is changed to 50.")
		threshold = 50
	# READ FILE AND CALCUL BOUNDARIES
	df, allRanks_array, resultsId, palette_cmap = orsum_readResultFile(inputFile = inputFile)
	boundariesCB = createBoundaries4Colorbar(df, step = 100)
	# CALCULATE QUARTILES
	allQuartiles_array = calculateQuartileFromRanks(allRanks_array = allRanks_array)
	# PLOTS
	orsum_linePlot(df = df, plotName = lineplotName)
	orsum_barplot(df = df, nbTerm = threshold, sizeMax = df['ranks'].max(), sizeMin = df['ranks'].min(), plotName = barplotName, ticks = boundariesCB)
	#orsum_heatmap(allRanks_array = allRanks_array, df = df, nbTerm = threshold, plotName = heatmapName, conditionName = resultsId, palette_cmap = palette_cmap, ticks = boundariesCB)
	orsum_heatmap_quartile_quantitative(quartiles_array = allQuartiles_array, df = df, nbTerm = threshold, plotName = heatmapQuartQuantitativeName, conditionName = resultsId)

