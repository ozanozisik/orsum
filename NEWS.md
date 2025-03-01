# orsum 1.8.0

## Changes, bug fixes
- Replaced "np.unicode_", which was removed in the NumPy 2.0 release, with "np.str_". orsum still works with previous NumPy versions.
- Changed os.mkdir() to os.makedirs() in output folder creation, which allows nested directory creation.
- Fixed issues related to the handling of a single representative term.
- Added handling of duplicate terms in the input files (the first appearance is used to determine the rank).


# orsum 1.7.0

## Bug fix
- When using multiple enrichment results, if one input enrichment result is empty it is not added to the results. Before the fix, the file alias of the empty enrichment result persisted, causing shifting of aliases and disappearance of last alias.


# orsum 1.6.0

## Major Changes
- Fixed an issue in clustered heatmap creation function

## Minor Changes
- Code refactoring


# orsum 1.5.0

## Major Changes
- Added the function that additionally creates a clustered heatmap for the top n representative terms.

## Minor Changes
- Added maxTermSize parameter.
- Fixed an issue in file alias parsing.
- Fixed the issue of exiting when one of the input files is empty or all the terms in it are discarded due to minTermSize. Now this file is discarded and analysis continues with the remaining files.

# orsum 1.4.0
This is the version used in our publication <a href=https://pubmed.ncbi.nlm.nih.gov/35870894/>PMID: 35870894</a>
