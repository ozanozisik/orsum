# orsum
orsum, which stands for "over-representation summary", is a tool for the summary of enriched terms resulting from one or more enrichment analyses. orsum works on hierarchical annotations, e.g. GO and REACTOME, and finds representative terms among the original list to summarize the results.<br>

As input, orsum requires a <a href=https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29>GMT file</a> contaning the gene sets, and one or more files containing the enriched term IDs.<br>

The subterm-superterm relations among terms are inferred from the GMT file. The GMT file for the annotations can be found from multiple sources, for consistency, if it is available, it is better to download the exact GMT file that is used in enrichment analysis. If you already use <a href=https://biit.cs.ut.ee/gprofiler/gost>g:Profiler</a> for enrichment or GMT file was not available from your enrichment tool, getting from <a href=https://biit.cs.ut.ee/gprofiler/static/gprofiler_hsapiens.name.zip>this g:Profiler data source link</a> is an option. This zip file contains GMT for Reactome, GO:BP, GO:MF and GO:CC which you can use with orsum.<br><br>
The input files containing the enriched term IDs (e.g. GO:0008150 or REAC:R-HSA-1640170) must contain one ID per line and must be sorted, with the most significant term at top. There should not be a header.<br>

The terms initially get their rank based on their order in the input file and each is a representative term (representing only itself).
Throughout the orsum run, the representative terms get their rank from the best ranked term they represent, and the resulting terms are ordered by this rank.<br>

orsum produces multiple files as output:<br>
<ul>
	<li> An HTML file which presents the summarized results - the list of representative terms - with an option to click on each term to see the represented terms.
	<li> Two TSV files (-Summary.tsv, -Detailed.tsv) which contain the information in the HTML file in different formats, with some extras, like term size. -Summary.tsv file contains only representative terms while -Detailed.tsv additionally contains representative terms. -Summary.tsv file is very useful for getting an overview when comparing multiple enrichment results.
	<li> A TSV file that consists of two columns, mapping representative term IDs to representing term IDs. This file is for programmatic access.
	<li> Three figures for<br>
		<ul>
			<li> the representative terms and how many terms they represent
			<li> the ranks of the representative terms (for each enrichment result if multiple enrichment results are input)
			<li> represented term size vs rank of representative terms <br>
		</ul>
</ul>
<br>

In order to use the tool you can either download the .py files from this repository and run orsum.py or you can download from bioconda.<br>
<code>
conda install -c bioconda orsum
</code><br>
<br>
<br>
Usage:
<br>
<code>
orsum.py [-h] [-v] --gmt GMT --files FILES [FILES ...]
                [--fileAliases FILEALIASES [FILEALIASES ...]]
                [--outputFolder OUTPUTFOLDER] [--maxRepSize MAXREPSIZE]
                [--minTermSize MINTERMSIZE]
                [--numberOfTermsToPlot NUMBEROFTERMSTOPLOT]
</code>
<br>
<ul>
<li>--gmt: Path of the GMT file. (required)
<li>--files: Paths of the enrichment result files. (required)
<li>--fileAliases: Aliases for input enrichment result files to be used in orsum results. (optional, default file names)
<li>--outputFolder: Path for the output result files. If it is not specified, results are written to the current directory. (optional, default=".")
<li>--maxRepSize: The maximum size of a representative term. Terms bigger than this size will not be discarded but also will not be able to represent other terms. (optional, default=1000000, bigger than any annotation term, which means that it has no effect.)
<li>--minTermSize: The minimum size of the terms to be processed. Smaller terms will be discarded. (optional, default=10)
<li>--numberOfTermsToPlot: The number of representative terms to be presented in barplot and heatmap. (optional, default=50)
</ul>
<br>

Example command:<br>
<code>
orsum.py --gmt 'hsapiens.GO:BP.name.gmt' --files 'Enrichment-GOBP.txt' --outputFolder 'OutputGOBP'
</code><br>


Example command:<br>
<code>
orsum.py --gmt 'hsapiens.REAC.name.gmt' --files 'Enrichment-Method1-Reac.txt' 'Enrichment-Method2-Reac.txt' 'Enrichment-Method3-Reac.txt' --fileAliases 'Method 1' 'Method 2' 'Method 3' --outputFolder 'OutputReac' --maxRepSize 2000 --minTermSize 20 --numberOfTermsToPlot 20
</code><br>

orsum first assigns each term a rank based on its order in the corresponding input file.<br>
If there are multiple input enrichment results, there are multiple instances of the term initially, each getting its rank based on the corresponding input. These terms are unified and the best rank is used as the representative term rank.<br>
Next, orsum starts filtering terms by the following principle: Superterms will represent their less significant (worse ranked) subterms. This includes equal terms, i.e. the terms that annotate exactly the same set of genes. This is performed by iterating over terms starting with the highest ranking term, and checking each worse ranked term to see whether it is a subterm of the considered term. When the condition is satisfied, the subterm is discarded and represented by the considered term (representative term).
<br>
<br>
<br>


This work is conducted in the Networks and Systems Biology for Diseases group of Anaïs Baudot (https://www.marseille-medical-genetics.org/a-baudot/).
