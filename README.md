# orsum
orsum, which stands for over-representation summary, is a tool for the summary of enriched terms resulting from one or more enrichment analyses. orsum works on hierarchical annotations, e.g. GO and REACTOME, and finds representative terms to summarize the results by applying rules based on subterm-superterm relations.<br>

As input, orsum takes a <a href=https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29>GMT file</a> contaning the gene sets, and one or more files containing the enriched term IDs.<br>

The subterm-superterm relations among terms are inferred from the GMT file. The GMT file for the annotations can be found from multiple sources, for consistency, if it is available, it is better to download the exact GMT file that is used in enrichment analysis. If you already use <a href=https://biit.cs.ut.ee/gprofiler/gost>g:Profiler</a> for enrichment or GMT file was not available from your enrichment tool, getting from <a href=https://biit.cs.ut.ee/gprofiler/static/gprofiler_hsapiens.name.zip>this g:Profiler data source link</a> is an option. This zip file contains GMT for Reactome, GO:BP, GO:MF and GO:CC which you can use with orsum.<br><br>
The file containing the enriched term IDs (e.g. GO:0008150 or REAC:R-HSA-1640170) must contain one ID per line and must be sorted, with the most significant term at top.<br>

The user can specify which rules to apply, their order, and the maximum size for a term to be used as a representative term.<br>

The terms initially get their rank based on their order in the input file and each is a representative term (representing only itself). 
Throughout the orsum run, the representative terms get their rank from the best ranked term they represent, and the resulting terms are ordered by this rank.<br>

orsum produces multiple files as output:<br>
<ul>
	<li> An HTML file which presents the summarized results - the list of representative terms - with an option to click on each term to see the represented terms.
	<li> Two TSV files (-Summary.tsv, -Detailed.tsv) which contain the information in the HTML file in different formats, with some extras, like term size. -Summary.tsv file contains only representative terms while -Detailed.tsv additionally contains representative terms. -Summary.tsv file is very useful for getting an overview when comparing multiple enrichment results.
	<li> Three figures for<br>
		<ul>
			<li> the representative terms and how many terms they represent
			<li> the ranks of the representative terms (for each enrichment result if multiple enrichment results are input)
			<li> represented term size vs rank. <br>
		</ul>
</ul>
<br>

In order to use the tool you can either download the .py files from this repository and run orsum.py or you can download from bioconda.<br>
<code>
conda install -c bioconda orsum
</code><br>
<br>
<br>
usage: orsum.py [-h] --gmt GMT --files FILES [FILES ...]
                [--outputFolder OUTPUTFOLDER] [--rules [RULES [RULES ...]]]
                [--maxRepSize MAXREPSIZE] [--outputAll]
<br>
<ul>
<li>--gmt: Path for the GMT file. (required)
<li>--files: Paths for the enrichment result files. (required)
<li>--outputFolder: Path for the output result files. If it is not specified, results are written to the current directory. (optional, default=".")
<li>--rules: List of ordered numbers of the rules to apply while summarizing. By default, all the rules are run from 1 to 5. (optional)
<li>--maxRepSize: The maximum size of a representative term. Terms bigger than this will not be discarded but also will not be used to represent other terms. (optional, default=2000)
<li>--outputAll: When this option is used, a summary file is created after applying each rule, otherwise only final summary is created. (optional)
</ul>
<br>

Example command:<br>
<code>
orsum.py --gmt 'hsapiens.GO:BP.name.gmt' --files 'Enrichment-GOBP.csv' --outputFolder 'OutputGOBP'
</code><br>


Example command:<br>
<code>
orsum.py --gmt 'hsapiens.REAC.name.gmt' --files 'Enrichment-Method1-Reac.csv' 'Enrichment-Method2-Reac.csv' 'Enrichment-Method3-Reac.csv' --rules 1 2 5 --outputFolder 'OutputReac' --maxRepSize 1000
</code><br>

The rules are given below:
<ol>
	<li> Super &lt;- Sub (w worse rank) || Superterms represent their less significant subterms. This includes terms with exactly the same genes.
	<li> Sub &lt;- Super (w worse rank, subterm is 75% of superterm) || Subterms with more significance represent their superterms whose gene list is at least 75% composed of the subterm.
	<li> Sub &lt;- Super (w worse rank less representative power) || Subterms represent their superterms with less significance / lower rank and less representative power. Representative power is the number of terms they represent.
	<li> Common Super &lt;- (Sub Sub) (w less or equal representative power) || A superterm represents its two subterms that have less representative power.
	<li> Large Super &lt;- Large Sub || Superterms larger than maxRepSize represent their less significant subterms which are also larger than maxRepSize.
</ol>
<br>
<br>

This work is conducted in the Networks and Systems Biology for Diseases group of Anaïs Baudot (https://www.marseille-medical-genetics.org/a-baudot/). 
