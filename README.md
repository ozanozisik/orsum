# orsum
orsum, which stands for over-representation summary, is a tool for the summary of enriched terms resulting from one or more enrichment analyses. orsum works on hierarchical annotations, e.g. GO and REACTOME, and finds representative terms to summarize the results by applying several rules like the ones below:</br>

<ul>
<li>Superterm covers less significant subterm</li>
<li>Sibling terms are represented by the parent</li>
<li>Subterm with more representative power covers less significant superterm</li>
</ul>

As input, orsum takes a <a href=https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29>GMT file</a> contaning the gene sets, and one or more files containing the enriched term IDs.</br>

The hierarchy among terms is inferred from the GMT file. The GMT file for the annotations might be found from multiple sources, for consistency, if it is available, it is better to download the exact GMT file that is used in enrichment analysis. If it is not convenient, another option is getting from <a href=https://biit.cs.ut.ee/gprofiler/static/gprofiler_hsapiens.name.zip>this g:Profiler data source link</a>. This zip file contains GMT for Reactome, GO:BP, GO:MF and GO:CC which you can use with orsum.</br></br>
The file containing the enriched term IDs must contain one ID per line and must be sorted, with the most significant term at top.</br>
</br>



The user can specify which rules to apply, their order, and the maximum size for a term to be used as a representative term.</br>

The method produces two files after applying each rule: one summarized list containing only the representative terms and one detailed list containing representative terms and the represented terms. The representative terms get their rank from the best ranked term they represent, and the resulting terms are ordered by this rank.</br>

In order to use the tool you can either download the .py files from this repository and run orsum.py or you can download from bioconda.</br>
<code>
conda install -c bioconda orsum
</code></br>
</br>
</br>
usage: orsum.py [-h] --gmt GMT --hierarchyFile HIERARCHYFILE [--createHF]
                --files FILES [FILES ...] --outputFolder OUTPUTFOLDER
                [--rules [RULES [RULES ...]]] [--maxRepSize MAXREPSIZE]
</br>
<ul>
<li>--gmt: Path for the GMT file. (required)
<li>--hierarchyFile: Path for the hierarchy file. It is created if the file does not exist. (required)
<li>--files: Paths for the enrichment result files. (required)
<li>--outputFolder: Path for the output result files. If it is not specified, results are written to the current directory. (optional, default=".")
<li>--createHF: Forces the creation of new hierarchy file, otherwise orsum tries to read, if absent creates it. (optional)
<li>--rules: List of ordered numbers of the rules to apply while summarizing. First rule is numbered 1. It should be run first. By default, all the rules are run from 1 to 10. (optiona)
<li>--maxRepSize: The maximum size of a representative term. Terms bigger than this will not be discarded but also will not be used to represent other terms. By default, maxRepSize = 5000 but we advise using a lower number, the default might change to 2000 soon. (optional, default=5000)
<li>--outputAll: When this option is used, a summary file is created after applying each rule, otherwise only final summary is created. (optional)
</ul>
</br>

Example command:</br>
<code>
orsum.py --gmt 'hsapiens.GO:BP.name.gmt' --hierarchy 'hierarchyDict-GOBP.tsv' --files 'Enrichment-GOBP.csv' --outputFolder 'DemoOutputGOBP'
</code></br>


Example command:</br>
<code>
orsum.py --gmt 'hsapiens.REAC.name.gmt' --hierarchy 'hierarchyDict-Reac.tsv' --files 'Enrichment-Method1-Reac.csv' 'Enrichment-Method2-Reac.csv' 'Enrichment-Method3-Reac.csv' --rules 1 2 3 4 8 9 10 --outputFolder 'DemoOutputReac' --maxRepSize 2000
</code></br>

Full list of rules are given below:
<ol>
<li>Different terms containing same genes are unified, terms from multiple results are matched</li>
<li>Superterm covers less significant subterm</li>
<li>Sibling terms whose parent is also in the list are represented by the parent</li>
<li>Cousin terms (they have a common grandparent) are represented by the grandparent term which is in the list</li>
<li>Aunt - nephew terms (one’s parent is other’s grandparent) are represented by the parent/grandparent term which is in the list</li>
<li>Subterm with more representative power covers less significant superterm. Representative power is currently the number of terms they represent</li>
<li>Subterm covers less significant superterm if rank difference is less than threshold (currently set to < 2)</li>
<li>Sibling terms are represented by the parent (even if the parent is not in the list)</li>
<li>Cousin terms (they have a common grandparent) are represented by the grandparent (even if the grandparent is not in the list)</li>
<li>Aunt - nephew terms (one’s parent is other’s grandparent) are represented by the parent/grandparent term (even if the parent/grandparent is not in the list)</li>
</ol>
</br>



