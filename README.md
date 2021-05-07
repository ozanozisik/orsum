# orsum
orsum, which stands for over-representation summary, is a tool for the summary of enriched terms resulting from one or more enrichment analyses. orsum works on hierarchical annotations, e.g. GO and REACTOME, and finds representative terms to summarize the results by applying several rules like the ones below:</br>

<ul>
<li>Superterm covers less significant subterm</li>
<li>Sibling terms are represented by the parent</li>
<li>Subterm with more representative power covers less significant superterm</li>
</ul>

As input, orsum takes a GMT file contaning the gene sets and file(s) containing the enriched term IDs.
The hierarchy among terms is inferred from the GMT file..</br>
The file containing the enriched term IDs must contain one ID per line and must be sorted, with the most significant term at top.</br>


The user can specify which rules to apply, their order, and the maximum size for a term to be used as a representative term.</br>

The method produces two files after applying each rule: one summarized list containing only the representative terms and one detailed list containing representative terms and the represented terms. The representative terms get their rank from the best ranked term they represent, and the resulting terms are ordered by this rank.</br>

In order to use the tool you can either download the .py files from this repository and run orsum.py or you can download from bioconda.</br>
<code>
conda install -c bioconda orsum
</code></br>

Example command:</br>
<code>
orsum.py --gmt 'hsapiens.REAC.name.gmt' --hierarchy 'hierarchyDict-Reac.tsv' --files 'Enrichment-Method1-Reac.csv' 'Enrichment-Method2-Reac.csv' 'Enrichment-Method3-Reac.csv' --rules 1 2 3 4 8 9 10 --outputFolder 'DemoOutput' --maxRepSize 2000
</code></br>

Full list of rules are given below:
<ol>
<li>Terms from multiple results are matched</li>
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



