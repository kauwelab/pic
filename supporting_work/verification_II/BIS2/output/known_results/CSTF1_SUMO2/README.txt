I======================== 
README file for summary/ 
======================== 

summary/ contains input and output files for BIS2.

alignment.fasta
---------------  

Multiple Sequence Alignment provided by the user, in FASTAformat.

tree.txt
--------  

Phylogenetic tree in Newick format, either computed with BioNJ (default), with 
PhyML, or provided by the user.

PF00000-CLUSTERFILE-5DD-dimD.html
---------------------------------  

Position of BIS2 clusters on the MSA, for dimension D. 
Only the first 30 sequences of the MSA are shown.
Amino acids in the MSA are coloured according to the physico-chemical class they 
belong to.
A histogram under the MSA representation reports the conservation level of a 
position (i.e., the frequency of the most abundant amino acid in a column).
Positions in clusters are labelled H (hit) or E (extension).

PF00000-CLUSTERFILE-5DD-dimD-table.html
---------------------------------------

Details on BIS2 clusters, for dimension D.
A table with six columns reports the following information:

Dim:                     dimension D
Cluster:                 cluster ID, 1-based
Sym:                     symmetrical score for the cluster, computed by CLAG
Env:                     environmental score for the cluster, computed by CLAG
Pvalue:                  p-value of the pattern, computed with a Fisher test
Hit patterns and blocks: 1-based position of hits in the MSA, with a table
                         representing the coevolution pattern, and blocks 
                         position (hits + extensions). If the pattern is fully
                         conserved (up to exceptions), the table is not shown.

PF00000-CLUSTERFILE-5DD-dimD-table.txt
--------------------------------------

Details on BIS2 clusters, for dimension D (tsv format).
A table with seven columns reports the following information:

list_of_positions:      1-based positions of hits in the MSA, separated by "/"
number_of_sequences:    the number of sequences for each distinguished word, 
                        separated by "/". 
aa_type:                coevolution pattern, namely, the distinguished words,
                        separated by "/". Together with "number_of_sequences", 
                        it contains the same information reported in the table 
                        of "Hit patterns and blocks" in 
                        PF00000-CLUSTERFILE-5DD-dimD-table.html
p_value:                p-value of the pattern, computed with a Fisher test
cluster_id:             cluster ID, 1-based
Ssymm:                  symmetrical score for the cluster, computed by CLAG
Senv:                   environmental score for the cluster, computed by CLAG


