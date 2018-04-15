
# Microbial Analysis 16S rRNA for gene sequencing by QIIME
****
16S gene sequencing has developed over the years to understand the microbial phylogeny and taxonomy using various diversity analysis techniques. This repository will focus on the usage of Quantitative Insights into Microbial Ecology (QIIME) in investigating how the environmental conditions such as pH, weather temperature and mineral levels affect the microbiome diversity in the soil samples collected from various location in the Gordon Square Park, London.


The data used is the result of Illumina sequencing of the genetic marker, 16S rRNA V4 region, that is found in the soil sample. For speedy analysis, high performance computer that is based in Amsterdam known as Cartesius is used for analysis by installing the QIIME environment into the server.

 
**Note**: The steps that will be listed may not be the best methods as everyone has their own way of defining what is the ‘best’ steps should be included and not included.

 
**Breakdown of the steps for this diversity analyses:**
1. Validate mapping file
2. Join paired ends
3. Demultiplexing and quality filtering
4. Picking operational taxonomical units (OTUs) using open and closed reference
5. Core diversity analysis
6. Comparing categories using Adonis and ANOSIM 
7. Compare distance matrices
*****

## Validate mapping file

Since mapping file is created by users, it is prones to errors. The mapping file contains necessary informations such as sampleID, barcode sequence, linker primer sequence, sample description, sample coordinates, sample depths, weather temperature, sample site descriptions and soil contents in terms of pH, Nitrogen (N), Phosphorus (P) and Potassium (K). In order to check for errors in the mapping file, the script [validate_mapping_file.py] (http://qiime.org/scripts/validate_mapping_file.html) would produce a `map.tsv.log` output which describes the error.   
```
# validating map file
echo "check id map"
time validate_mapping_file.py -m map.tsv -o mapping_output	
```
Output containing:
```
overlib.js   map.tsv_corrected.txt   map.tsv.log  map.tsv.html

```
*****
## Joining paired ends

Before the data is demultiplexed, we used the [join_paired_end.py](http://qiime.org/scripts/join_paired_ends.html) to join the forward and reverse fastq files Illumina reads using fastq-join method by [Erik Aronesty](https://expressionanalysis.github.io/ea-utils/). 
```
#joining paired ends
echo "join paired ends"
time join_paired_ends.py -m fastq-join -f Read1.fastq -r Read2.fastq -b Index.fastq -o Read12/fastq-join_joined
```
The output directory will contain the join paired ends files:
```
fastqjoin.un1.fastq   fastqjoin.un2.fastq   fastqjoin.join.fastq   fastqjoin.join_barcodes
```
*******
## Demultilpexing and quality filtering

Demultiplexing assigns the multiplexed read from Illumina reads to the sample based of their barcode sequences and simultaneously perform quality filtering `-q` to remove ambiguity. The [‘split_libraries_fastq.py](http://qiime.org/scripts/split_libraries_fastq.html) is used:
```
# splitting libraries
echo "splitting libraries"
time split_libraries_fastq.py \
--barcode_type 12 \
-m mapping_output/map.tsv -i Read12/fastq-join_joined/fastqjoin.join.fastq \
-b Read12/fastq-join_joined/fastqjoin.join_barcodes.fastq \
-o dem1 \
-q 3 \

--rev_comp_barcode \
--rev_comp_mapping_barcodes
```
which gives the output of:
```
seqs.fna  split_library_log.txt   histograms.txt

``` 
The script created three files at the output directory known as split_library_log.txt, histograms.txt and seqs.fna. Both tab-delimited files show summary of splitting and the number of reads at regular intervals respectively. Furthermore, the fasta formatted file contains each sequence renamed according to the sample where it came from. 
*****
## Picking OTUs

The resulting demultiplexed read consists of unknown sequences. These sequences are clustered together with other related sequences based on their DNA sequence similarity of 16S rRNA gene and then quantified. The clusters of similar sequencing reads are defined as operational taxonomic units (OTUs).  This pragmatic computational approach caters for sequence analysis and taxonomic profiling using the approach of assigning known sequences to OTUs clusters of unknown sequences through sequence similarity by referring to existing rRNA databases, which in this experiment we used Greengenes (GG) only. There are 3 different strategies to pick OTUs which are de novo, closed and open reference. [Open referencing] (http://qiime.org/scripts/pick_open_reference_otus.html) is deemed as the best approach as the script involves a total of 6 steps to pick OTUs which suggest the thouroughness of this method. In this data analysis, I am going to use both open and [closed reference] (http://qiime.org/scripts/pick_closed_reference_otus.html) for comparison.

Closed reference script:
```
# picking OTUs
echo "Picking OTUs with closed reference"
time pick_closed_reference_otus.py \
-i dem1/seqs.fna \
-o otus \
```
producing output:
```
uclust_ref_picked_otus   otu_table.biom     log_20180406235444.txt     97_otus.tree
```
Open reference script:
```
# picking OTUs
echo "Picking OTUs with open reference"
time pick_open_reference_otus.py \
-i dem1/seqs.fna \
-o otusopen \
-a \
-O 16
```
producing output:
```
step1_otus
step2_otus
step3_otus
final_otu_map.txt
step4_otus
final_otu_map_mc2.txt
new_refseqs.fna
rep_set.fna
otu_table_mc2.biom
uclust_assigned_taxonomy
otu_table_mc2_w_tax.biom 
pynast_aligned_seqs
rep_set.tre
otu_table_mc2_w_tax_no_pynast_failures.biom
log_20180407234530.txt
index.html
```
The OTU table is summarized using `biom summarize-table -i otusopen/otu_table_mc2.biom` to determine the sampling depth parameter to include in the core diversity analyses.

An example of summarized BIOM table for open reference:
```
Num samples: 30
Num observations: 68847
Total count: 5343652
Table density (fraction of non-zero values): 0.112

Counts/sample summary:
 Min: 765.0
 Max: 1827022.0
 Median: 123337.500
 Mean: 178121.733
 Std. dev.: 310684.360
 Sample Metadata Categories: None provided
 Observation Metadata Categories: taxonomy

Counts/sample detail:
515rcbc20: 765.0
515rcbc36: 39701.0
515rcbc8: 47167.0
515rcbc13: 66726.0
515rcbc37: 69375.0
515rcbc34: 71216.0
515rcbc9: 72835.0
515rcbc27: 77928.0
515rcbc23: 82008.0
515rcbc35: 100004.0
515rcbc32: 102661.0
515rcbc11: 108531.0
515rcbc33: 115639.0
515rcbc17: 120272.0
515rcbc15: 120409.0
515rcbc28: 126266.0
515rcbc18: 128548.0
515rcbc21: 135789.0
515rcbc14: 137995.0
515rcbc25: 142976.0
515rcbc31: 143392.0
515rcbc30: 151779.0
515rcbc19: 163258.0
515rcbc26: 176923.0
515rcbc29: 185631.0
515rcbc16: 193824.0
515rcbc10: 203236.0
515rcbc24: 206755.0
515rcbc22: 225021.0
515rcbc12: 1827022.0
```

************
## Core diversity analysis

By running the [core_diversity_analyses.py](http://qiime.org/scripts/core_diversity_analyses.html), the script would calculate several diversity measurements within the sample, known as the alpha diversity metric and beta diversity metric which is between the samples. For better understanding of the metrices, these numbers are visualized using plots to show the measurement relativity to the sample diversity in the community. This scripts is a combination of other scripts such as: [alpha_rarefaction.py](http://qiime.org/scripts/alpha_rarefaction.html), [beta_diversity_through_plots.py](http://qiime.org/scripts/beta_diversity_through_plots.html) and [summarize_taxa_through_plots.py](http://qiime.org/scripts/summarize_taxa_through_plots.html).
```
# core diversity analyses
echo "Core diversity analyses"
time core_diversity_analyses.py \
--recover_from_failure \
-i otusopen1/otu_table_mc2_w_tax_no_pynast_failures.biom \
-m map.tsv \
-t otusopen1/rep_set.tre \
-e 764 \
-o cda-open \
-a \
-O 16
```
This command gives output:
```
arare_max765   
bdiv_even765   
biom_table_summary.txt   
index.html   
log_20180406211900.txt   
table_even616.biom.gz   
table_mc616.biom.gz   
taxa_plots
```

Alpha diversity (arare_max765 in this case) describes the numbers and diversity of species within each sample and there are different metrics that have been devised to measure alpha diversity with emphasis on different aspects of community structure. This metric can be measured using 3 strategies which are observed OTUs, PD whole tree and Chao1 estimator. These diversity metrices can be further visualised using rarefaction plots (number of species against sequencing depth). The seqeuncing depth parameter is determined by the lowest OTU counts from `biom summarize-table` command.


Furthermore, the diversity between samples is measured using the beta diversity metrics (in this example bdiv_even765) which analyzed based on distance scores using UniFrac distances (unique fraction). UniFrac is a method that measures the phylogenetic distance between sets of taxa in a phylogenetic tree  where it measures the distance between communities based on the their lineages from either one sample environment or the other but not both. UniFrac distance can be further quantified by using the relative abundance of OTU counts to reveal the community differences that are due to these abundances. This measured is identified as weighted UniFrac distance. Besides, UniFrac distance can also be qualitatively measured by using the presence or absence of data to emphasis rare OTUs. This is referred to as unweighted UniFrac distance. Applying qualitative and quantitative measure to the same sample can provide insights on the nature of community differences. Unweighted and weighted UniFrac distance metrices can be visualised using the principal coordinate analysis (PCoA) plot. 
*********
## Comparing categories

The [compare_categories.py](http://qiime.org/scripts/compare_categories.html) allows for the analysis of the strength and statistical significance of sample groupings using a distance matrix as the primary input. We will use Adonis, a nonparametric statistical method that takes QIIME distance matrix file such as UniFrac distance matrix, a mapping file and a category in the mapping file and a category in the mapping file to determine sample grouping from.

Adonis creates a set by first identifying the relevant representative OTU sample of the data and then calculating the squared deviations from these points. After that, significance tests are performed using F-tests based on sequential sums of squares from permutations of the raw data. ANOVA uses the F-test to determine whether the variability between group means is larger than the variability of the observations within the groups. If that ratio is sufficiently large, you can conclude that not all the means are equal. Because the F-distribution assumes that the null hypothesis is true, we can place the F-value from our study in the F-distribution to determine how consistent our results are with the null hypothesis and to calculate probabilities. The probability that we want to calculate is the probability of observing an F-statistic that is at least as high as the value that our study obtained. That probability allows us to determine how common or rare our F-value is under the assumption that the null hypothesis is true. If the probability is low enough, we can conclude that our data is inconsistent with the null hypothesis. The evidence in the sample data is strong enough to reject the null hypothesis for the entire population.

To compare category (sample pH) with the distance metrics:
```
# comparing categoris for statistical significant
echo "compare categories"
time compare_categories.py \
--method adonis \
-i cda-open1/bdiv_even765/unweighted_unifrac_dm.txt \
-m map.tsv \
-c SamplePh \
-o adonis_open_out \
-n 999
```
will give an output:
```
Call:
adonis(formula = as.dist(qiime.data$distmat) ~ qiime.data$map[[opts$category]],      permutations = opts$num_permutations)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                                Df SumsOfSqs MeanSqs F.Model    R2 Pr(>F)
qiime.data$map[[opts$category]]  1    0.2310 0.23096  1.2275 0.042  0.044 *
Residuals                       28    5.2685 0.18816         0.958
Total                           29    5.4994                 1.000
```
In this example, R^2 = 0.042 = 4% of the sums of squares can be explained by the SamplePh. The null hypothesis is there is no significant difference between samples in Category, any observed difference being due to sampling or experimental error. The p-value = 0.044, if the statistical significant is defined as p<0.05, then null hypothesis is rejected.

Analysis of similarity (ANOSIM), another nonparametric statistical significance test, is to tests whether two or more groups of samples are significantly different based on a categorical variable found in the metadata mapping file. You can specify a category in the metadata mapping file to separate samples into groups and then test whether there are significant differences between those groups. 
An example of the script, comparing the same category as Adonis, Sample pH:
```
# comparing categoris fro statistical significant
echo "compare categories"
time compare_categories.py \
--method anosim \
-i cda-open1/bdiv_even765/unweighted_unifrac_dm.txt \
-m map.tsv \
-c SamplePh \
-o anosim_open_out \
-n 999
```
This will give output as:
```
method name     ANOSIM
test statistic name     R
sample size     	30
number of groups        9
test statistic  -0.029120700843608099
p-value 	0.58599999999999997
number of permutations  999
```
An R value near +1 means there is dissimilarity between groups for mapping file category whilst R value near 0 indicates no significant dissimilarity between the groups.In this example, for sample pH, R statistic of -0.029 means there is no statistical difference in microbial composition between samples and p-value of 0.585, showed no statistical significant given we define the statistical difference is (p<0.05). This could mean for ANOSIM method, the sample pH does not affect the microbial composition between samples.  

******
##Other diversity analysis methods that could be considered 
### Comparing distance matrices through Mantel test
The input for this script is a mapping file and the name of a column, it has to be numeric, from which a distance matrix will be created. The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison. Therefore, the sample N, K and P categories is not able to create distance matrices. The 5 categories (Longitude, Latitude, Depth, pH and weather temperature) is run with [distance_matrix_from_mapping.py](http://qiime.org/scripts/distance_matrix_from_mapping.html) to create .dm.txt files.
By using  [compare_distance_matrices.py](http://qiime.org/tutorials/distance_matrix_comparison.html), all of the currently available comparison techniques are based on the Mantel test, which is a non-parametric statistical method that computes the correlation between two distance matrices. In addition to this statistical method, QIIME also provides the partial Mantel test and Mantel correlogram. One common application of distance matrix comparison techniques is to determine if correlation exists between a community distance matrix (e.g. UniFrac distance matrix) and a second matrix derived from an environmental parameter that is numeric/continuous in nature (e.g. pH, Latitute, Longitude, depth and weather temperature). Mantel test to test the correlation between 2 matrices as shown below:
```
# distance matrix from mapping
echo " distance matrix from mapping"
time compare_distance_matrices.py \
--method=mantel \
-i unweighted_unifrac_dm.txt,PH_dm.txt, Lat_dm.txt,Long_dm.txt,Weather_dm.txt,Depth_dm.txt \
-o comparedm/all_open \
-n 999
```
This will give an output of Mantel statistic:
```
DM1     DM2     Number of entries       Mantel r statistic      p-value Number of permutations  Tail type
unweighted_unifrac_dm.txt       weighted_unifrac_dm.txt 30      0.82452 0.001   999     two-sided
unweighted_unifrac_dm.txt       PH_dm.txt       30      0.02905 0.407   999     two-sided
unweighted_unifrac_dm.txt       Lat_dm.txt      30      0.12229 0.370   999     two-sided
unweighted_unifrac_dm.txt       Long_dm.txt     30      0.04409 0.734   999     two-sided
unweighted_unifrac_dm.txt       Depth_dm.txt    30      0.10678 0.429   999     two-sided
unweighted_unifrac_dm.txt       Weather_dm.txt  30      -0.08082        0.435   999     two-sided
weighted_unifrac_dm.txt PH_dm.txt       30      0.01666 0.577   999     two-sided
weighted_unifrac_dm.txt Lat_dm.txt      30      0.10321 0.413   999     two-sided
weighted_unifrac_dm.txt Long_dm.txt     30      -0.01557        0.897   999     two-sided
weighted_unifrac_dm.txt Depth_dm.txt    30      0.11316 0.382   999     two-sided
weighted_unifrac_dm.txt Weather_dm.txt  30      -0.11522        0.254   999     two-sided
PH_dm.txt       Lat_dm.txt      30      -0.00772        0.756   999     two-sided
PH_dm.txt       Long_dm.txt     30      0.01544 0.650   999     two-sided
PH_dm.txt       Depth_dm.txt    30      -0.02707        0.401   999     two-sided
PH_dm.txt       Weather_dm.txt  30      0.01413 0.678   999     two-sided
Lat_dm.txt      Long_dm.txt     30      0.92115 0.002   999     two-sided
Lat_dm.txt      Depth_dm.txt    30      0.17292 0.161   999     two-sided
Lat_dm.txt      Weather_dm.txt  30      0.28072 0.054   999     two-sided
Long_dm.txt     Depth_dm.txt    30      0.00705 0.980   999     two-sided
Long_dm.txt     Weather_dm.txt  30      0.36869 0.028   999     two-sided
Depth_dm.txt    Weather_dm.txt  30      -0.08915        0.407   999     two-sided```
The Mantel r statstic determines the correlation between the two distance matrices. An r value +1 indicates strong positive correlation while -1 means strong negative correlation. 
********
### Mapping phylogenetic tree using Cytoscape
Using [make bipartite network] (http://qiime.org/scripts/make_bipartite_network.html), we could visualise the phylogenetic tree of the microbial composition in the soil other than summarizing the taxa through plots in core diversity analysis step.





