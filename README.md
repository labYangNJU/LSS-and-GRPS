# LSS-and-GRPS
A novel approach to comprehensively explore genetic contribution and nominate reliable causal genes for complex diseases

## **About**

We introduce a locus-specific stratification (LSS) and gene regulatory prioritization score (GRPS) approach to comprehensively explore genetic contribution and nominate reliable causal genes for complex diseases. LSS efficiently extracts candidate causal variants from multi-signals based on the unique distribution of locus-specific background and overcomes the challenge of multi-signals in complex diseases. GRPS uncovered candidate causal genes with high credibility based on cumulative-regulation from multi-signals by integrating candidate causal variants from LSS with the regulatory network.

## **Getting Started**

This program requires R 4.0, several R packages. Versions the R and R packages has been tested on R (version = 4.3.2), data.table (version = 1.16.4), GenomicRanges (version = 1.58.0), ieugwasr (version = 1.0.3)

## **Demo**

We have provided example input data in the 'input' folder. You can run the example file using the code provided in the 'code' folder. The outputs generated from running above are provided in the 'output' folder.

***Input***\
locus_example_gwas_information.txt (GWAS summary statistics)\
kidney_scATACseq_peaks.txt (Functional annotation)\
kidney_regulation.rds (Regulatory network)

***output***\
locus_example_lss_profile (Quantile-quantile plot for each example locus)\
locus_example_lss_candidate_variants.csv (Identification of candidate causal variants)\
locus_example_lss_functional_variants.csv (Identification of functional SNPs based on functional annotation)\
locus_example_lss_variants_regulated_genes.csv (Identification of candidate causal genes based on regulatory network)\
locus_example_grps_result.txt（GRPS results）

Expected runtime for the demo on a normal desktop computer is approximately 10 to 15 minutes

## **Contact**

Jing Zhang （DG20350040@smail.nju.edu.cn）
