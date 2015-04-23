# Ahya-White-Syndromes

##DESeq2## for gene expression analysis requires counts data file "25july_diseasecounts.tab"
This file is a table of counts per gene: rows are genes, columns are samples
Sample name abbreviations are as follows: AL = ahead of the lesion, D = disease, H = healthy
Numbers following sample names are the genotypes of the individuals: Genotypes 1-8 were sampled twice (disease and ahead of the lesion tissues were sampled), genotypes 9-16 were sampled once (healthy colonies)

WGCNA  requires variance stabilized gene expression data for the most significantly differentially expressed genes (unadjusted pvalue < 0.1) "genes4WGCNA.csv"
Rows are genes, values are variance stabilized data. Columns are samples with the same abbreviations as stated above.

WGCNA also requires a table of "traits", in our case simply healthy state and individual genotype ID: "diseasetraits_reps.csv"
Columns are traits, rows are samples. 0 and 1 denote the presence and absence of these "traits", respectively.

PCoA script requires variance stabilized gene expression for all genes "VSDandPVALS_disease.csv". 
