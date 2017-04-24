# GSEA
Gene Set Enrichment Analysis (GSEA) is a method and a tool for functional interpretation of groups of
genes. These groups often share a common biological function, chromosomal location, or regulatory
schema. Gene sets have been defined over years of biomedical research. The goal of GSEA is to guide
the interpretation of gene expression data by linking important genes with known gene sets. You can
learn more about the method from the original scientific paper, available here:
http://www.pnas.org/content/102/43/15545.abstract . This article should also help you design the solution.

Input Data: 
1) leukemia.txt - This file includes gene expression profiles for two groups of leukemia patients (in columns). This example data includes a selected subset of genes from the originaldata set described here: https://github.com/ramhiser/datamicroarray/wiki/Golub(1999) . 
2) pathways.txt - This file is a list of gene sets that describe the metabolic pathways (pathways.txt).

Output:
This R scripts defines a function to calculate the gene set enrichment statistics as described in the above paper. It
returns the list of gene sets with their corresponding normalized enrichment scores and nominal pvalues.
Results are sorted by the normalized enrichment score (descending).



