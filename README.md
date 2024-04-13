# PAR1_MRdrugtarget_Project

The contents of this repository relate to the manuscript 'Association between genetically proxied PAR-1 inhibition and kidney functions: a drug-target Mendelian Randomization study' by Haotian Tang, Venexia M Walker and Tom R Gaunt.

# Link to the manuscript

Will copy it here once we submit to arxiv

# Abstract

Background: Recently nephrotic syndrome (NS) patients have received anti-protease-activated receptor-1 (PAR-1) therapies in clinical trials. This study aimed to investigate the causal relationship of PAR-1 with NS and chronic kidney disease (CKD) using drug-target Mendelian randomization (MR).

Methods: We assessed the effect of genetically proxied PAR-1 on NS and CKD, using SNPs in and near the \textit{PAR1} gene (p-value $<$ 5e-8) from UK Biobank-Pharma Proteomics Project (UKB-PPP), eQTLGen and Genotype-Tissue Expression (GTEx) portal as instruments. Outcome data were obtained from a GWAS meta-analysis of FinnGen and UKB for NS and a CKDGen GWAS for CKD. PAR-1 has a known effect on platelet aggregation, so we used venous thromboembolism (VTE) from FinnGen as a positive control outcome.

Results: We present our findings as odds ratios (ORs) with 95\% confidence intervals (CIs) per unit change in plasma protein or gene expression level. We found genetically proxied PAR-1 reduced the risk of CKD with ORs of 0.85 [0.72 to 1.01] (UKB-PPP), 0.89 [0.81 to 0.98] (eQTLGen) and 0.80 [0.68 to 0.95] (GTEx). Genetically proxied PAR-1 increased the risk of VTE with ORs of 1.18 [ 1.00 to 1.40] (UKB-PPP), 1.06 [0.99 to 1.14] (eQTLGen) and 1.12 [0.94 to 1.33] (GTEx). The effect of PAR1 on NS was inconclusive due to a lack of power.

Interpretation: Our study provides support for genetically proxied PAR-1 activation causing an increased risk of VTE and a reduced risk of CKD. Evidence regarding NS was inconclusive.

# Questions
Please send any questions to Haotian Tang (haotian.tang@bristol.ac.uk).

# Script documentation

Each script within the script repo includes a documentation at the beginning that clarifies the purpose of the script. The numbers in each script name indicate the order in which the code should be run.

Version of R: 4.2.2

Main packages that need to be installed for running scripts: 

[TwoSampleMR 0.5.11](https://github.com/MRCIEU/TwoSampleMR)

[coloc 5.2.3](https://chr1swallace.github.io/coloc/)

[ieugwasr 0.1.8](https://mrcieu.github.io/ieugwasr)

[forestplot 3.1.3](https://github.com/gforge/forestplot)

