# Script for doing coloc -------------------------------------------------------

# nephQTL2 is 37 build ---------------------------------------------------------
# liftove for 38 build ---------------------------------------------------------

f2r_Glomneph2 = fread(paste0(rdsf_personal,"data/format_data/f2r_Glomneph2_GWAS_tidy_outcome.csv"))
f2r_Tubeneph2 = fread(paste0(rdsf_personal,"data/format_data/f2r_Tubeneph2_GWAS_tidy_outcome.csv"))

f2r_Glomneph2_marker = f2r_Glomneph2[,c("SNP","chr.outcome","pos.outcome")]
colnames(f2r_Glomneph2_marker) = c("rsid","chrom","pos")
f2r_Glomneph2_marker_hg38 = perform_liftover(f2r_Glomneph2_marker)
f2r_Glomneph2_hg3738 <- merge(f2r_Glomneph2, f2r_Glomneph2_marker_hg38, by.x = "SNP", by.y = "rsid_b37")

f2r_Tubeneph2_marker = f2r_Tubeneph2[,c("SNP","chr.outcome","pos.outcome")]
colnames(f2r_Tubeneph2_marker) = c("rsid","chrom","pos")
f2r_Tubeneph2_marker_hg38 = perform_liftover(f2r_Tubeneph2_marker)
f2r_Tubeneph2_hg3738 <- merge(f2r_Tubeneph2, f2r_Tubeneph2_marker_hg38, by.x = "SNP", by.y = "rsid_b37")

# ukb is 38 build
# eqtlgen is 37 build
# gtex is 38 build

# coloc between ukb-ppp and nephqtl2 -------------------------------------------
# leading SNP: rs168753
# chr pos 5 76732299
76732299 + 300000
76732299 - 300000

f2r_ukb = fread(paste0(rdsf_personal,"data/format_data/f2r_ukb_GWAS_tidy_outcome.csv")) %>% data.frame()
f2r_ukb_coloc = subset(f2r_ukb, chr.outcome ==5 & pos.outcome >= 76432299 & pos.outcome <= 77032299)
f2r_ukb_coloc <- convert_outcome_to_exposure_local(f2r_ukb_coloc,ignore_samplesize = T)
f2r_ukb_coloc$exposure = "F2R (UKB-PPP)"

f2r_ukb_Glomneph2_merge = merge(f2r_ukb_coloc, f2r_Glomneph2_hg3738, by = "SNP")
# check if position is right
# summary(f2r_ukb_Glomneph2_merge$pos.exposure == f2r_ukb_Glomneph2_merge$pos_b38)
f2r_ukb_Tubeneph2_merge = merge(f2r_ukb_coloc, f2r_Tubeneph2_hg3738, by = "SNP")
# summary(f2r_ukb_Tubeneph2_merge$pos.exposure == f2r_ukb_Tubeneph2_merge$pos_b38)

df = f2r_ukb_Glomneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
ukb_Glomneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 240,df$SNP) 
ukb_Glomneph2_result
# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2        H3      H4
# 0.5795156 0.1841 0.1841 0.0338744 0.01841
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1.841000e+03 1.528251e-45 7.646278e-01 3.673479e-46 1.837430e-01 5.162921e-02 

df = f2r_ukb_Tubeneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
ukb_Tubeneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 311,df$SNP) 
ukb_Tubeneph2_result

# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2         H3      H4
# 0.5866647 0.1812 0.1812 0.03281532 0.01812
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1.812000e+03 1.087527e-45 5.441206e-01 8.415268e-46 4.210049e-01 3.487449e-02 

# coloc between eqtlgen and nephqtl2 -------------------------------------------
# leading SNP: rs250735
# chr pos 5 76732299
75990144+300000
75990144-300000
# f2r_eqtlg_dat_format is from script 05-prepare_instruments_ukb_eqtlgen

f2r_eqtlg_coloc = subset(f2r_eqtlg_dat_format, pos.outcome >= 75690144 & pos.outcome <= 76290144)
colnames(f2r_eqtlg_coloc)[grep(".outcome", colnames(f2r_eqtlg_coloc))] <- gsub(".outcome", ".exposure", colnames(f2r_eqtlg_coloc)[grep(".outcome", colnames(f2r_eqtlg_coloc))])
f2r_eqtlg_coloc$exposure = "F2R (eQTLGen)"

f2r_eqtlg_Glomneph2_merge = merge(f2r_eqtlg_coloc, f2r_Glomneph2_hg3738, by = "SNP")
# check if position is right
# summary(f2r_eqtlg_Glomneph2_merge$pos.exposure == f2r_eqtlg_Glomneph2_merge$pos_b37)
f2r_eqtlg_Tubeneph2_merge = merge(f2r_eqtlg_coloc, f2r_Tubeneph2_hg3738, by = "SNP")
# summary(f2r_eqtlg_Tubeneph2_merge$pos.exposure == f2r_eqtlg_Tubeneph2_merge$pos_b37)

df = f2r_eqtlg_Glomneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
eqtlg_Glomneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 240,df$SNP) 
eqtlg_Glomneph2_result
# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2         H3      H4
# 0.5802559 0.1838 0.1838 0.03376406 0.01838
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1838.0000000    0.0000000    0.6570926    0.0000000    0.1561137    0.1867937 


df = f2r_eqtlg_Tubeneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
eqtlg_Tubeneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 311,df$SNP) 
eqtlg_Tubeneph2_result
# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2         H3      H4
# 0.5866647 0.1812 0.1812 0.03281532 0.01812
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1.812000e+03 0.000000e+00 6.136843e-02 0.000000e+00 4.766816e-02 8.909634e-01 



# coloc between gtex and nephqtl2 ----------------------------------------------
# leading SNP:rs250753
# chr pos: 5 76741726
76741726 + 300000
76741726 - 300000

f2r_gtexbld_gwas_mformat = fread(paste0(rdsf_personal,"data/par1/f2r_gtexbld_mformat.tsv.gz"))

f2r_gtexbld_format = format_data(data.frame(f2r_gtexbld_gwas_mformat),
                                           type = "exposure",
                                           snp_col = "SNP",
                                           beta_col = "BETA",
                                           se_col = "SE",
                                           eaf_col = "FRQ",
                                           effect_allele_col = "A2",
                                           other_allele_col = "A1",
                                           pval_col = "P",
                                           chr_col = "CHR",
                                           pos_col = "BP") %>% mutate(exposure = "F2R (GTEx)")

f2r_gtexbld_coloc = subset(f2r_gtexbld_format, pos.exposure >= 76441726 & pos.exposure <= 77041726)


f2r_gtexbld_Glomneph2_merge = merge(f2r_gtexbld_coloc, f2r_Glomneph2_hg3738, by = "SNP")
# check if position is right
# summary(f2r_gtexbld_Glomneph2_merge$pos.exposure == f2r_gtexbld_Glomneph2_merge$pos_b38)
f2r_gtexbld_Tubeneph2_merge = merge(f2r_gtexbld_coloc, f2r_Tubeneph2_hg3738, by = "SNP")
# summary(f2r_gtexbld_Tubeneph2_merge$pos.exposure == f2r_gtexbld_Tubeneph2_merge$pos_b38)

df = f2r_gtexbld_Glomneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
gtexbld_Glomneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 240,df$SNP) 
gtexbld_Glomneph2_result
# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2         H3      H4
# 0.5629344 0.1908 0.1908 0.03638556 0.01908
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1.908000e+03 4.116099e-04 7.255307e-01 1.019009e-04 1.795229e-01 9.443289e-02 

df = f2r_gtexbld_Tubeneph2_merge[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
df <- df[complete.cases(df), ]
df <- maf_format(df)
gtexbld_Tubeneph2_result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, 33784, 311,df$SNP) 
gtexbld_Tubeneph2_result

# SNP Priors
# p1    p2   p12 
# 1e-04 1e-04 1e-05 
# 
# Hypothesis Priors
# H0     H1     H2         H3      H4
# 0.5693796 0.1882 0.1882 0.03540042 0.01882
# 
# Posterior
# nsnps           H0           H1           H2           H3           H4 
# 1.882000e+03 1.676066e-05 2.954343e-02 9.589828e-06 1.594918e-02 9.544810e-01 

