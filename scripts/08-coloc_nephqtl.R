# Script for doing coloc where Trait 2 is from neph2 ---------------------------
# window size == ±300kb --------------------------------------------------------
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

df <- maf_format(f2r_ukb_Glomneph2_merge)
ukb_Glomneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome, 
                                             se1 = df$se.exposure, se2 = df$se.outcome, 
                                             MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                             N1 = 33784, N2 = 240, SNP = df$SNP) 
ukb_Glomneph2_result
# nsnps    PP.H0.abf PP.H1.abf    PP.H2.abf PP.H3.abf  PP.H4.abf
# 1  1841 1.528251e-45 0.7646278 3.673479e-46  0.183743 0.05162921

df <- maf_format(f2r_ukb_Tubeneph2_merge)
ukb_Tubeneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome, 
                                             se1 = df$se.exposure, se2 = df$se.outcome, 
                                             MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                             N1 = 33784, N2 = 311, SNP = df$SNP) 
ukb_Tubeneph2_result
# nsnps    PP.H0.abf PP.H1.abf    PP.H2.abf PP.H3.abf  PP.H4.abf
# 1  1812 1.087527e-45 0.5441206 8.415268e-46 0.4210049 0.03487449

# coloc between eqtlgen and nephqtl2 -------------------------------------------
# leading SNP: rs250733
# chr pos 5 76732299
75990144+300000
75990144-300000
# f2r_eqtlg_dat_format is from script 05-prepare_instruments_ukb_eqtlgen

f2r_eqtlg_coloc = subset(f2r_eqtlg_dat_format, pos.outcome >= 75690144 & pos.outcome <= 76290144)
colnames(f2r_eqtlg_coloc)[grep(".outcome", colnames(f2r_eqtlg_coloc))] <- 
  gsub(".outcome", ".exposure", colnames(f2r_eqtlg_coloc)[grep(".outcome", colnames(f2r_eqtlg_coloc))])
f2r_eqtlg_coloc$exposure = "F2R (eQTLGen)"

f2r_eqtlg_Glomneph2_merge = merge(f2r_eqtlg_coloc, f2r_Glomneph2_hg3738, by = "SNP")
# check if position is right
# summary(f2r_eqtlg_Glomneph2_merge$pos.exposure == f2r_eqtlg_Glomneph2_merge$pos_b37)
f2r_eqtlg_Tubeneph2_merge = merge(f2r_eqtlg_coloc, f2r_Tubeneph2_hg3738, by = "SNP")
# summary(f2r_eqtlg_Tubeneph2_merge$pos.exposure == f2r_eqtlg_Tubeneph2_merge$pos_b37)

df <- maf_format(f2r_eqtlg_Glomneph2_merge)
eqtlg_Glomneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome, 
                                               se1 = df$se.exposure, se2 = df$se.outcome, 
                                               MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                               N1 = 31684, N2 = 240, SNP = df$SNP) 
eqtlg_Glomneph2_result
# nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf
# 1  1838         0 0.6570929         0 0.1561138 0.1867933

df <- maf_format(f2r_eqtlg_Tubeneph2_merge)
eqtlg_Tubeneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome, 
                                               se1 = df$se.exposure, se2 = df$se.outcome, 
                                               MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                               N1 = 31684, N2 = 311, SNP = df$SNP) 
eqtlg_Tubeneph2_result

# nsnps PP.H0.abf  PP.H1.abf PP.H2.abf  PP.H3.abf PP.H4.abf
# 1  1812         0 0.06136816         0 0.04766794 0.8909639

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

df <- maf_format(f2r_gtexbld_Glomneph2_merge)
gtexbld_Glomneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome, 
                                                 se1 = df$se.exposure, se2 = df$se.outcome, 
                                                 MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                                 N1 = 670, N2 = 240, SNP = df$SNP) 
gtexbld_Glomneph2_result
# nsnps    PP.H0.abf PP.H1.abf    PP.H2.abf PP.H3.abf PP.H4.abf
# 1  1908 0.0008230391 0.7181235 0.0002037571 0.1776804 0.1031693

df <- maf_format(f2r_gtexbld_Tubeneph2_merge)
gtexbld_Tubeneph2_result <- coloc.analysis.quant(beta1 = df$beta.exposure, beta2 = df$beta.outcome,
                                                 se1 = df$se.exposure, se2 = df$se.outcome,
                                                 MAF1 = df$MAF1, MAF2 = df$MAF2, 
                                                 N1 = 670, N2 = 311, SNP = df$SNP) 
gtexbld_Tubeneph2_result
# nsnps    PP.H0.abf  PP.H1.abf    PP.H2.abf  PP.H3.abf PP.H4.abf
# 1  1882 3.512633e-05 0.03064857 2.009798e-05 0.01658326 0.9527129
