# extract instruments from gwas data -------------------------------------------
# F2R gene location
# flanking region 100kb

# 38 build Chromosome 5: 76,716,126-76,735,770 
# 38build location to be extracted: 5: 76616126 - 76835770

# 37 build Chromosome 5: 76,011,951-76,031,595
# 37 build location to be extracted: 5: 75911951 - 76131595

source("fn-ld_clump_local")

# f2r gene expression from ukb-ppp ----------------------------------------
# ukb-ppp is 38 build

f2r_ukb = fread(paste0(rdsf_personal,"data/format_data/f2r_ukb_GWAS_tidy_outcome.csv"))

f2r_ukb_str_exp = data.frame(f2r_ukb) %>% 
  filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R ukb str")

f2r_ukb_wk_exp =  data.frame(f2r_ukb) %>% 
  filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.1) %>% 
  mutate(exposure = "F2R ukb wk")

# f2r gene expression from eqtlgen ---------------------------------------------
# filter with Linux commands first ---------------------------------------------
# eqtl gen 37 build

# awk -F' ' '$8 == "ENSG00000181104" { print }' eqtlgen_full_cis_summary_statistics.txt > f2r_eqtlgen.txt
# awk -F' ' '$8 == "ENSG00000164251" { print }' eqtlgen_full_cis_summary_statistics.txt > f2rl1_eqtlgen.txt
# awk -F' ' '$8 == "ENSG00000164220" { print }' eqtlgen_full_cis_summary_statistics.txt > f2rl2_eqtlgen.txt
# awk -F' ' '$8 == "ENSG00000127533" { print }' eqtlgen_full_cis_summary_statistics.txt > f2rl3_eqtlgen.txt

# Load and format the f2r dataset from eqtlg and  ------------------------------

f2r_eqtlg_dat = fread(file = paste0(rdsf_personal,
                                    "data/f2r_eqtlgen.txt"),header = F)
colnames(f2r_eqtlg_dat) = c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele","Zscore","Gene","GeneSymbol","GeneChr",
                            "GenePos","NrCohorts","NrSamples","FDR","BonferroniP")

# Load allele frequency dataset ------------------------------------------------

eqtlg_freq = fread(file = paste0(rdsf_personal,"data/eqtlg_allelefreq.txt"))

f2r_eqtlg_dat = merge(f2r_eqtlg_dat,eqtlg_freq,by = "SNP")

# caluculate eaf, beta and se --------------------------------------------------

f2r_eqtlg_dat$eaf = ifelse(f2r_eqtlg_dat$AlleleB == f2r_eqtlg_dat$AssessedAllele,f2r_eqtlg_dat$AlleleB_all,1-f2r_eqtlg_dat$AlleleB_all)
f2r_eqtlg_dat$beta = f2r_eqtlg_dat$Zscore/sqrt(2*f2r_eqtlg_dat$eaf*(1-f2r_eqtlg_dat$eaf)*(f2r_eqtlg_dat$NrSamples + f2r_eqtlg_dat$Zscore^2))
f2r_eqtlg_dat$se = 1/sqrt(2*f2r_eqtlg_dat$eaf*(1-f2r_eqtlg_dat$eaf)*(f2r_eqtlg_dat$NrSamples + f2r_eqtlg_dat$Zscore^2))

f2r_eqtlg_dat_format = TwoSampleMR::format_data(
  data.frame(f2r_eqtlg_dat),
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  samplesize_col = "NrSamples",
  chr_col = "SNPChr",
  pos_col = "SNPPos",
  min_pval = 1e-500)

f2r_eqtl_str_exp =  f2r_eqtlg_dat_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76131595 & pos.outcome>=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R eqtl str")

f2r_eqtl_wk_exp =  f2r_eqtlg_dat_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76131595 & pos.outcome>=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.1) %>% 
  mutate(exposure = "F2R eqtl wk")

# save data --------------------------------------------------------------------

write.table(rbind(f2r_ukb_str_exp,f2r_ukb_wk_exp), file = paste0(rdsf_personal,"data/par1/f2r_ukb_exp.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(rbind(f2r_eqtl_str_exp,f2r_eqtl_wk_exp), file = paste0(rdsf_personal,"data/par1/f2r_eqtl_exp.csv"),
            sep= ',', row.names = F,col.names= T)
