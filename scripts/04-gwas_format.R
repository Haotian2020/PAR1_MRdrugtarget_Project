rm(list=ls())
graphics.off()
`%notin%` <- Negate(`%in%`)

# Load libraries ---------------------------------------------------------------

source("lib setting.R")

# Specify paths ----------------------------------------------------------------

source("specify_paths.R")

# format gwas data from FinnGen ------------------------------------------------

# R10 Finngen NS
ns_gwas_format = FinngenFormat("summary_stats_finngen_R10_N14_NEPHROTICSYND.gz","NS Finngen R10",917,428292)

data.table::fwrite(ns_gwas_format, paste0(rdsf_personal,"./data/format_data/ns_GWAS_tidy_outcome.csv"))

# R10 Finngen VTE
vte_gwas_format = FinngenFormat("summary_stats_finngen_R10_I9_VTE.gz","VTE Finngen R10",21702,407507)

data.table::fwrite(vte_gwas_format, paste0(rdsf_personal,"./data/format_data/vte_GWAS_tidy_outcome.csv"))

# R10 Finngen AET
aet_gwas_format = FinngenFormat("summary_stats_finngen_R10_I9_ARTEMBTHR.gz","AET Finngen R10",1883,427326)

data.table::fwrite(aet_gwas_format, paste0(rdsf_personal,"./data/format_data/aet_GWAS_tidy_outcome.csv"))

# R10 Finngen DVT
dvt_gwas_format = FinngenFormat("summary_stats_finngen_R10_I9_PHLETHROMBDVTLOW.gz","DVT Finngen R10",6501,422708)

data.table::fwrite(dvt_gwas_format, paste0(rdsf_personal,"./data/format_data/dvt_GWAS_tidy_outcome.csv"))

# R10 Finngen MI
mi_gwas_format = FinngenFormat("summary_stats_finngen_R10_I9_MI_STRICT.gz","MI Finngen R10",26060,403149)

data.table::fwrite(mi_gwas_format, paste0(rdsf_personal,"./data/format_data/mi_GWAS_tidy_outcome.csv"))

# format gwas data from meta-analysis of NS UKB and FinnGen R10-----------------

ns_meta = fread(paste0(rdsf_personal,"data/meta_ns.txt")) %>% data.frame()

# calculate lambuda
p_value=ns_meta$`P-value`
z = qnorm(p_value/ 2)
round(median(z^2, na.rm = TRUE) / 0.454, 3)
# 0.989
colnames(ns_meta)[colnames(ns_meta) == "P-value"] <- "pval"

# check how many SNPs are common
summary(ns_meta$MarkerName%in%ns_ukb$SNP)
#     Mode    FALSE     TRUE 
# logical  8032040 10920780 
summary(ns_meta$MarkerName %in% c(ns_fg$rsids,ns_ukb$SNP))
# Mode    FALSE     TRUE 
# logical  1001787 17951033

ns_meta_format = TwoSampleMR::format_data(
  dat = ns_meta,
  type = "outcome",
  header = TRUE,
  snp_col = "MarkerName",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "pval") %>% mutate(outcome = "Nephrotic syndrome UKB+FinnGen R10")

data.table::fwrite(ns_meta_format, paste0(rdsf_personal,"./data/format_data/ns_meta_GWAS_tidy_outcome.csv"))

# format gwas data from CKDGen -------------------------------------------------

# ensure there is no error during reading data ---------------------------------
ckd_gwas = fread(paste0(rdsf_personal,"data/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz"))
egfr_gwas = fread(paste0(rdsf_personal,"data/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt"))
egfr_meta_gwas = fread(paste0(rdsf_personal,"data/metal_eGFR_meta_ea1.TBL.map.annot.gc"))
bun_gwas = fread(paste0(rdsf_personal,"data/BUN_overall_EA_YL_20171108_METAL1_nstud24.dbgap.txt.gz"))
uacr_gwas = vroom(paste0(rdsf_personal,"data/formatted_20180517-UACR_overall-EA-nstud_18-SumMac_400.tbl.rsid"))
ma_gwas = vroom(paste0(rdsf_personal,"data/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid"))

# CKD
ckd_gwas_format = format_data(
  ckd_gwas,
  type = "outcome",
  header = TRUE,
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  min_pval = 1e-1000,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  log_pval = FALSE
)
ckd_gwas_format$outcome = "CKD"

data.table::fwrite(ckd_gwas_format, paste0(rdsf_personal,"./data/format_data/ckd_GWAS_tidy_outcome.csv"))

# egfr sd
egfr_gwas_outcome_format = format_data(egfr_gwas,
                                       type = 'outcome',
                                       snp_col = "RSID",
                                       beta_col = "Effect",
                                       se_col = "StdErr",
                                       effect_allele_col = "Allele1",
                                       other_allele_col = "Allele2",
                                       eaf_col = "Freq1",
                                       pval_col = "P-value",
                                       samplesize_col = "n_total_sum")

egfr_gwas_outcome_format$outcome = "eGFR"

# SD unite transformation
egfr_sd_gwas_outcome_format = egfr_gwas_outcome_format
egfr_sd_gwas_outcome_format$beta.outcome = egfr_sd_gwas_outcome_format$beta.outcome/0.13
egfr_sd_gwas_outcome_format$se.outcome = egfr_sd_gwas_outcome_format$se.outcome/0.13

data.table::fwrite(egfr_sd_gwas_outcome_format, paste0(rdsf_personal,"./data/format_data/egfr_sd_GWAS_tidy_outcome.csv"))

# BUN
bun_gwas_outcome_format = format_data(bun_gwas,
                                      type = 'outcome',
                                      snp_col = "RSID",
                                      beta_col = "Effect",
                                      se_col = "StdErr",
                                      effect_allele_col = "Allele1",
                                      other_allele_col = "Allele2",
                                      eaf_col = "Freq1",
                                      pval_col = "P-value",
                                      samplesize_col = "n_total_sum")
bun_gwas_outcome_format$outcome = "BUN"
bun_gwas_outcome_format$samplesize.outcome = 243031

# SD unite transformation
bun_sd_gwas_outcome_format = bun_gwas_outcome_format
bun_sd_gwas_outcome_format$beta.outcome = bun_sd_gwas_outcome_format$beta.outcome/0.24
bun_sd_gwas_outcome_format$se.outcome = bun_sd_gwas_outcome_format$se.outcome/0.24

data.table::fwrite(bun_sd_gwas_outcome_format, paste0(rdsf_personal,"./data/format_data/bun_sd_GWAS_tidy_outcome.csv"))

# uacr
uacr_gwas_format = format_data(
  uacr_gwas,
  type = "outcome",
  header = TRUE,
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  min_pval = 1e-1000,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  log_pval = FALSE
)
uacr_gwas_format$outcome = "UACR"

data.table::fwrite(uacr_gwas_format, paste0(rdsf_personal,"./data/format_data/uacr_GWAS_tidy_outcome.csv"))

# ma
ma_gwas_format = format_data(
  dat = ma_gwas,
  type = "outcome",
  header = TRUE,
  snp_col = "RSID",
  beta_col = "Effect",
  se_col = "StdErr",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
  samplesize_col = "n_total_sum",
  min_pval = 1e-1000,
  chr_col = "Chr",
  pos_col = "Pos_b37",
  log_pval = FALSE
)
ma_gwas_format$outcome = "Microalbuminuria"
ma_gwas_format$ncase.outcome = 51861
ma_gwas_format$ncontrol.outcome = 297093
data.table::fwrite(ma_gwas_format, paste0(rdsf_personal,"./data/format_data/ma_GWAS_tidy_outcome.csv"))

# f2r gwas from ukb-ppp --------------------------------------------------------

f2r_gwas = data.frame()
for(i in 1:22){
  print(i)
  f2r_chr = vroom(paste0(rdsf_personal,"data/F2R_P25116_OID20691_v1_Inflammation/discovery_chr",i,"_F2R:P25116:OID20691:v1:Inflammation.gz"))
  match_chr = vroom(paste0(rdsf_shared,"data/UKB-PPP/snp-rsid-maps/olink_rsid_map_mac5_info03_b0_7_chr",i,"_patched_v2.tsv.gz"))
  f2r_chr_rsid = merge(f2r_chr,match_chr)
  f2r_gwas = rbind(f2r_gwas,f2r_chr_rsid)
}

f2r_gwas$p = 10^(-f2r_gwas$LOG10P)

f2r_gwas_format = format_data(
  f2r_gwas,
  type = "outcome",
  header = TRUE,
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "p",
  samplesize_col = "N",
  min_pval = 1e-1000,
  chr_col = "CHROM",
  pos_col = "POS38")

f2r_gwas_format$outcome = "F2R whole ukb"

data.table::fwrite(f2r_gwas_format, paste0(rdsf_personal,"./data/format_data/f2r_GWAS_tidy_outcome.csv"))
