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
f2r_ukb = vroom(paste0(rdsf_personal,"data/format_data/f2r_GWAS_tidy_outcome.csv"))

f2r_ukb_str_exp = f2r_ukb %>% 
  filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R ukb str")

f2r_ukb_wk_exp =  f2r_ukb %>% 
  filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.1) %>% 
  mutate(exposure = "F2R ukb wk")

# f2r gene expression from eqtlgen ---------------------------------------------
# filter with Linux commands first ---------------------------------------------
# eqtl gen 37 build

# awk -F' ' '$8 == "ENSG00000181104" { print }' eqtlgen_full_cis_summary_statistics.txt > f2r_eqtlgen.txt

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
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  samplesize_col = "NrSamples",
  chr_col = "SNPChr",
  pos_col = "SNPPos")

f2r_eqtl_str_exp =  f2r_eqtlg_dat_format %>% 
  filter(chr.exposure == 5 & pos.exposure<=76131595 & pos.exposure>=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R eqtl str")

f2r_eqtl_wk_exp =  f2r_eqtlg_dat_format %>% 
  filter(chr.exposure == 5 & pos.exposure<=76131595 & pos.exposure>=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.1) %>% 
  mutate(exposure = "F2R eqtl wk")

# gtex blood dataset -----------------------------------------------------------

f2r_gtexbld_gwas = fread(paste0(rdsf_personal,"data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"))
f2r_gtexbld_gwas = f2r_gtexbld_gwas[grepl("ENSG00000181104", f2r_gtexbld_gwas$gene_id), ]
head(f2r_gtexbld_gwas)
f2r_gtexbld_gwas = subset(f2r_gtexbld_gwas,pval_nominal<=5e-8)
f2r_gtexbld_gwas$SNP = c("rs250737","rs250735","rs250734","rs250753")
# check with reference panel that each alt is the minor allele
split_cols <- strsplit(f2r_gtexbld_gwas$variant_id, "_")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos", "ref","alt","build")
f2r_gtexbld_gwas <- cbind(f2r_gtexbld_gwas, split_df)

f2r_gtexbld_gwas_format = format_data(
  dat = f2r_gtexbld_gwas,
  type = "exposure",
  snp_col = "SNP",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  beta_col = "slope",
  se_col = "slope_se",
  pval_col = "pval_nominal",
  chr_col = "chr",
  pos_col = "pos")

# no siginificant signal in gtex v8 kidney tissue ------------------------------

# Create the EAF variable by reasoning that if the minor allele is the references allele, then
# 1) the EAF is 1-MAF for our data 
# 2) otherwise, the EAF is the MAF

# my_data2$Minor.Allele.Global.Frequency=as.numeric(my_data2$Minor.Allele.Global.Frequency)
# my_data2$EAF <- my_data2$Minor.Allele.Global.Frequency
# my_data2$EAF <- ifelse(my_data2$match=="1", 1-my_data2$Minor.Allele.Global.Frequency, 
#                        my_data2$Minor.Allele.Global.Frequency)
# 
# my_data2$dbSNP=as.character(my_data2$dbSNP)

f2r_gtex_str_clump_local =  ld_clump(
  dplyr::tibble(rsid=f2r_gtexbld_gwas_format$SNP,
                pval=f2r_gtexbld_gwas_format$pval.exposure),
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p = 0.99,
  plink_bin = genetics.binaRies::get_plink_binary(),
  bfile = paste0(rdsf_personal,"data/1kg_eur/EUR"))
# str clump = weak clump

f2r_gtexbld_exp = f2r_gtexbld_gwas_format %>% subset(SNP == f2r_gtex_str_clump_local$rsid)%>% mutate("F2R eQTL GTex")




f2r_all_instruments = data.frame()
f2r_all_instruments = rbind()

write.table(f2r_all_instruments, file = paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"),
            sep= ',', row.names = F,col.names= T)
