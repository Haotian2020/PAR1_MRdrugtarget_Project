# extract instruments from gwas data -------------------------------------------
# F2R gene location
# flanking region 100kb

# 38 build Chromosome 5: 76,716,126-76,735,770 
# 38build location to be extracted: 5: 76616126 - 76835770

# 37 build Chromosome 5: 76,011,951-76,031,595
# 37 build location to be extracted: 5: 75911951 - 76131595

source("fn-ld_clump_local")

# gtex blood dataset -----------------------------------------------------------
# 38 build

f2r_gtexbld_gwas = fread(paste0(rdsf_personal,"data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"))
f2r_gtexbld_gwas = f2r_gtexbld_gwas[grepl("ENSG00000181104", f2r_gtexbld_gwas$gene_id), ]
# f2rl1_gtexbld_gwas = f2r_gtexbld_gwas[grepl("ENSG00000164251", f2r_gtexbld_gwas$gene_id), ]
f2r_gtexbld_gwas = subset(f2r_gtexbld_gwas,pval_nominal<=5e-8)

# manually add SNP -------------------------------------------------------------

f2r_gtexbld_gwas$SNP = c("rs250737","rs250735","rs250734","rs250753")

# check with reference panel that each alt is the minor allele -----------------

split_cols <- strsplit(f2r_gtexbld_gwas$variant_id, "_")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos", "ref","alt","build")
f2r_gtexbld_gwas <- cbind(f2r_gtexbld_gwas, split_df)
f2r_gtexbld_gwas$chr <- 5

f2r_gtexbld_gwas_format = format_data(
  dat = data.frame(f2r_gtexbld_gwas),
  type = "outcome",
  snp_col = "SNP",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  beta_col = "slope",
  se_col = "slope_se",
  pval_col = "pval_nominal",
  chr_col = "chr",
  pos_col = "pos") %>% mutate(samplesize.outcome = 755)

f2r_gtex_str_exp = f2r_gtexbld_gwas_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(. ,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R Gtex str")

# weakly correlated SNP is same as the strongly independent one ----------------

# no siginificant signal in gtex v8 kidney tissue ------------------------------

write.table(f2r_gtex_str_exp, file = paste0(rdsf_personal,"data/par1/f2r_gtex_str_exp.csv"),
            sep= ',', row.names = F,col.names= T)


# write a file for coloc use ---------------------------------------------------

f2r_gtexbld_gwas = fread(paste0(rdsf_personal,"data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"))
f2r_gtexbld_gwas = f2r_gtexbld_gwas[grepl("ENSG00000181104", f2r_gtexbld_gwas$gene_id), ]
data.table::fwrite(f2r_gtexbld_gwas, paste0(rdsf_personal,"./data/par1/f2r_gtexbld_gwas.txt"))

f2r_gtexbld_gwas_mformat = MungeSumstats::format_sumstats(path=paste0(rdsf_personal,"data/par1/f2r_gtexbld_gwas.txt"),
                                                       ref_genome="GRCh38",
                                                       log_folder_ind=TRUE,
                                                       imputation_ind=TRUE,
                                                       log_mungesumstats_msgs=TRUE,
                                                       log_folder = paste0(rdsf_personal,"data/"),
                                                       dbSNP = 144,
                                                       return_data = F,
                                                       force_new=TRUE,
                                                       save_path = paste0(rdsf_personal,"data/par1/f2r_gtexbld_mformat.tsv.gz"))
