# extract instruments from gwas data -------------------------------------------
# F2R gene location
# flanking region 100kb

# 38 build Chromosome 5: 76,716,126-76,735,770 
# 38build location to be extracted: 5: 76616126 - 76835770

# 37 build Chromosome 5: 76,011,951-76,031,595
# 37 build location to be extracted: 5: 75911951 - 76131595

# Susztaklab kidney meta and  Tubule-specifc -----------------------------------
# url: https://susztaklab.com/Kidney_eQTL/download.php
# F2R gene id ENSG00000181104
# b37/hg19
# kiney_glom = vroom(paste0(rdsf_personal,"data/par1/Kidney_eQTL.GlomsigeQTLsFormated.txt.gz"))
# f2r_kiney_glom = kiney_glom[grep(pattern = "ENSG00000181104",kiney_glom$GeneID),] %>% data.frame()
# no data from glom

source("fn-ld_clump_local")

kiney_meta = fread(paste0(rdsf_personal,"data/par1/Kidney_eQTL_Meta_S686_Significant.q0.01.txt.gz"))

kiney_tube = vroom(paste0(rdsf_personal,"data/par1/Kidney_eQTL.TubsigeQTLsFormated.txt.gz"))

f2r_kiney_meta = kiney_meta[grep(pattern = "ENSG00000181104",kiney_meta$GeneID),] %>% data.frame()

f2r_kiney_tube = kiney_tube[grep(pattern = "ENSG00000181104",kiney_tube$GeneID),] %>% data.frame()

split_cols <- strsplit(f2r_kiney_meta$SNP_POS, ":")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos")
f2r_kiney_meta <- cbind(f2r_kiney_meta, split_df)

# find these 1618 SNPs by using VEP in European pop ----------------------------

write.table(f2r_kiney_meta$RSID, file = paste0(rdsf_personal,"data/par1/f2r_kiney_meta_rsid.txt"),
            sep= ' ', row.names = F,col.names= F,quote = F)

kidneymeta_rsid_maf = fread(paste0(rdsf_personal,"data/par1/findeurmafforpar1kidney.txt"))

# qc step

kidneymeta_rsid_maf = subset(kidneymeta_rsid_maf,EUR_AF!="-")
kidneymeta_rsid_maf = kidneymeta_rsid_maf[,c("#Uploaded_variation","Allele","EUR_AF")]
colnames(kidneymeta_rsid_maf) = c("RSID","Allele","EUR_AF")
kidneymeta_rsid_maf = kidneymeta_rsid_maf[!duplicated(kidneymeta_rsid_maf$RSID),]

f2r_kiney_meta = left_join(f2r_kiney_meta,kidneymeta_rsid_maf,by = "RSID")
f2r_kiney_meta$eaf = as.numeric(f2r_kiney_meta$EUR_AF)
f2r_kiney_meta$chr = as.numeric(f2r_kiney_meta$chr)
f2r_kiney_meta$pos = as.numeric(f2r_kiney_meta$pos)

f2r_kiney_meta_format = format_data(
  dat = f2r_kiney_meta %>% data.frame(),
  type = "outcome",
  snp_col = "RSID",
  beta_col = "Effect2ALT",
  se_col = "SE",
  pval_col = "PVAL",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  eaf_col = "eaf",
  chr_col = "chr",
  pos_col = "pos")

f2r_kidney_str_exp = f2r_kiney_meta_format %>% 
  filter(chr.outcome == 5 & pos.outcome <=76131595 & pos.outcome>=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001, ignore_samplesize = T) %>% 
  mutate(exposure = "F2R kidney meta str")

# Susztaklab Tubule-specifc ----------------------------------------------------
split_cols <- strsplit(f2r_kiney_tube$SNP_Location, ":")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos")
f2r_kiney_tube <- cbind(f2r_kiney_tube, split_df)
f2r_kiney_tube$chr = 5

# find these 22 SNPs by using VEP ----------------------------------------------
write.table(f2r_kiney_tube$rsID, file = paste0(rdsf_personal,"data/par1/f2r_kiney_tube_rsid.txt"),
            sep= ' ', row.names = F,col.names= F,quote = F)

kidneytube_rsid_maf = fread(paste0(rdsf_personal,"data/par1/kidneytube_rsid_eurmaf.txt"))
kidneytube_rsid_maf = subset(kidneytube_rsid_maf,EUR_AF!="-")
kidneytube_rsid_maf = kidneytube_rsid_maf[,c("#Uploaded_variation","Allele","EUR_AF")]
colnames(kidneytube_rsid_maf) = c("rsID","Allele","EUR_AF")
kidneytube_rsid_maf = kidneytube_rsid_maf[!duplicated(kidneytube_rsid_maf$rsID),]

f2r_kiney_tube = left_join(f2r_kiney_tube,kidneytube_rsid_maf,by = "rsID")
f2r_kiney_tube$eaf = as.numeric(f2r_kiney_tube$EUR_AF)
f2r_kiney_tube$pos = as.numeric(f2r_kiney_tube$pos)

f2r_kiney_tube_format = format_data(
  dat = f2r_kiney_tube %>% data.frame(),
  type = "outcome",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "Std",
  pval_col = "Pvalue",
  effect_allele_col = "Alt",
  other_allele_col = "Ref",
  eaf_col = "eaf",
  chr_col = "chr",
  pos_col = "pos")

f2r_kidneytube_str_exp = f2r_kiney_tube_format %>% 
  filter(chr.outcome == 5 & pos.outcome <=76131595 & pos.outcome >=75911951) %>% 
  ld_clump_local(.,threshold = 5e-8, r2 = 0.001, ignore_samplesize = T) %>% 
  mutate(exposure = "F2R tubule meta str")

# No. weakly correlated SNP = 1; same as strongly independent ------------------

write.table(f2r_kidney_str_exp, file = paste0(rdsf_personal,"data/par1/f2r_kidney_str_exp.csv"),
            sep= ',', row.names = F,col.names= T)

write.table(f2r_kidneytube_str_exp, file = paste0(rdsf_personal,"data/par1/f2r_kidneytube_str_exp.csv"),
            sep= ',', row.names = F,col.names= T)
