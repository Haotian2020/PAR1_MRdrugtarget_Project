# save all instruments in the same file ----------------------------------------

col_names = c(
  "SNP",
  "exposure",
  "chr.exposure",
  "pos.exposure",
  "effect_allele.exposure",
  "other_allele.exposure",
  "eaf.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "mr_keep.exposure",
  "pval_origin.exposure",
  "id.exposure"
)

f2r_all_instruments = data.frame(rbind(fread(paste0(rdsf_personal,"data/par1/f2r_kidney_str_exp.csv")) %>% dplyr::select(all_of(col_names)),
                                       fread(paste0(rdsf_personal,"data/par1/f2r_kidneytube_str_exp.csv")) %>% dplyr::select(all_of(col_names)),
                                       fread(paste0(rdsf_personal,"data/par1/f2r_ukb_exp.csv")) %>% dplyr::select(all_of(col_names)),
                                       fread(paste0(rdsf_personal,"data/par1/f2r_eqtl_exp.csv")) %>% dplyr::select(all_of(col_names)),
                                       fread(paste0(rdsf_personal,"data/par1/f2r_gtex_str_exp.csv")) %>% dplyr::select(all_of(col_names))))

write.table(f2r_all_instruments, file = paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"),
            sep= ',', row.names = F,col.names= T)

