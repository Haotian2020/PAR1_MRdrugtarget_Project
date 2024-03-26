# doing LD clump in local

ld_clump_local =  function(out_dat,
                           threshold = 5e-8,
                           r2 = 0.001,
                           ignore_sampleszie = T) {
  if ("pval.outcome" %in% colnames(out_dat)) {
    tmp = subset(out_dat, pval.outcome < threshold) %>% TwoSampleMR::convert_outcome_to_exposure()
  } else
  {
    tmp = subset(out_dat, pval.exposure < threshold)
  }
  
  snps = ieugwasr::ld_clump(
    dplyr::tibble(rsid = tmp$SNP,
                  pval = tmp$pval.exposure),
    clump_kb = 10000,
    clump_r2 = r2,
    clump_p = 0.99,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = paste0(rdsf_personal, "data/1kg_eur/EUR")
  )
  
  ins = subset(tmp, SNP %in% snps$rsid)
  if (ignore_sampleszie == F) {
    ss = subset(out_dat, SNP %in% snps$rsid) %>%  dplyr::select(c(
      "SNP",
      starts_with("chr"),
      starts_with("pos"),
      starts_with("samplesize")
    ))
    colnames(ss) = c("SNP",
                     "chr.exposure",
                     "pos.exposure",
                     "samplesize.exposure")
  } else {
    ss = subset(out_dat, SNP %in% snps$rsid) %>%  dplyr::select(c("SNP", starts_with("chr"), starts_with("pos")))
    colnames(ss) = c("SNP", "chr.exposure", "pos.exposure")
  }
  
  ins = merge(ins, ss, by = "SNP")
  
  return(ins)
}
