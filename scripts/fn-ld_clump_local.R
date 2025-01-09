ld_clump_local <- function(out_dat,
                           threshold = 5e-8,
                           r2 = 0.001,
                           ignore_samplesize = FALSE) {
  # Step 1: Filter based on p-value and convert outcome if applicable
  if ("pval.outcome" %in% colnames(out_dat)) {
    tmp <- subset(out_dat, pval.outcome < threshold) %>% 
      TwoSampleMR::convert_outcome_to_exposure()
  } else if ("pval.exposure" %in% colnames(out_dat)) {
    tmp <- subset(out_dat, pval.exposure < threshold)
  } else {
    stop("There is no p-value column in the provided dataset.")
  }
  
  # Step 2: Perform LD clumping
  snps <- ieugwasr::ld_clump(
    dplyr::tibble(rsid = tmp$SNP, pval = tmp$pval.exposure),
    clump_kb = 10000,
    clump_r2 = r2,
    clump_p = 0.99,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = paste0(rdsf_personal, "data/1kg_eur/EUR")
  )
  
  # Step 3: Filter for SNPs that survived clumping
  ins <- subset(tmp, SNP %in% snps$rsid)
  
  # Step 4: Select additional columns based on ignore_samplesize flag
  if (ignore_samplesize == FALSE) {
    chr_cols <- any(grepl("^chr", colnames(out_dat)))
    pos_cols <- any(grepl("^pos", colnames(out_dat)))
    samplesize_cols <- any(grepl("^samplesize", colnames(out_dat)))
    
    # Check for the presence of required columns
    if (chr_cols && pos_cols && samplesize_cols) {
      ss <- subset(out_dat, SNP %in% snps$rsid) %>%
        dplyr::select(SNP, starts_with("chr"), starts_with("pos"), starts_with("samplesize"))
      colnames(ss) <- c("SNP", "chr.exposure", "pos.exposure", "samplesize.exposure")
      ins <- merge(ins, ss, by = "SNP")
    }
  } else {
    chr_cols <- any(grepl("^chr", colnames(out_dat)))
    pos_cols <- any(grepl("^pos", colnames(out_dat)))
    
    # Check for the presence of required columns
    if (chr_cols && pos_cols) {
      ss <- subset(out_dat, SNP %in% snps$rsid) %>%
        dplyr::select(SNP, starts_with("chr"), starts_with("pos"))
      colnames(ss) <- c("SNP", "chr.exposure", "pos.exposure")
      ins <- merge(ins, ss, by = "SNP")
    }
  }
  
  return(ins)
}
