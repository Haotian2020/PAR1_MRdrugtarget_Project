convert_outcome_to_exposure_local <- function(out_dat, ignore_samplesize=F){
  
  tmp_dat = convert_outcome_to_exposure(out_dat)
  
  if (ignore_samplesize == F) {
    ss = out_dat %>%  dplyr::select(c(
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
    ss = out_dat %>%  dplyr::select(c("SNP", starts_with("chr"), starts_with("pos")))
    colnames(ss) = c("SNP", "chr.exposure", "pos.exposure")
  }
  
  exp_dat = merge(tmp_dat, ss, by = "SNP")
  return(exp_dat)
}