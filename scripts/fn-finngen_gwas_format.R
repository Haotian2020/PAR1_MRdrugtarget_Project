FinngenFormat = function(finngwas_name,outcome_name,ncase,ncontrol){
  finngwas = vroom(paste0(rdsf_personal,"data/",finngwas_name))
  head(finngwas)
  finngwas_format = format_data(
    dat = finngwas,
    type = "outcome",
    header = TRUE,
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    eaf_col = "af_alt",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    min_pval = 1e-1000,
    chr_col = "#chrom",
    pos_col = "pos",
    log_pval = FALSE
  )
  
  finngwas_format$outcome = outcome_name
  finngwas_format$ncase.outcome = ncase
  finngwas_format$ncontrol.outcome = ncontrol
  
  print(nrow(finngwas_format))
  print(paste0("returning formated data of ",outcome_name," to dataframe"))
  return(finngwas_format)
}
