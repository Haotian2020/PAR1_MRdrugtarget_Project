ao <- available_outcomes()

uvmr <- function(exp_name, out_name){
  
  
  instruments <- data.table::fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"), data.table = FALSE)
  
  tmp_exp = instruments %>% dplyr::filter(exposure == exp_name)
  
  print(paste0("Exposure name is ", unique(tmp_exp$exposure)))
  
  if(exp_name %in% c("f2r","f2rl1")){
    
    tmp_exp = fread(paste0(rdsf_personal, "data/format_data/",exp_name,"_tophits.tsv"))
    
    print("data from whole set of eqtlgen")
  }
  
  
  mr_bin = data.frame()
  
  for(each in out_name){
    
    if(each %in% ao$id){
      
      print(paste0("Outcome name is ",each))
      print("Reading from IEU open GWAS database")
      
      out_dat <- extract_outcome_data(snps = tmp_exp$SNP, outcomes = each, proxies = F)
      
    } else if(each %in% c("ns_meta","ns","egfr_sd","ckd","bun_sd","uacr","ma",
                          "vte","dvt","aet",
                          "cra","iga")){
      
      print(paste0("outcome name is ",each))
      
      out_dat <- TwoSampleMR::read_outcome_data(snps = tmp_exp$SNP,
                                                filename = paste0(rdsf_personal,"data/format_data/",each,"_GWAS_tidy_outcome.csv"),
                                                sep = ",",
                                                phenotype_col = "outcome",
                                                snp_col = "SNP",
                                                beta_col = "beta.outcome",
                                                se_col = "se.outcome",
                                                eaf_col = "eaf.outcome",
                                                effect_allele_col = "effect_allele.outcome",
                                                other_allele_col = "other_allele.outcome",
                                                pval_col = "pval.outcome",
                                                samplesize_col = "samplesize.outcome")
      
    }
    
    if(!is.null(out_dat)){
      
      dat <- TwoSampleMR::harmonise_data(exposure_dat = tmp_exp, outcome_dat = out_dat,action = 2)
    
      set.seed(123)
    
      mr <- mr(dat)
      
      mr_bin <- rbind(mr_bin,mr)
      
    }
    

  }
  
  return(mr_bin)
}
