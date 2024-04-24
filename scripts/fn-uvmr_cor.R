# functions for MR with correlated instruments ---------------------------------

# function formatting ivw and egger --------------------------------------------

mr_ld_format = function(mrivw, mregger){
  dfivw <- c("outcome" = mrivw@Outcome, "exposure"= mrivw@Exposure,"method" = "Inverse variance weighted (correlated)", "b"=mrivw@Estimate, "se"=mrivw@StdError, 
             "pval"=mrivw@Pvalue, "nsnp"=mrivw@SNPs,"id.exposure" = "MRaccountLD","id.outcome" = "MRaccountLD")
  
  dfegger <- c("outcome" = mregger@Outcome, "exposure"= mregger@Exposure,"method" = "MR Egger (correlated)", "b"=mregger@Estimate, "se"=mregger@StdError.Est, 
               "pval"=mregger@Pvalue.Est, "nsnp"=mregger@SNPs,"id.exposure" = "MRaccountLD","id.outcome" = "MRaccountLD")
  
  res_ld_mat =  rbind(as.data.frame(t(as.matrix(dfivw))),
                      as.data.frame(t(as.matrix(dfegger))))%>% 
    mutate(b = as.numeric(b), 
           se = as.numeric(se))
  res_ld <- res_ld_mat
  return(res_ld)}

# function using MR package to do correlated MR --------------------------------

MRforCorrelated = function(dat_har,i){
  mrivw = MendelianRandomization::mr_ivw(dat_har[[i]], correl = TRUE)
  mregger = MendelianRandomization::mr_egger(dat_har[[i]], correl = TRUE)
  mr_ld = mr_ld_format(mrivw,mregger)
  return(mr_ld)
}

# function doing MR with correlated instruments --------------------------------

uvmr_cor <- function(exp_name, out_name){
  
  instruments <- data.table::fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"), data.table = FALSE)
  
  tmp_exp = instruments %>% filter(exposure == exp_name)
  
  print(paste0("Exposure name is ", unique(tmp_exp$exposure)))
  
  mr_ld_bin = data.frame()
  
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
  
    dat <- TwoSampleMR::harmonise_data(exposure_dat = tmp_exp, outcome_dat = out_dat,action = 2)
    
    if(nrow(dat) >= 2){
      
      dat2 = dat_to_MRInput(dat,get_correlations = T,pop = "EUR")
      
      set.seed(123)
      
      mr_ld <- MRforCorrelated(dat2,1)
      
      mr_ld_bin <- rbind(mr_ld_bin,mr_ld)
      
    }else{
      
      mr_ld <- mr(dat)
      
      mr_ld_bin <- rbind(mr_ld_bin,mr_ld)
      
    }
    
  }
  
  return(mr_ld_bin)
}