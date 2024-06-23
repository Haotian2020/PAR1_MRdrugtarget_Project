# coloc function for quantitative trait ----------------------------------------

coloc.analysis.quant <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,SNP){
  
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study 
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1, snp=SNP)
  
  #type, quant (quantitative) for quantitative trait study 
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="quant",N=N2, snp=SNP)
  
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)
  
  #Format the data to save out.
  
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  
  #Make the filename and save out.
  return(df)
}

# coloc function for binary trait ----------------------------------------------
coloc.analysis <- function(beta1,beta2,se1,se2,MAF1,MAF2,N1,N2,s,SNP){
  
  #Convert the inputs in order to run in coloc function.
  #type, quant (quantitative) for pQTL study
  dataset1 <- list(beta=beta1, varbeta=se1^2, MAF=MAF1,type="quant", N=N1,snp=SNP)
  
  #type, cc (case coontrol) for binary study
  dataset2 <- list(beta=beta2, varbeta=se2^2,MAF=MAF2, type="cc",s=s,N=N2,snp=SNP)
  
  #Run the coloc analysis, setting the prior probabilities for association with each trait (p1, p2) and both traits together (p12) as 1E-5.
  #p1 prior probability a SNP is associated with trait 1, default 1e-4
  #p2 prior probability a SNP is associated with trait 2, default 1e-4
  #p12 prior probability a SNP is associated with both traits, default 1e-5
  result <- coloc.abf(dataset1, dataset2, p1=1e-4, p2=1e-4, p12=1e-5)  
  
  #Format the data to save out.
  
  #List into data frame.
  df <- data.frame(matrix(unlist(result$summary), nrow=1, byrow=T))
  
  #Label the columns in the data frame.
  names(df) <- c("nsnps", "PP.H0.abf",    "PP.H1.abf",    "PP.H2.abf",    "PP.H3.abf",    "PP.H4.abf")
  
  #Make the filename and save out.
  return(df)
  
}

# minor allele frequency format ------------------------------------------------
maf_format <- function(dat){
  dat = dat[,c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome","eaf.exposure","eaf.outcome")]
  dat <- dat[complete.cases(dat), ]
  dat$MAF1 <- NULL
  dat$MAF2 <- NULL
  for (j in 1:nrow(dat)) {
    if (dat$eaf.exposure[j] > 0.5) {
      dat$MAF1[j] = 1 - dat$eaf.exposure[j]
    } else {
      dat$MAF1[j] = dat$eaf.exposure[j]
    }
    if (dat$eaf.outcome[j] > 0.5) {
      dat$MAF2[j] = 1 - dat$eaf.outcome[j]
    } else {
      dat$MAF2[j] = dat$eaf.outcome[j]
    }
  }
  return(dat)
}