# function for lines in the figures --------------------------------------------

generate_lines <- function(line_number) {
  
  lines <- list()
  
  if(all(is.na(line_number))){
    
    lines = NA
    
  }else{
    for (i in line_number) {
      key <- as.character(i * 1 + 2)
      lines[[key]] <- gpar(lty = 1, lwd = 1, col = "black")
    }
  }
  
  return(lines)
  
}

# function to plot mr results --------------------------------------------------

uvmr_plot <- function(dat, exp, out, line_number, box_size = 0.15, line_width = 1, 
                      xlabel, x_ticks, intervals, type, order, make_na) {
  # remove id.exposure and id.outcome columns ----------------------------------
  
  mydata  <- data.frame(dat) %>%
    select(-c("id.exposure", "id.outcome")) %>%
    subset(exposure %in% exp & outcome %in% out) %>%
    generate_odds_ratios()
  
  # format data ----------------------------------------------------------------
  
  mydata$beta <-
    format(round(mydata$b, digits = 3),
           nsmall = 3)
  mydata$CI <-
    paste0(format(round(mydata$lo_ci, digits = 3),
                  nsmall = 3),
           ", ",
           format(round(mydata$up_ci, digits = 3),
                  nsmall = 3))
  mydata$pvalue <-
    format(round(mydata$pval, digits = 3))
  
  mydata$OR <-
    format(round(mydata$or, digits = 3),
           nsmall = 3)
  mydata$ORCI <-
    paste0(format(round(mydata$or_lci95, digits = 3),
                  nsmall = 3),
           ", ",
           format(round(mydata$or_uci95, digits = 3),
                  nsmall = 3))
  
  
  mydata$method <-
    factor(
      mydata$method,
      levels = c(
        "IVW",
        "IVW (correlated)",
        "WR"
      )
    )
  
  mydata$outcome <-
    factor(
      mydata$outcome,
      levels = c(
        "NS (Finngen)",
        "NS (Meta-analyzed)",
        "eGFR",
        "CKD",
        "uACR",
        "MA",
        "Microalbumin in urine",
        "Albumin",
        "AET",
        "VTE",
        "DVT"
      )
    )
  
  mydata$exposure <-
    factor(
      mydata$exposure,
      levels = c(
        "PAR1 (UKB-PPP)",
        "F2R (eQTLGen)",
        "F2R (GTEx blood)",
        "F2R (Susztaklab Kidney)",
        "F2R (Susztaklab Tubule)"
      )
    )
  
  if(order == "exposure"){
    
    sorted_index <- order(mydata$exposure, mydata$outcome, mydata$method)
    
  }else if(order == "outcome"){
    
    sorted_index <- order(mydata$outcome, mydata$exposure, mydata$method)
    
  }else{
    sorted_index <- order(mydata$method, mydata$exposure, mydata$outcome)
  }
  
  mydata = mydata[sorted_index, ]
  
  rownames(mydata) <- seq_len(nrow(mydata))
    
  # make the table shown with the figure ---------------------------------------
  
  if(any(is.na(make_na))){
    
    print("no change to dataframe")
    
  }else if(length(make_na) == 1){
    
    if(make_na %in% c("exposure","outcome")){
      
      print(paste0("Leave one variable in column of ", make_na))
      
      mydata[-1, make_na] <- NA
    }
    
  }else{
    
    print("Delete specific place in datasets")
    
    mydata[c(make_na[2:length(make_na)]), make_na[1]] <- NA
    
  }
  
  mydata = data.frame(mydata)
  
  print(mydata)
  
  tabletext <- if(type == "conti") {
    cbind(
      c("Expsoure", as.character(mydata[, 'exposure'])),
      c("Outcome", as.character(mydata[, 'outcome'])),
      c("Approach", as.character(mydata[, 'method'])),
      c('Number of SNPs', as.character(mydata[, 'nsnp'])),
      c("Beta", as.character(mydata[, 'beta'])),
      c("95% CI", as.character(mydata[, 'CI'])),
      c("p-value", as.character(mydata[, 'pvalue']))
    )
  } else if(type == "binary"){
    cbind(
      c("Expsoure", as.character(mydata[, 'exposure'])),
      c("Outcome", as.character(mydata[, 'outcome'])),
      c("Approach", as.character(mydata[, 'method'])),
      c('Number of SNPs', as.character(mydata[, 'nsnp'])),
      c("OR", as.character(mydata[, 'OR'])),
      c("95% CI", as.character(mydata[, 'ORCI'])),
      c("p-value", as.character(mydata[, 'pvalue']))
    )
  }
  
  mean_values <- if(type == "conti") {
    as.numeric(rbind(NA, cbind(mydata[, 'beta'])))
  } else if(type == "binary") {
    as.numeric(rbind(NA, cbind(mydata[, 'or'])))
  }
  
  lower_values <- if(type == "conti") {
    as.numeric(rbind(NA, cbind(mydata[, 'lo_ci'])))
  } else if(type == "binary") {
    as.numeric(rbind(NA, cbind(mydata[, 'or_lci95'])))
  }
  
  upper_values <- if(type == "conti") {
    as.numeric(rbind(NA, cbind(mydata[, 'up_ci'])))
  } else if(type == "binary") {
    as.numeric(rbind(NA, cbind(mydata[, 'or_uci95'])))
  }
  
  zero_value <- if(type == "conti") 0.0 else if(type == "binary") 1.0
    
  p <- forestplot(
    tabletext,
    graph.pos = 4,
    mean = mean_values,
    lower = lower_values,
    upper = upper_values,
    new_page = F,
    psignif = 0.05,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1),
      title = gpar(cex = 1)
    ),
    hrzl_lines = if(all(is.na(line_number))) FALSE else generate_lines(line_number),
    boxsize = box_size,
    line.margin = 0.1,
    lty.ci = 1,
    col = fpColors(box = "black", lines = "darkgray"),
    lwd.ci = line_width,
    ci.vertices = T,
    ci.vertices.height = 0.15,
    graphwidth = unit(150, "mm"),
    is.summary = c(T, rep(F, nrow(tabletext))),
    colgap = unit (5, "mm"),
    zero = zero_value,
    xticks = x_ticks,
    clip = intervals,
    xlab = xlabel
  )
  
  return(p)
    
}
