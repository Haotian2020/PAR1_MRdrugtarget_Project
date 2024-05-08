# plot for thrombosis outcomes -------------------------------------------------
# reading results and format data ----------------------------------------------

res_pos = fread( paste0(rdsf_personal,"./results/par1_pos_res.csv")) %>% data.frame() %>% 
  subset(method %in% c( "Inverse variance weighted","Inverse variance weighted (correlated)","Wald ratio"))

# change the direction to the drug effect --------------------------------------
res_pos$b = -res_pos$b
res_pos = res_pos %>% generate_odds_ratios()

# format res_pos dataset -------------------------------------------------------
res_pos$outcome[res_pos$outcome == "VTE Finngen R10" ] <- "VTE"
res_pos$outcome[res_pos$outcome == "AET Finngen R10" ] <- "AET"
res_pos$outcome[res_pos$outcome == "DVT Finngen R10" ] <- "DVT"

res_pos$exposure[res_pos$exposure == "F2R ukb str"] <- "F2R (UKB-PPP)"
res_pos$exposure[res_pos$exposure == "F2R eqtl str"] <- "F2R (eQTLGen)"
res_pos$exposure[res_pos$exposure == "F2R Gtex str"] <- "F2R (GTEx)"
res_pos$exposure[res_pos$exposure == "F2R Gtex str"] <- "F2R (GTEx)"
res_pos$exposure[res_pos$exposure ==  "F2R ukb wk"] <- "F2R (UKB-PPP)"
res_pos$exposure[res_pos$exposure ==  "F2R eqtl wk"] <- "F2R (eQTLGen)"
res_pos$exposure[res_pos$exposure ==  "F2R kidney meta str"] <- "F2R (Susztaklab Kidney)"
res_pos$exposure[res_pos$exposure ==  "F2R tubule meta str"] <- "F2R (Susztaklab Tubule)"

res_pos$method[res_pos$method ==  "Wald ratio"] <- "WR"
res_pos$method[res_pos$method ==  "Inverse variance weighted"] <- "IVW"
res_pos$method[res_pos$method ==  "Inverse variance weighted (correlated)"] <- "IVW (correlated)"

data.table::fwrite(res_pos, paste0(rdsf_personal,"./results/par1_pos_res_drugdir.csv"))

# draw -------------------------------------------------------------------------
# Options: 
# "F2R kidney meta str" "F2R tubule meta str" "F2R ukb str"        
# "F2R ukb wk"          "F2R eqtl str"        "F2R eqtl wk"        
# "F2R Gtex str"        "F2R Gtex cross str"  "F2R Gtex cross wk"  

str_exp = c(
  "F2R ukb str",
  "F2R kidney meta str",
  "F2R tubule meta str",
  "F2R eqtl str",
  "F2R Gtex str",
  "F2R Gtex cross str"
)

wk_exp = c(
  "F2R ukb wk" ,
  "F2R eqtl wk",
  "F2R Gtex cross wk"
)

out_name = unique(res_pos[res_pos$exposure == "F2R ukb str","outcome"])



for(k in str_exp){
  print(k)
  
  p = uvmr_plot(dat = res_pos,
                exp = k,
                out = out_name,
                line_number = 0,
                xlabel = "Beta (with 95% CI) for each binary outcome per SD unit change in protein level or mRNA level",
                x_ticks = c(-0.1,0,0.1,0.2,0.3,0.4),
                intervals = c(-0.1,0.4))
  
  pdf(paste0(rdsf_personal,"results/", k," on pos forestplot.pdf"),width = 18, height = 5)
  plot.new()
  print(p)
  dev.off()
}

for(l in wk_exp){
  print(l)
  
  p = uvmr_plot(dat = res_pos,
                exp = l,
                out = out_name,
                line_number = 0,
                xlabel = "Beta (with 95% CI) for each binary outcome per SD unit change in protein level or mRNA level",
                x_ticks = c(-0.1,0,0.1,0.2,0.3,0.4,0.5),
                intervals = c(-0.1,0.5))
  
  pdf(paste0(rdsf_personal,"results/", l," on pos forestplot.pdf"),width = 18, height = 5)
  plot.new()
  print(p)
  dev.off()
}
