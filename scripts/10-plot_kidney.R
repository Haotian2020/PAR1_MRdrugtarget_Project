# plot for kidney function outcomes --------------------------------------------
# reading results and format data ----------------------------------------------

res = fread(paste0(rdsf_personal,"results/par1_res.csv")) %>% data.frame() %>% 
  subset(method %in% c( "Inverse variance weighted","Inverse variance weighted (correlated)","Wald ratio"))

# change the direction to the drug effect
res$b = -res$b
res = res %>% generate_odds_ratios()

data.table::fwrite(res, paste0(rdsf_personal,"./results/par1_res_drugdir.csv"))

# draw -------------------------------------------------------------------------
# Options: 
# "F2R kidney meta str" "F2R tubule meta str" "F2R ukb str"        
# "F2R ukb wk"          "F2R eqtl str"        "F2R eqtl wk"        
# "F2R Gtex str"        "F2R Gtex cross str"  "F2R Gtex cross wk"  

# define outcomes --------------------------------------------------------------

unique_out <- unique(res[res$exposure == "F2R ukb str","outcome"])
`%notin%` <- Negate(`%in%`)
out_names <- unique_out[unique_out %notin% c("NS Finngen R10")]
out_noukb <- unique_out[unique_out %notin% c("Nephrotic syndrome UKB+FinnGen R10","Microalbumin in urine || id:ukb-d-30500_irnt","Albumin || id:ukb-d-30600_irnt")]

# ------------------------------------------------------------------------------

ukb_exp = c("F2R ukb str","F2R ukb wk")

for(i in ukb_exp){
  print(i)
  p = uvmr_plot(dat = res,
                exp = i,
                out = out_noukb,
                line_number = 0,
                xlabel = "Beta (with 95% CI) for each binary outcome per SD unit change in protein level",
                x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2),
                intervals = c(-0.2,0.2))
  
  pdf(paste0(rdsf_personal,"results/",i," on kidney function forestplot.pdf"),width = 18, height = 5)
  plot.new()
  print(p)
  dev.off()
}

# ------------------------------------------------------------------------------

noukb_exp = c("F2R kidney meta str", "F2R tubule meta str", 
              "F2R eqtl str","F2R eqtl wk",
              "F2R Gtex str","F2R Gtex cross str","F2R Gtex cross wk")
              
for(j in noukb_exp){
  print(j)
  p = uvmr_plot(dat = res,
                exp = j,
                out = out_names,
                line_number = 0,
                xlabel = "Beta (with 95% CI) for each binary outcome per SD unit change in mRNA level",
                x_ticks = c(-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2),
                intervals = c(-0.2,0.2))
  
  pdf(paste0(rdsf_personal,"results/", j," on kidney function forestplot.pdf"),width = 18, height = 7)
  plot.new()
  print(p)
  dev.off()
}
