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
exp = c("F2R (UKB-PPP)","F2R (eQTLGen)","F2R (GTEx)","F2R (Susztaklab Kidney)", "F2R (Susztaklab Tubule)")

out = c("VTE","DVT","AET")

p = uvmr_plot(dat = res_pos %>% subset(method %in% c("WR","IVW")),
              exp = exp,
              out = out,
              line_number = c(3,6,9,12),
              xlabel = "OR (with 95% CI) for each binary outcome per SD unit change in protein level or mRNA level",
              x_ticks = c(0.5,1.0,1.5),
              intervals = c(0.5,1.5),
              type = "binary",
              order = "outcome",
              make_na = c("exposure",2,3,5,6,8,9,11,12,14,15))

pdf(paste0(rdsf_personal,"results/all on pos forestplot.pdf"),width = 15, height = 6)
plot.new()
print(p)
dev.off()

p = uvmr_plot(dat = res_pos %>% subset(method == "IVW (correlated)" & exposure %in% exp),
              exp = exp,
              out = out,
              line_number = 3,
              xlabel = "OR (with 95% CI) for VTE per SD unit change in protein level or mRNA level",
              x_ticks = c(0.6,0.8,1.0,1.2),
              intervals = c(0.6,1.2),
              type = "binary",
              order = "outcome",
              make_na = c("exposure",2,3,5,6))
pdf(paste0(rdsf_personal,"results/correlated F2R on pos forestplot.pdf"),width = 15, height = 3)
plot.new()
print(p)
dev.off()

