# plot for kidney function outcomes --------------------------------------------
# reading results and format data ----------------------------------------------

res = fread(paste0(rdsf_personal,"results/par1_res.csv")) %>% data.frame() %>% 
  subset(method %in% c( "Inverse variance weighted","Inverse variance weighted (correlated)","Wald ratio"))

# change the direction to the drug effect
res$b = -res$b
res = res %>% generate_odds_ratios()

# format res dataset -----------------------------------------------------------

res = split_outcome(res)
res$outcome[res$outcome == "ckd"] <- "CKD"
res$outcome[res$outcome == "eGFR (CKDGen2019)"] <- "eGFR"
res$outcome[res$outcome == "UACR"] <- "uACR"
res$outcome[res$outcome == "Microalbuminuria"] <- "MA"
res$outcome[res$outcome == "Nephrotic syndrome UKB+FinnGen R10" ] <- "NS (Meta-analyzed)"
res$outcome[res$outcome == "NS Finngen R10" ] <- "NS (Finngen)"

res$exposure[res$exposure == "F2R ukb str"] <- "PAR1 (UKB-PPP)"
res$exposure[res$exposure == "F2R eqtl str"] <- "F2R (eQTLGen)"
res$exposure[res$exposure == "F2R Gtex str"] <- "F2R (GTEx blood)"
res$exposure[res$exposure ==  "F2R ukb wk"] <- "PAR1 (UKB-PPP)"
res$exposure[res$exposure ==  "F2R eqtl wk"] <- "F2R (eQTLGen)"
res$exposure[res$exposure ==  "F2R kidney meta str"] <- "F2R (Susztaklab Kidney)"
res$exposure[res$exposure ==  "F2R tubule meta str"] <- "F2R (Susztaklab Tubule)"

res$method[res$method ==  "Wald ratio"] <- "WR"
res$method[res$method ==  "Inverse variance weighted"] <- "IVW"
res$method[res$method ==  "Inverse variance weighted (correlated)"] <- "IVW (correlated)"

data.table::fwrite(res, paste0(rdsf_personal,"./results/par1_res_drugdir.csv"))

# draw -------------------------------------------------------------------------
# define outcomes --------------------------------------------------------------

exp = c("PAR1 (UKB-PPP)","F2R (eQTLGen)","F2R (GTEx blood)","F2R (Susztaklab Kidney)", "F2R (Susztaklab Tubule)")

# all outcomes
conti_out = c("eGFR","uACR")

binary_out = c("CKD","MA","NS (Meta-analyzed)","NS (Finngen)")

sec_out = c("Albumin")

# main binary ------------------------------------------------------------------
# CKD MA NS
p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
              exp = exp,
              out = binary_out,
              line_number = c(5,9,10),
              xlabel = expression(OR ~ "(with 95% CI) for kidney diseases per SD unit change in lower PAR1 protein or" ~ italic("F2R") ~"expression level"),
              x_ticks = c(0.7,1,1.3,1.6),
              intervals = c(0.7,1.6),
              type = "binary",
              order = "outcome",
              make_na = c("outcome",2,3,4,5,7,8,9,12,13,14))

pdf(paste0(rdsf_personal,"results/F2R on kidney binary forestplot.pdf"), width = 15, height = 7)
plot.new()
print(p)
dev.off()

# main conti -------------------------------------------------------------------

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
                exp = exp,
                out = conti_out,
                line_number = 5,
                xlabel = expression("Beta (with 95% CI) for SD unit change in log(eGFR) or log(uACR) per SD unit change in lower PAR1 protein or" ~ italic("F2R")~"expression level"),
                x_ticks = c(-0.075,-0.05,-0.025,0,0.025,0.05,0.075),
                intervals = c(-0.075,0.075),
                type = "conti",
                order = "outcome",
                make_na = c("outcome",2,3,4,5,7,8,9,10))

pdf(paste0(rdsf_personal,"results/F2R on conti forestplot.pdf"),width = 15, height = 6)
plot.new()
print(p)
dev.off()

# for 2nd outcomes -------------------------------------------------------------

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
              exp = exp,
              out = sec_out,
              line_number = NA,
              xlabel = expression("Beta (with 95% CI) for SD unit change in albumin per SD unit change in" ~ italic("F2R")~"expression level"),
              x_ticks = c(-0.1,-0.05,0,0.05,0.1),
              intervals = c(-0.1,0.1),
              type = "conti",
              order = "outcome",
              make_na = "outcome")

pdf(paste0(rdsf_personal,"results/F2R on albumin forestplot.pdf"),width = 15, height = 3)
plot.new()
print(p)
dev.off()

# for weakly associated SNP figures --------------------------------------------
# binary
p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = c("CKD", "MA","NS (Meta-analyzed)","NS (Finngen)"),
              line_number = c(2, 3),
              xlabel = expression(OR ~ "(with 95% CI) for kidney diseases per SD unit change in lower PAR1 protein or" ~ italic("F2R") ~"expression level"),
              x_ticks = c(0.8,1,1.2,1.4),
              intervals = c(0.8,1.4),
              type = "binary",
              order = "outcome",
              make_na = c("outcome",2))

pdf(paste0(rdsf_personal,"results/correlated F2R on kidney binary forestplot.pdf"),width = 15, height = 5)
plot.new()
print(p)
dev.off()

# conti 
p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = conti_out,
              line_number = 2,
              xlabel = expression("Beta (with 95% CI) for SD unit change in log(eGFR) or log(uACR) per SD unit change in lower PAR1 protein or" ~ italic("F2R") ~"expression level"),
              x_ticks = c(-0.1,-0.05,0,0.05),
              intervals = c(-0.1,0.2),
              type = "conti",
              order = "outcome",
              make_na = c("outcome",2,4))

pdf(paste0(rdsf_personal,"results/correlated F2R on conti forestplot.pdf"),width = 15, height = 3)
plot.new()
print(p)
dev.off()

# weakly associated on albumin

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = c("Albumin"),
              line_number = NA,
              xlabel = expression("Beta (with 95% CI) for SD unit change in albumin per SD unit change in lower PAR1 protein or" ~ italic("F2R") ~"expression level"),
              x_ticks = c(-0.1,-0.05,0,0.05,0.1),
              intervals = c(-0.1,0.1),
              type = "conti",
              order = "outcome",
              make_na = "outcome")

pdf(paste0(rdsf_personal,"results/correlated F2R on albumin forestplot.pdf"),width = 16, height = 2)
plot.new()
print(p)
dev.off()