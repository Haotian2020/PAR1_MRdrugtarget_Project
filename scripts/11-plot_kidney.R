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

res$exposure[res$exposure == "F2R ukb str"] <- "F2R (UKB-PPP)"
res$exposure[res$exposure == "F2R eqtl str"] <- "F2R (eQTLGen)"
res$exposure[res$exposure == "F2R Gtex str"] <- "F2R (GTEx)"
res$exposure[res$exposure == "F2R Gtex str"] <- "F2R (GTEx)"
res$exposure[res$exposure ==  "F2R ukb wk"] <- "F2R (UKB-PPP)"
res$exposure[res$exposure ==  "F2R eqtl wk"] <- "F2R (eQTLGen)"
res$exposure[res$exposure ==  "F2R kidney meta str"] <- "F2R (Susztaklab Kidney)"
res$exposure[res$exposure ==  "F2R tubule meta str"] <- "F2R (Susztaklab Tubule)"

res$method[res$method ==  "Wald ratio"] <- "WR"
res$method[res$method ==  "Inverse variance weighted"] <- "IVW"
res$method[res$method ==  "Inverse variance weighted (correlated)"] <- "IVW (correlated)"

data.table::fwrite(res, paste0(rdsf_personal,"./results/par1_res_drugdir.csv"))

# draw -------------------------------------------------------------------------
# define outcomes --------------------------------------------------------------

exp = c("F2R (UKB-PPP)","F2R (eQTLGen)","F2R (GTEx)","F2R (Susztaklab Kidney)", "F2R (Susztaklab Tubule)")

main_out = c("CKD","eGFR")

sec_out = c("NS (Finngen)",
            "uACR",
            "MA",
            "NS (Meta-analyzed)",
            "Microalbumin in urine" ,
            "Albumin")

# main results 1 ---------------------------------------------------------------

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
              exp = exp,
              out = "CKD",
              line_number = 0,
              xlabel = "OR (with 95% CI) for CKD per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(0.8,1,1.2,1.4),
              intervals = c(0.8,1.4),
              type = "binary",
              order = "outcome",
              make_na = "outcome")

pdf(paste0(rdsf_personal,"results/F2R on CKD forestplot.pdf"), width = 15, height = 3)
plot.new()
mtext("A)",side = 3,line = 2,adj = -0.06, cex = 1.5,padj = 0)
print(p)
dev.off()

# main results 2 ---------------------------------------------------------------

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
                exp = exp,
                out = "eGFR",
                line_number = 0,
                xlabel = "Beta (with 95% CI) for SD unit change in log(eGFR) per SD unit change in lower PAR1 protein/expression level",
                x_ticks = c(-0.1,-0.05,0,0.05),
                intervals = c(-0.1,0.2),
                type = "conti",
                order = "outcome",
                make_na = "outcome")

pdf(paste0(rdsf_personal,"results/F2R on eGFR forestplot.pdf"),width = 15, height = 3)
plot.new()
mtext("B)",side = 3,line = 2,adj = -0.06, cex = 1.5,padj = 0)
print(p)
dev.off()

# for 2nd outcomes -------------------------------------------------------------

conti_out = c("uACR", "Microalbumin in urine", "Albumin")
binary_out = c("NS (Finngen)", "NS (Meta-analyzed)", "MA")

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
              exp = exp,
              out = binary_out,
              line_number = 5,
              xlabel = "OR (with 95% CI) for binary outcomes per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(0.8,1,1.2,1.4),
              intervals = c(0.8,1.4),
              type = "binary",
              order = "outcome",
              make_na = c("outcome",3,4,5,7,8,9,10))

pdf(paste0(rdsf_personal,"results/F2R on other binary forestplot.pdf"),width = 16, height = 4)
plot.new()
mtext("A)",side = 3,line = 2,adj = -0.05, cex = 1.5,padj = 0)
print(p)
dev.off()

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method != "IVW (correlated)"),
              exp = exp,
              out = conti_out,
              line_number = c(5,9),
              xlabel = "Beta (with 95% CI) for continuous outcomes per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(-0.1,-0.05,0,0.05,0.1),
              intervals = c(-0.1,0.1),
              type = "conti",
              order = "outcome",
              make_na = c("outcome",2,3,4,5,7,8,9,11,12,13))

pdf(paste0(rdsf_personal,"results/F2R on other conti forestplot.pdf"),width = 16, height = 5)
plot.new()
mtext("B)",side = 3,line = 2,adj = -0.05, cex = 1.5,padj = 0)
print(p)
dev.off()


# for weakly associated SNP figures --------------------------------------------

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = "CKD",
              line_number = 0,
              xlabel = "OR (with 95% CI) for CKD per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(0.8,1,1.2,1.4),
              intervals = c(0.8,1.4),
              type = "binary",
              order = "outcome",
              make_na = "outcome")

pdf(paste0(rdsf_personal,"results/correlated F2R on CKD forestplot.pdf"),width = 15, height = 2)
plot.new()
mtext("A)",side = 3,line = 2,adj = -0.06, cex = 1.5,padj = 0)
print(p)
dev.off()

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = "eGFR",
              line_number = 0,
              xlabel = "Beta (with 95% CI) for SD unit change in log(eGFR) per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(-0.1,-0.05,0,0.05),
              intervals = c(-0.1,0.2),
              type = "conti",
              order = "outcome",
              make_na = "outcome")

pdf(paste0(rdsf_personal,"results/correlated F2R on eGFR forestplot.pdf"),width = 15, height = 2)
plot.new()
mtext("B)",side = 3,line = 2,adj = -0.06, cex = 1.5,padj = 0)
print(p)
dev.off()

# weakly associated on other kidney outcomes

p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = binary_out,
              line_number = 2,
              xlabel = "OR (with 95% CI) for binary outcomes per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(0.8,1,1.2,1.4),
              intervals = c(0.8,1.4),
              type = "binary",
              order = "outcome",
              make_na = c("outcome",4))

pdf(paste0(rdsf_personal,"results/correlated F2R on other binary forestplot.pdf"),width = 16, height = 4)
plot.new()
mtext("A)",side = 3,line = 2,adj = -0.05, cex = 1.5,padj = 0)
print(p)
dev.off()


p = uvmr_plot(dat = res %>% subset(exposure %in% exp & method == "IVW (correlated)"),
              exp = exp,
              out = conti_out,
              line_number = 1,
              xlabel = "Beta (with 95% CI) for continuous outcomes per SD unit change in lower PAR1 protein/expression level",
              x_ticks = c(-0.1,-0.05,0,0.05,0.1),
              intervals = c(-0.1,0.1),
              type = "conti",
              order = "outcome",
              make_na = c("exposure",3,4))

pdf(paste0(rdsf_personal,"results/correlated F2R on other conti forestplot.pdf"),width = 16, height = 5)
plot.new()
mtext("B)",side = 3,line = 2,adj = -0.05, cex = 1.5,padj = 0)
print(p)
dev.off()