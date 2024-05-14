# doing MR for all exp on kidney function --------------------------------------

f2r_all_instruments = fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"))

kidney_outcomes = c("ns_meta","ns","egfr_sd","ckd","uacr","ma",
                    "ukb-d-30500_irnt","ukb-d-30600_irnt")

# for strongly independent exposure --------------------------------------------
`%notin%` <- Negate(`%in%`)

ukb_eqtl_mr = rbind(uvmr(exp_name = "F2R ukb str", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns_meta","ukb-d-30500_irnt","ukb-d-30600_irnt")]),
                    uvmr(exp_name = "F2R eqtl str", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns")]) )

kidney_mr = rbind(uvmr(exp_name = "F2R kidney meta str", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns")]),
                  uvmr(exp_name = "F2R tubule meta str", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns")]))

gtex_mr = uvmr(exp_name = "F2R Gtex str", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns")])

# for weakly independent exposure ----------------------------------------------
# multiple correlated SNPs from ukb-ppp, eqtlgen and gtex cross-tissue----------

ukb_eqtl_mr_ld = rbind(uvmr_cor(exp_name = "F2R ukb wk", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns_meta","ukb-d-30500_irnt","ukb-d-30600_irnt")]),
                       uvmr_cor(exp_name = "F2R eqtl wk", out_name = kidney_outcomes[kidney_outcomes %notin% c("ns")]))

# save results -----------------------------------------------------------------

par1_res =  rbind(ukb_eqtl_mr %>% filter(method %in% c("Inverse variance weighted","Wald ratio")),
      kidney_mr %>% filter(method %in% c("Inverse variance weighted","Wald ratio")),
      gtex_mr %>% filter(method %in% c("Inverse variance weighted","Wald ratio")),
      ukb_eqtl_mr_ld
      ) %>% generate_odds_ratios()

# write results ----------------------------------------------------------------

data.table::fwrite(par1_res, paste0(rdsf_personal,"./results/par1_res.csv"))