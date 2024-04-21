# doing MR for all exp on kindey function --------------------------------------

f2r_all_instruments = fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"))

kidney_outcomes = c("ns_meta","egfr_sd","ckd",
                    "bun_sd","uacr","ma",
                    "ukb-d-30500_irnt","ukb-d-30600_irnt")

# for strongly independent exposure --------------------------------------------

ukb_eqtl_mr = rbind(uvmr(exp_name = "F2R ukb str", out_name = kidney_outcomes),
                    uvmr(exp_name = "F2R eqtl str", out_name = kidney_outcomes) )

kidney_mr = rbind(uvmr(exp_name = "F2R kidney meta str", out_name = kidney_outcomes),
                  uvmr(exp_name = "F2R tubule meta str", out_name = kidney_outcomes))

gtex_mr = uvmr(exp_name = "F2R Gtex str", out_name = kidney_outcomes)

# for weakly independent exposure ----------------------------------------------

ukb_eqtl_mr_ld = rbind(uvmr_cor(exp_name = "F2R ukb wk", out_name = kidney_outcomes),
                       uvmr_cor(exp_name = "F2R eqtl wk", out_name = kidney_outcomes) )

ukb_eqtl_mr = rbind(uvmr_cor(exp_name = "F2R ukb str", out_name = kidney_outcomes),
                    uvmr_cor(exp_name = "F2R eqtl str", out_name = kidney_outcomes) )

kidney_mr = rbind(uvmr_cor(exp_name = "F2R kidney meta str",out_name = kidney_outcomes),
                  uvmr_cor(exp_name = "F2R tubule meta str",out_name = kidney_outcomes) )

gtex_mr = uvmr_cor(exp_name = "F2R Gtex str",out_name = kidney_outcomes)

# write results ----------------------------------------------------------------

data.table::fwrite(rbind(ukb_eqtl_mr%>%filter(method %in% c("Inverse variance weighted","Wald ratio"))%>%generate_odds_ratios,
                         gtex_mr%>%filter(method %in% c("Inverse variance weighted","Wald ratio"))%>%generate_odds_ratios), 
                   paste0(rdsf_personal,"./results/par1_kidney_abstract.csv"))


uvmr(exp_name = "F2R ukb str", out_name = c("ieu-b-38", "ieu-b-39"))

uvmr(exp_name = "F2R eqtl str", out_name = c("ieu-b-38","ieu-b-39"))

uvmr(exp_name = "F2R kidney meta str", out_name = c("ieu-b-38","ieu-b-39"))