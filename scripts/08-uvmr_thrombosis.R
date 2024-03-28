# doing MR for all exp on thrombosis function ----------------------------------

f2r_all_instruments = fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"))

thromb_outcomes = c("vte","aet","dvt")

# for strongly independent exposure --------------------------------------------

ukb_eqtl_pos_mr = rbind(uvmr(exp_name = "F2R ukb str", out_name = thromb_outcomes),
                    uvmr(exp_name = "F2R eqtl str", out_name = thromb_outcomes) )

kidney_pos_mr = rbind(uvmr(exp_name = "F2R kidney meta str",out_name = thromb_outcomes),
                  uvmr(exp_name = "F2R tubule meta str",out_name = thromb_outcomes) )

gtex_pos_mr = uvmr(exp_name = "F2R Gtex str",out_name = thromb_outcomes)

# for weakly independent exposure ----------------------------------------------

ukb_eqtl_pos_mr_ld = rbind(uvmr_cor(exp_name = "F2R ukb wk", out_name = thromb_outcomes),
                       uvmr_cor(exp_name = "F2R eqtl wk", out_name = thromb_outcomes) )



data.table::fwrite(rbind(ukb_eqtl_pos_mr%>%filter(method %in% c("Inverse variance weighted","Wald ratio"))%>%generate_odds_ratios,
                         gtex_pos_mr%>%filter(method %in% c("Inverse variance weighted","Wald ratio"))%>%generate_odds_ratios), 
                   paste0(rdsf_personal,"./results/par1_pos_abstract.csv"))