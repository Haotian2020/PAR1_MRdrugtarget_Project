# doing MR for all exp on atherosclerosis --- ----------------------------------

f2r_all_instruments = fread(paste0(rdsf_personal,"data/par1/f2r_all_instruments.csv"))

outcomes = c("cra","iga")

ukb_test_mr = uvmr(exp_name = "F2R ukb str", out_name = outcomes)

eqtl_test_mr = uvmr(exp_name = "F2R eqtl str", out_name = outcomes)

gtex_pos_mr = uvmr(exp_name = "F2R Gtex str",out_name = outcomes)

uvmr(exp_name = "F2R ukb str", out_name = c("ieu-a-302","ieu-a-301","ieu-a-299"))

uvmr(exp_name = "F2R eqtl str", out_name = c("ieu-a-302","ieu-a-301","ieu-a-299"))

uvmr(exp_name = "F2R Gtex str",out_name = c("ieu-a-302","ieu-a-301","ieu-a-299"))