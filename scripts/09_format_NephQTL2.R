# Script for formatting dataset from NephQTL2
# GRCh37 build in neph2
# chr5:76,011,868-76,031,606

76011868 - 500000
76031606 + 500000

# command
# gunzip eQTLs_chrAll_40_peers_cis1000kb_Glom.NephQTL2.txt.gz
# gunzip eQTLs_chrAll_50_peers_cis1000kb_Tube.NephQTL2.txt.gz 
# awk -F' ' '$1 == "ENSG00000181104" { print }' eQTLs_chrAll_40_peers_cis1000kb_Glom.NephQTL2.txt > f2r_Glom.neph2.txt

# awk -F' ' '$1 == "ENSG00000181104" { print }' eQTLs_chrAll_50_peers_cis1000kb_Tube.NephQTL2.txt > f2r_Tube.neph2.txt
# read f2r_Glom.neph2.txt

f2r_Glomneph2 = fread(paste0(rdsf_personal,"data/par1/download_res/f2r_Glomneph2.txt"))
f2r_Tubeneph2 = fread(paste0(rdsf_personal,"data/par1/download_res/f2r_Tube.neph2.txt"))
colnames(f2r_Glomneph2) = c("gene","SNP_ID","CHR","POS","REF","ALT","AltFreq","beta","se","t-stat","pval","FDR")
colnames(f2r_Tubeneph2) = c("gene","SNP_ID","CHR","POS","REF","ALT","AltFreq","beta","se","t-stat","pval","FDR")

write.table(f2r_Glomneph2, file = paste0(rdsf_personal,"data/par1/download_res/f2r_Glomneph2_head.txt"),
            sep= ',', row.names = F,col.names= T)
write.table(f2r_Tubeneph2, file = paste0(rdsf_personal,"data/par1/download_res/f2r_Tubeneph2_head.txt"),
            sep= ',', row.names = F,col.names= T)


f2r_Glomneph2_mformat = MungeSumstats::format_sumstats(path=paste0(rdsf_personal,"data/par1/download_res/f2r_Glomneph2_head.txt"),
                                                       ref_genome="GRCh37",
                                                       log_folder_ind=TRUE,
                                                       imputation_ind=TRUE,
                                                       log_mungesumstats_msgs=TRUE,
                                                       log_folder = paste0(rdsf_personal,"data/format_data"),
                                                       dbSNP = 144,
                                                       return_data = F,
                                                       force_new=TRUE,
                                                       save_path = paste0(rdsf_personal,"data/par1/download_res/f2r_Glomneph2_mformat.txt.gz"))

f2r_Tubeneph2_mformat = MungeSumstats::format_sumstats(path=paste0(rdsf_personal,"data/par1/download_res/f2r_Tubeneph2_head.txt"),
                                                       ref_genome="GRCh37",
                                                       log_folder_ind=TRUE,
                                                       imputation_ind=TRUE,
                                                       log_mungesumstats_msgs=TRUE,
                                                       log_folder = paste0(rdsf_personal,"data/format_data"),
                                                       dbSNP = 144,
                                                       return_data = F,
                                                       force_new=TRUE,
                                                       save_path = paste0(rdsf_personal,"data/par1/download_res/f2r_Tubeneph2_mformat.txt.gz"))

f2r_Glomneph2_mformat = fread(paste0(rdsf_personal,"data/par1/download_res/f2r_Glomneph2_mformat.txt.gz"))
f2r_Tubeneph2_mformat = fread(paste0(rdsf_personal,"data/par1/download_res/f2r_Tubeneph2_mformat.txt.gz"))

f2r_Glomneph2_format = format_data(f2r_Glomneph2_mformat,
                                   type = "outcome",
                                   header = TRUE,
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "ALTFREQ",
                                   effect_allele_col = "A2",
                                   other_allele_col = "A1",
                                   pval_col = "P",
                                   min_pval = 1e-1000,
                                   chr_col = "CHR",
                                   pos_col = "BP") %>% mutate(outcome = "F2R GLOM")

f2r_Tubeneph2_format = format_data(f2r_Tubeneph2_mformat,
                                   type = "outcome",
                                   header = TRUE,
                                   snp_col = "SNP",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "ALTFREQ",
                                   effect_allele_col = "A2",
                                   other_allele_col = "A1",
                                   pval_col = "P",
                                   min_pval = 1e-1000,
                                   chr_col = "CHR",
                                   pos_col = "BP")%>% mutate(outcome = "F2R TUBE")

write.table(f2r_Glomneph2_format, file = paste0(rdsf_personal,"data/format_data/f2r_Glomneph2_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)
write.table(f2r_Tubeneph2_format, file = paste0(rdsf_personal,"data/format_data/f2r_Tubeneph2_GWAS_tidy_outcome.csv"),
            sep= ',', row.names = F,col.names= T)