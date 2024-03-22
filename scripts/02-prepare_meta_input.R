# Calculate effective sample size for ukb + fg meta analysis -------------------

# UKB
# case 349
# control 459687
# calculate effective sample size ----------------------------------------------
4/(1/349 + 1/459687)
# 1394.941

# FG
# case 917
# control 428292
# calculate effective sample size ----------------------------------------------
4/(1/917 + 1/428292)
# 3660.163

# ukb is build 37; need liftover -----------------------------------------------
source("fn-perform_liftover.R")
ns_ukb = fread(paste0(rdsf_personal,"data/GWAS_imputed/df_ns_imputed.txt.gz"))
ns_ukb$Neff = 1394.941
ns_ukb_marker = ns_ukb[,c("SNP","CHR","BP")]
colnames(ns_ukb_marker) = c("rsid","chrom","pos")
ns_ukb_marker_hg38 = perform_liftover(ns_ukb_marker)
ns_ukb_hg3738 <- merge(ns_ukb, ns_ukb_marker_hg38, by.x = "SNP", by.y = "rsid_b37")
ns_ukb_hg3738 = dplyr::select(ns_ukb_hg3738,-c(CHR,BP,chr_b37,pos_b37))
ns_ukb_hg3738 = subset(ns_ukb_hg3738,nchar(ALLELE1) == 1& nchar(ALLELE0)==1)

# correct beta and se for ukb GWAS ---------------------------------------------
# mu = Ncase/(Ncase + Ncontrol)
459687+349
349/460036
0.0007586363*(1-0.0007586363)
0.0007580608
ns_ukb_hg3738$BETA = ns_ukb_hg3738$BETA/0.0007580608
ns_ukb_hg3738$SE = ns_ukb_hg3738$SE/0.0007580608

# finngen is build 38 ----------------------------------------------------------

ns_fg = fread(paste0(rdsf_personal,"data/summary_stats_finngen_R10_N14_NEPHROTICSYND.gz"))
ns_fg$Neff = 3660.163
ns_fg = subset(ns_fg,nchar(ref) == 1 & nchar(alt) ==1)

# only keep rsid start with rs -------------------------------------------------

ns_ukb_hg3738 = ns_ukb_hg3738[ns_ukb_hg3738$SNP %like% "^rs", ]
ns_fg = ns_fg[ns_fg$rsids %like% "^rs", ]

# write input files for meta analysis ------------------------------------------

write.table(ns_ukb_hg3738,paste0(rdsf_personal,"./data/ns_ukb.txt"),
            sep= ' ', row.names = F,col.names= T,quote = F)

write.table(ns_fg,paste0(rdsf_personal,"./data/ns_fg.txt"),
            sep= ' ', row.names = F,col.names= T,quote = F)
