devtools::install_github("myles-lewis/locuszoomr")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
# library(AnnotationHub)
library(LDlinkR)
library(locuszoomr)

source("ld_token.R")

# F2R location
# chr5:76716126-76735770
# only draw 300kb window -------------------------------------------------------
# ukb
f2r_ukb = fread(paste0(rdsf_personal,"data/format_data/f2r_ukb_GWAS_tidy_outcome.csv")) %>% data.frame()
f2r_ukb_ch5 = subset(f2r_ukb, chr.outcome ==5)
f2r_ukb_ch5 <- convert_outcome_to_exposure_local(f2r_ukb_ch5,ignore_samplesize = T)
f2r_ukb_coloc_snp = f2r_ukb_ch5[grep("^rs", f2r_ukb_ch5$SNP), ]

loc_f2r_ukb <- locus(data = f2r_ukb_coloc_snp,
                     gene = 'F2R',
                     flank = 1.5e5,
                     chrom = "chr.outcome",
                     pos = "pos.outcome",
                     p = "pval.outcome",
                     labs = "SNP",
                     index_snp = "rs168753",
                     ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_ukb)
loc_f2r_ukb <- link_LD(loc_f2r_ukb, token = link_ld_token)
summary(loc_f2r_ukb)

# eqtlg data from script 10-coloc ----------------------------------------------
f2r_eqtlg_coloc_snp = f2r_eqtlg_coloc[grep("^rs", f2r_eqtlg_coloc$SNP), ]
# need liftover for EnsDb.Hsapiens.v86 (hg38)
colnames(f2r_eqtlg_coloc_snp)[colnames(f2r_eqtlg_coloc_snp) == "chr.exposure"] <- "chrom"
colnames(f2r_eqtlg_coloc_snp)[colnames(f2r_eqtlg_coloc_snp) == "pos.exposure"] <- "pos"
colnames(f2r_eqtlg_coloc_snp)[colnames(f2r_eqtlg_coloc_snp) == "SNP"] <- "rsid"

f2r_eqtlg_coloc_snp <- perform_liftover(f2r_eqtlg_coloc_snp)

loc_f2r_eqtlg <- locus(data = f2r_eqtlg_coloc_snp,
                     gene = 'F2R',
                     flank = 1.5e5,
                     chrom = "chr_b38",
                     pos = "pos_b38",
                     p = "pval.exposure",
                     labs = "rsid_b37",
                     index_snp = "rs250735",
                     ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_eqtlg)
loc_f2r_eqtlg <- link_LD(loc_f2r_eqtlg, token = link_ld_token)
summary(loc_f2r_eqtlg)

# gtex data from script 10-coloc -----------------------------------------------
f2r_gtex_coloc_snp = f2r_gtexbld_format[grep("^rs", f2r_gtexbld_format$SNP), ] %>% subset(chr.exposure == 5)

loc_f2r_gtex <- locus(data = f2r_gtex_coloc_snp,
                       gene = 'F2R',
                       flank = 1.5e5,
                       chrom = "chr.exposure",
                       pos = "pos.exposure",
                       p = "pval.exposure",
                       index_snp = "rs250753",
                       labs = "SNP",
                       ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_gtex)
loc_f2r_gtex <- link_LD(loc_f2r_gtex, token = link_ld_token)
summary(loc_f2r_gtex)

# neph2 data -------------------------------------------------------------------
f2r_Glomneph2_hg3738
f2r_Tubeneph2_hg3738

loc_f2r_Glomneph2 <- locus(data = f2r_Glomneph2_hg3738,
                       gene = 'F2R',
                       flank = 1.5e5,
                       chrom = "chr_b38",
                       pos = "pos_b38",
                       p = "pval.outcome",
                       labs = "SNP",
                       index_snp = "rs6890835",
                       ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_Glomneph2)
loc_f2r_Glomneph2 <- link_LD(loc_f2r_Glomneph2, token = link_ld_token)
summary(loc_f2r_Glomneph2)

loc_f2r_Tubeneph2 <- locus(data = f2r_Tubeneph2_hg3738,
                           gene = 'F2R',
                           flank = 1.5e5,
                           chrom = "chr_b38",
                           pos = "pos_b38",
                           p = "pval.outcome",
                           labs = "SNP",
                           ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_Tubeneph2)
loc_f2r_Tubeneph2 <- link_LD(loc_f2r_Tubeneph2, token = link_ld_token)
summary(loc_f2r_Tubeneph2)

# plot.new()
# pdf(paste0(rdsf_personal,"results/Tubeneph2 locus.pdf"),width = 18, height = 5)
# print(locus_plot(loc_f2r_Tubeneph2, legend_pos = "topright", heights = 2))
# dev.off()

# combine plots in one figure for neph2 ----------------------------------------

pf <- quote({
  v <- loc_f2r_ukb$TX[loc_f2r_ukb$TX$gene_name == "F2R", c("start", "end")]
  abline(v = v, col = "orange")
})

pdf(paste0(rdsf_personal,"results/comb_coloc_neph2.pdf"), width = 10, height = 10)
# set up layered plot with 2 plots & a gene track; store old par() settings
oldpar <- set_layers(5)
scatter_plot(loc_f2r_ukb, xticks = FALSE, legend_pos = NULL, ylab = bquote(-log[10] * P ~ " (UKB-PPP)"), panel.first = pf)
scatter_plot(loc_f2r_eqtlg, xticks = FALSE, legend_pos = NULL, ylab = bquote(-log[10] * P ~ "(eQTLGen)"), panel.first = pf)
scatter_plot(loc_f2r_gtex, xticks = FALSE, legend_pos = NULL, ylab = bquote(-log[10] * P ~ "(GTEx blood)"), panel.first = pf)
scatter_plot(loc_f2r_Glomneph2, xticks = FALSE, legend_pos = NULL, ylab = bquote(-log[10] * P ~ "(NephQTL2 GLOM)"), panel.first = pf)
scatter_plot(loc_f2r_Tubeneph2, ylab = bquote(-log[10] * P ~ "(NephQTL2 TUBE)"), panel.first = pf)
genetracks(loc_f2r_ukb, highlight = "F2R", maxrows = 3, filter_gene_biotype = 'protein_coding')
           # gene_col = 'grey', exon_col = 'orange', exon_border = 'darkgrey')
par(oldpar)  # revert par() settings
dev.off()

