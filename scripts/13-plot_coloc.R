devtools::install_github("myles-lewis/locuszoomr")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(LDlinkR)
library(locuszoomr)

source("ld_token.R")

f2r_ukb_coloc_snp = f2r_ukb_coloc[grep("^rs", f2r_ukb_coloc$SNP), ]

loc_f2r_ukb <- locus(data = f2r_ukb_coloc_snp,
                     gene = 'F2R',
                     flank = 3e5,
                     labs = "SNP",
                     chrom = "chr.exposure",
                     pos = "pos.exposure",
                     p = "pval.exposure",
                     # yvar = "pval.outcome",
                     ens_db = "EnsDb.Hsapiens.v86")
summary(loc_f2r_ukb)
loc_f2r_ukb <- link_LD(loc_f2r_ukb, token = link_ld_token)
summary(loc_f2r_ukb)
p = locus_plot(loc_f2r_ukb)
pdf(paste0(rdsf_personal,"results/ukb locus.pdf"),width = 18, height = 5)
plot.new()
print(p)
dev.off()