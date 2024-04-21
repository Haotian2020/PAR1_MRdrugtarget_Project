# Script for comparison between GWAS of ukb and finnGen around PAR1 gene region
# 38 
# 5: 76616126 - 76835770

devtools::install_github("myles-lewis/locuszoomr")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)
library(AnnotationHub)
library(LDlinkR)
library(locuszoomr)
library(qqman)
library(ggplot2)

ns_ukb = fread(paste0(par1_local,"data/ns_ukb.txt"))
head(ns_ukb)
ns_ukb <- ns_ukb[grep("^rs", ns_ukb$SNP), ]
ns_ukb = subset(ns_ukb,chr_b38 == 5) %>% subset(pos_b38 >= 76616126 & pos_b38 <= 76835770)
head(ns_ukb)

ns_fg = fread(paste0(par1_local,"data/ns_fg.txt"))
head(ns_fg)
ns_fg <- ns_fg[grep("^rs", ns_fg$rsids), ]
ns_fg = subset(ns_fg,`#chrom` == 5) %>% subset(pos >= 76616126 & pos <= 76835770)

meta = fread(paste0(par1_local,"data/meta_ns.txt"))
head(meta)

merged_snps <- unique(union(ns_ukb$SNP, ns_fg$rsids))
meta_local = subset(meta, MarkerName%in%merged_snps)
common_snps <- intersect(ns_ukb$SNP, ns_fg$rsids)
head(meta_local)
summary(meta_local)

meta_local_direct = subset(meta_local, Direction %in% c("++","--","+-","-+"))
meta_local_direct

meta_local_direct = left_join(meta_local_direct,
                              ns_ukb[ns_ukb$SNP%in% common_snps,c("SNP","chr_b38","pos_b38")],
                                                       by = c("MarkerName" = "SNP"))
summary(meta_local_direct)

loc_f2r_meta_h <- locus(data = meta_local_direct,
                     gene = 'F2R',
                     flank = 5e5,
                     labs = "MarkerName",
                     chrom = "chr_b38",
                     pos = "pos_b38",
                     p = "HetPVal",
                     ens_db = "EnsDb.Hsapiens.v86")

summary(loc_f2r_meta_h)
loc_f2r_meta_h <- link_LD(loc_f2r_meta_h, token = link_ld_token)
summary(loc_f2r_meta_h)
locus_plot(loc_f2r_meta_h)


p = manhattan(loc_f2r_meta_h$data, chr="chr_b38", bp="pos_b38", snp="MarkerName", p="HetPVal",
          suggestiveline = -log10(0.05)) 

ggplot_build(p)$plot%>% ggplot2::xlim(76616126, 76835770)

dat = loc_f2r_meta_h$data
dat$yn = ifelse(test = dat$HetPVal<0.05, yes = "<0.05", no = ">=0.05")

manhattan_plot <- ggplot(dat, aes(x = pos_b38, y = -log10(HetPVal))) +
  geom_point(aes(color = as.factor(yn)), alpha = 0.6) +
  scale_color_manual(values = c("blue", "red")) +
  theme_minimal() +
  ylim(0, 3) +
  xlim(76616126, 76835770) +
  coord_cartesian(ylim = c(0, 3))+ 
  annotate("text", x = dat$pos_b38[dat$HetPVal <= 0.05], 
                                          y = -log10(dat$HetPVal[dat$HetPVal <= 0.05]), 
                                          label = dat$MarkerName[dat$HetPVal <= 0.05], 
                                          color = "black", vjust = 2+runif(length(dat$MarkerName[dat$HetPVal <= 0.05]), -2, 2))
  

manhattan_plot

