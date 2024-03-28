# extract instruments from gwas data -------------------------------------------
# F2R gene location
# flanking region 100kb

# 38 build Chromosome 5: 76,716,126-76,735,770 
# 38build location to be extracted: 5: 76616126 - 76835770

# 37 build Chromosome 5: 76,011,951-76,031,595
# 37 build location to be extracted: 5: 75911951 - 76131595

source("fn-ld_clump_local")

# gtex blood dataset -----------------------------------------------------------
# 38 build

f2r_gtexbld_gwas = fread(paste0(rdsf_personal,"data/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz"))
f2r_gtexbld_gwas = f2r_gtexbld_gwas[grepl("ENSG00000181104", f2r_gtexbld_gwas$gene_id), ]
f2r_gtexbld_gwas = subset(f2r_gtexbld_gwas,pval_nominal<=5e-8)

# manually add SNP -------------------------------------------------------------

f2r_gtexbld_gwas$SNP = c("rs250737","rs250735","rs250734","rs250753")

# check with reference panel that each alt is the minor allele -----------------

split_cols <- strsplit(f2r_gtexbld_gwas$variant_id, "_")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos", "ref","alt","build")
f2r_gtexbld_gwas <- cbind(f2r_gtexbld_gwas, split_df)
f2r_gtexbld_gwas$chr <- 5

f2r_gtexbld_gwas_format = format_data(
  dat = data.frame(f2r_gtexbld_gwas),
  type = "outcome",
  snp_col = "SNP",
  effect_allele_col = "alt",
  other_allele_col = "ref",
  eaf_col = "maf",
  beta_col = "slope",
  se_col = "slope_se",
  pval_col = "pval_nominal",
  chr_col = "chr",
  pos_col = "pos") %>% mutate(samplesize.outcome = 755)

f2r_gtex_str_exp = f2r_gtexbld_gwas_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(. ,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R Gtex str")

# weakly correlated SNP is same as the strongly independent one

# no siginificant signal in gtex v8 kidney tissue ------------------------------

write.table(f2r_gtex_str_exp, file = paste0(rdsf_personal,"data/par1/f2r_gtex_str_exp.csv"),
            sep= ',', row.names = F,col.names= T)

# cross-tissue gtex instruments ------------------------------------------------

drug_targets<- NULL
drug_targets$gene_id = "ENSG00000181104"
drug_targets$gene = "F2R"
gtex_path <- "data/GTEx_Analysis_v8_eQTL/"
gtex_data <- ".v8.signif_variant_gene_pairs.txt.gz"
# gtex_data <- ".v8.egenes.txt.gz"
tissues <- list.files(path = paste0(rdsf_personal,gtex_path),pattern = paste0("*",gtex_data))

input_exp <- data.frame(gene = rep(unique(drug_targets$gene),each = length(tissues)),
                        tissue = rep(tissues,times = length(unique(drug_targets$gene))),
                        stringsAsFactors = FALSE)

tmp_input_exp <- NULL

for (i in tissues) {
  tmp <- vroom(paste0(rdsf_personal,gtex_path,i))
  subset_tmp <- tmp[grepl("ENSG00000181104", tmp$gene_id), ]
  subset_tmp$tissue <- gsub(gtex_data,"",i)
  tmp_input_exp <- rbind(tmp_input_exp,subset_tmp)
}

tmp_input_exp <- data.frame(tmp_input_exp)
tmp_input_exp <- tmp_input_exp[,c("variant_id","gene_id","tissue","slope","slope_se","maf",
                                  "pval_nominal","pval_nominal_threshold")] %>% data.frame()

tmp_input_exp$samplesize = NA

# Load and format sample size file from gtex portal ----------------------------

gtex_ss = fread(paste0(rdsf_personal,"data/GTEx Portal sample size.csv"))
gtex_ss$Tissue <- gsub(" - ", "_", gtex_ss$Tissue)
gtex_ss$Tissue <- gsub(" ", "_", gtex_ss$Tissue)
gtex_ss$Tissue <- gsub("\\(", "", gtex_ss$Tissue)
gtex_ss$Tissue <- gsub("\\)", "", gtex_ss$Tissue)

# add sample size to tmp_input_exp ---------------------------------------------

for(row in 1:nrow(tmp_input_exp)){
  tmp_tissue <- tmp_input_exp[row,"tissue"]
  if(any(gtex_ss$Tissue == tmp_tissue)){
    tmp_input_exp[row,"samplesize"] <- gtex_ss[gtex_ss$Tissue == tmp_tissue,"# RNASeq Samples"]
  } else {
    tmp_input_exp[row,"samplesize"] <- NA
  }
}

tmp_input_exp_sig = filter(tmp_input_exp, pval_nominal<5e-8)

split_cols <- strsplit(tmp_input_exp_sig$variant_id, "_")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos", "ref","alt","build")
tmp_input_exp_sig <- cbind(tmp_input_exp_sig, split_df)

# exclude SNPs with ref or alt more than 1 characters --------------------------

tmp_input_exp_sig_keep = subset(tmp_input_exp_sig,nchar(tmp_input_exp_sig$ref)==1 & 
                                  nchar(tmp_input_exp_sig$alt)==1)

# format dataset ---------------------------------------------------------------

colnames(tmp_input_exp_sig_keep)[colnames(tmp_input_exp_sig_keep) == "slope"] <- "beta"
colnames(tmp_input_exp_sig_keep)[colnames(tmp_input_exp_sig_keep) == "slope_se"] <- "se"
tmp_input_exp_sig_keep$pos = as.numeric(tmp_input_exp_sig_keep$pos)
tmp_input_exp_sig_keep$chr = 5
tmp_input_exp_sig_keep$chr_pos = paste0(tmp_input_exp_sig_keep$chr,":",tmp_input_exp_sig_keep$pos)

# use SNPnexus -----------------------------------------------------------------

snpnexus_input = tmp_input_exp_sig_keep[,c("chr","pos","pos","ref","alt")]
colnames(snpnexus_input) = c("Chromosome","Start","End","ref","alt")
snpnexus_input <- snpnexus_input[rep(row.names(snpnexus_input), each = 2), ]

snpnexus_input$strand <- 1
snpnexus_input$strand[seq(2, nrow(snpnexus_input), by = 2)] <- -1
rownames(snpnexus_input) <- NULL
snpnexus_input$type = "Chromosome"
snpnexus_input = snpnexus_input[,c("type","Chromosome","Start","End","ref","alt","strand")]

write.table(snpnexus_input, file = paste0(rdsf_personal,"data/par1/crosstissue_snpnexus.txt"),
            sep= '\t', row.names = F, col.names= T, quote = F)



# using esemble database -------------------------------------------------------

snp = useEnsembl(biomart="snp")
grch38 = useEnsembl(biomart="snp",dataset = "hsapiens_snp")
x = listDatasets(grch38)
y = listAttributes(grch38)
z = listFilters(grch38)
attributes = c("refsnp_id","chr_name","chrom_start","chrom_end","allele","minor_allele","minor_allele_freq")
filters = c("chr_name","start","end")

# request maf from esemble -----------------------------------------------------
maf_store <- NULL

for (i in 1:length(tmp_input_exp_sig_keep$chr)) {
  success <- FALSE
  while (!success) {
    finding_list <- lapply(tmp_input_exp_sig_keep[i, c("chr", "pos", "pos")], as.numeric)
    print(i)
    print(finding_list)
    # Attempt to get data from the server
    rsid_maf <- tryCatch(
      expr = {
        getBM(attributes = attributes,
              filters = filters,
              values = finding_list, mart = grch38)
      },
      error = function(e) {
        # Print an error message if needed
        
        print(paste("Error in iteration", i, ":", conditionMessage(e)))
        
        # Return a placeholder value (e.g., an empty data frame)
        
        return(data.frame())
      }
    )
    
    print(rsid_maf)
    
    # Check if rsid_maf is not NULL
    if (!is.null(rsid_maf) && !identical(rsid_maf, data.frame())) {
      print("find maf is not NULL")
      maf_store <- rbind(maf_store, rsid_maf)
      success <- TRUE  # Set success flag to exit the while loop
      rsid_maf <- NULL
    }
  }
}

# QC step ----------------------------------------------------------------------
maf_store_qc1 = subset(maf_store,chrom_start == chrom_end)
maf_store_qc1$pos = maf_store_qc1$chrom_start
maf_store_qc1 = maf_store_qc1[,c("refsnp_id","chr_name","pos","allele","minor_allele","minor_allele_freq")]
maf_store_qc1 = subset(maf_store_qc1,!is.na(minor_allele)&!is.na(minor_allele_freq))
# 1459

unique_refsnp_rows <- maf_store_qc1[!duplicated(maf_store_qc1$refsnp_id), ]
unique_refsnp_rows$minor_allele[unique_refsnp_rows$minor_allele == TRUE] <- "T"

# merge tmp_input_exp_sig_keep with maf_store_qc1 by pos colume ----------------
merge_dat = left_join(x = tmp_input_exp_sig_keep,unique_refsnp_rows,by = "pos")
nrow(merge_dat)
# 1490
table(is.na(merge_dat$refsnp_id))
# 31
table(is.na(merge_dat$minor_allele_freq))
# 31

# check if maf = eaf -----------------------------------------------------------

table(merge_dat$alt == merge_dat$minor_allele)
# FALSE  TRUE
# 699   760
merge_dat$eaf = ifelse(merge_dat$alt == merge_dat$minor_allele,merge_dat$maf,1-merge_dat$maf)
merge_dat_keep = subset(merge_dat,!is.na(refsnp_id))
merge_dat_miss = subset(merge_dat,is.na(refsnp_id))
# check merge_dat_miss one by one
merge_dat_miss$pos
# duplicated 

snps = merge_dat_miss[,c("chr","pos")]
start_time <- Sys.time()
info = data.frame()
for(i in 1:nrow(snps)){
  # snps is in chr:pos style
  print(snps[i,])
  url = paste0("https://www.ncbi.nlm.nih.gov/snp/?term=",snps[i,"chr"],"%3A",snps[i,"pos"])
  web = read_html(url)
  # context is the information for each title
  title = web %>% html_nodes("div dt")%>% html_text()
  context = web %>% html_nodes("div dd")%>% html_text()
  if(length(context)>0){
    index_style = which("Variant type:" ==title)
    index_alleles = which("Alleles:" ==title)
    merge_info = web %>% html_nodes("div div div span")%>% html_text()
    merge_info = merge_info[startsWith(merge_info,"rs")]
    rsid_infor = web %>% html_nodes("dd span") %>% html_attrs()
    rsid = c()
    dat = data.frame()
    print(paste0("For SNP ",snps$chr[i],":",snps$pos[i],", How many results we can get from dbSNP: ",length(context)/length(unique(title))))
    for(j in 1:length(rsid_infor)){
      if(length(rsid_infor[[j]]>0)){
        print(j)
        print("something in list")
        for(k in 1:length(rsid_infor[[j]]))
          print(k)
        extracted_value <- str_extract(rsid_infor[[j]][k], "^rs[^:]+")
        if (!is.na(extracted_value)) {
          print("Extracting rsid")
          rsid <- c(rsid, extracted_value)
        }}}
    print(rsid)
    snp_style = context[index_style]
    allele_type = gsub("\\n| \\[Show Flanks\\]| ", "",context[index_alleles])
    dat = data.frame(rsid) %>% na.omit()
    dat$Variant_type = snp_style
    dat$Alleles = allele_type
    dat$merge_info = merge_info
    GRCh3738 = context[grepl("GRCh38",context)]
    clean_GRCh3738 <- gsub("\\(GRCh38\\)|\\(GRCh37\\)", "", GRCh3738)
    cleaned_chrpos <- strsplit(clean_GRCh3738, "\n")
    dat$GRCh38 <- vector("character", length(cleaned_chrpos))
    dat$GRCh37 <- vector("character", length(cleaned_chrpos))
    for (x in 1:length(cleaned_chrpos)) {
      dat$GRCh38[x] <- cleaned_chrpos[[x]][1]
      dat$GRCh37[x] <- cleaned_chrpos[[x]][2]
    }
    info = rbind(info,dat)
  }else{miss_rsid = rbind(miss_rsid,snps[i,])}
}

end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)
info
nrow(info)

filtered_data <- info %>% filter(!str_detect(merge_info, "has merged into"))
filtered_data <- subset(filtered_data,!duplicated(rsid)) %>% subset(Variant_type=="SNV")
# 9

split_cols <- strsplit(filtered_data$GRCh38, ":")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("chr", "pos")
filtered_data = cbind(filtered_data,split_df) %>% .[c("rsid","Alleles","chr","pos")]
filtered_data$pos = as.numeric(filtered_data$pos)

merge_dat_miss <- left_join(merge_dat_miss,filtered_data, by = "pos")
merge_dat_miss$refsnp_id = merge_dat_miss$rsid
merge_dat_miss$allele = merge_dat_miss$Alleles
split_cols <- strsplit(merge_dat_miss$allele, ">")
split_df <- as.data.frame(do.call(rbind, split_cols))
colnames(split_df) <- c("major_allele", "minor_allele")
merge_dat_miss = cbind(merge_dat_miss,split_df)
table(merge_dat_miss$major_allele == merge_dat_miss$ref)
# TRUE 31
merge_dat_miss$eaf = ifelse(merge_dat_miss$major_allele == merge_dat_miss$ref,merge_dat_miss$maf,1-merge_dat_miss$maf)
merge_dat_miss$minor_allele = merge_dat_miss$minor_allele.1
merge_dat_miss$chr = merge_dat_miss$chr.x
final_merge = rbind(merge_dat_keep,merge_dat_miss[,c(colnames(merge_dat_keep))])
nrow(final_merge)
# 1490
# if one SNP shows in multiple tissue, keep the one with lowest p
# final_merge <- final_merge[order(final_merge$pval_nominal), ]
final_merge <- final_merge[order(final_merge$pval_nominal,
                                 final_merge$refsnp_id), ]

final_merge_group <- final_merge %>%
  group_by(refsnp_id) %>%
  filter(pval_nominal == min(pval_nominal)) %>% data.frame()
# 247

write.table(final_merge_group$refsnp_id, file = paste0(rdsf_personal,"data/par1/f2r_gtexcrosstissue_rsid.txt"),
            sep= ' ', row.names = F,col.names= F,quote = F)

# use ensembl VEP to get Eur AF ------------------------------------------------

cross_tissue_rsid_maf = fread(paste0(rdsf_personal,"data/par1/crosstissue_rsid_eurmaf.txt"))
cross_tissue_rsid_maf = cross_tissue_rsid_maf[,c("#Uploaded_variation","Allele","EUR_AF")]
colnames(cross_tissue_rsid_maf) = c("refsnp_id","Allele","EUR_AF")
cross_tissue_rsid_maf = subset(cross_tissue_rsid_maf,EUR_AF != "-")
cross_tissue_rsid_maf = cross_tissue_rsid_maf[!duplicated(cross_tissue_rsid_maf$refsnp_id),]
final_merge_group = left_join(final_merge_group,cross_tissue_rsid_maf,by = "refsnp_id")

table(is.na(final_merge_group$Allele))
# 240 + 7
table(final_merge_group$alt == final_merge_group$Allele)
#   2 + 238
# 2 FALSE contain AF = 0
table(is.na(final_merge_group$Allele))
# 7 SNPs no maf information
# need to find 9 SNPs manually
final_merge_group_keep = final_merge_group[!is.na(final_merge_group$Allele)&final_merge_group$EUR_AF != 0,]
# all allele = alt
final_merge_group_keep$eaf = ifelse(final_merge_group_keep$EUR_AF<0.5,final_merge_group_keep$maf,1-final_merge_group_keep$maf)

final_merge_group_missing = final_merge_group[is.na(final_merge_group$Allele)|final_merge_group$EUR_AF == 0,]
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs250733",c("Allele","EUR_AF")] = c("G","0.73734")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs458380",c("Allele","EUR_AF")] = c("T","0.15008")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs6880329",c("Allele","EUR_AF")] = c("T","0.4277")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs11954573",c("Allele","EUR_AF")] = c("A","0.280649")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs2227750",c("Allele","EUR_AF")] = c("C","0.14909")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs463188",c("Allele","EUR_AF")] = c("T","0.7030")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs6884442",c("Allele","EUR_AF")] = c("C","0.22377")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs10514070",c("Allele","EUR_AF")] = c("A","0.202996")
final_merge_group_missing[final_merge_group_missing$refsnp_id == "rs153308",c("Allele","EUR_AF")] = c("T","0.2523")
final_merge_group_missing$eaf = ifelse(final_merge_group_missing$EUR_AF<0.5,final_merge_group_missing$maf,1-final_merge_group_missing$maf)

final_cross_tissue = rbind(final_merge_group_missing,final_merge_group_keep)

final_cross_tissue_format = format_data(final_cross_tissue,
                                        type = "exposure",
                                        chr_col = "chr",
                                        pos_col = "pos",
                                        snp_col = "refsnp_id",
                                        beta_col = "beta",
                                        se_col = "se",
                                        eaf_col = "eaf",
                                        effect_allele_col = "alt",
                                        other_allele_col = "ref",
                                        phenotype_col = "tissue",
                                        pval_col = "pval_nominal")
colnames(final_cross_tissue_format)[colnames(final_cross_tissue_format) == "exposure"] <- "tissue"
final_cross_tissue_format$exposure = "F2R cross tissues Gtex"

