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

tmp_unique_pos = tmp_input_exp_sig_keep[!duplicated(tmp_input_exp_sig_keep$pos),]

# using esemble database -------------------------------------------------------

grch38 = useEnsembl(biomart="snp",dataset = "hsapiens_snp")
x = listDatasets(grch38)
y = listAttributes(grch38)
z = listFilters(grch38)
attributes = c("refsnp_id","chr_name","chrom_start","chrom_end","allele")
filters = c("chr_name","start","end")

# request rsid from ensembl ----------------------------------------------------
# MAF from ensembl is global MAF, we need MAF for EUR --------------------------

rsid_store <- NULL

for (i in 1:length(tmp_unique_pos$chr)) {
  success <- FALSE
  while (!success) {
    finding_list <- lapply(tmp_unique_pos[i, c("chr", "pos", "pos")], as.numeric)
    print(i)
    print(finding_list)
    # Attempt to get data from the server
    rsid <- tryCatch(
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
    
    print(rsid)
    
    # Check if rsid is not NULL
    
    if (!is.null(rsid) && !identical(rsid, data.frame())) {

      rsid_store <- rbind(rsid_store, rsid)
      success <- TRUE  # Set success flag to exit the while loop
      rsid <- NULL
    }
  }
}

# QC step ----------------------------------------------------------------------
rsid_store_qc1 = subset(rsid_store,chrom_start == chrom_end)
rsid_store_qc1$pos = rsid_store_qc1$chrom_start
rsid_store_qc1 = rsid_store_qc1[,c("refsnp_id","chr_name","pos","allele")]
rsid_store_qc1 = subset(rsid_store_qc1,!is.na(allele))

unique_refsnp_rows <- rsid_store_qc1[!duplicated(rsid_store_qc1$refsnp_id), ]

merge_dat = left_join(x = tmp_unique_pos,unique_refsnp_rows,by = "pos")
nrow(merge_dat)

res = left_join(x = tmp_input_exp_sig_keep,unique_refsnp_rows,by = "pos")

# use SNPnexus -----------------------------------------------------------------

snpnexus_input = data.frame(merge_dat$refsnp_id)
snpnexus_input$type = "dbsnp"

snpnexus_input <- snpnexus_input[,c("type","merge_dat.refsnp_id")]

write.table(snpnexus_input, file = paste0(rdsf_personal,"data/par1/crosstissue_snpnexus.txt"),
            sep= '\t', row.names = F, col.names= F, quote = F)

# load maf from SNPnexus -------------------------------------------------------

maf = fread(paste0(paste0(rdsf_personal,"data/par1/crosstissue_maf.txt")))

maf = maf[,c("dbSNP","Chromosome","Position","Minor Allele","EUR Frequency")]

# merge data -------------------------------------------------------------------

final_merge = left_join(tmp_input_exp_sig_keep,maf,by = c("pos" = "Position"))

final_merge_missing = filter(final_merge,is.na(dbSNP))

final_merge_keep = filter(final_merge,!is.na(dbSNP))

# use spider to find missing rsid-----------------------------------------------

snps = final_merge_missing[,c("chr","pos")]
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

df <- data.frame(
  dbSNP = c("rs6884442","rs10514070","rs2227750","rs763692945","rs11954573","rs6880329","rs1743131620","rs458380","rs250733"),
  Chromosome = c(5,5,5,5,5,5,5,5,5),
  Position = c(76692609,76693219,76717639,77421814,76739242,76839299,77543474,76795672,76739176),
  `Minor Allele` = c("C","A","C","C","A","T",NA,"T","A"),
  `EUR Frequency` = c(0.22377,0.202996,0.14909,0,0.280649,0.4277,NA,0.15008,0.25688)
)

df = rbind(df,data.frame(maf))
df$EUR.Frequency = as.numeric(df$EUR.Frequency)

res = left_join(tmp_input_exp_sig_keep,df,by = c("pos" = "Position"))

res$eaf = ifelse(test = res$alt == res$Minor.Allele, yes = res$EUR.Frequency, no = 
                   ifelse(test = res$ref == res$Minor.Allele, yes = 1 - res$EUR.Frequency, no = NA))

res = filter(res,!is.na(eaf))

result <- res %>%
  group_by(dbSNP) %>%
  arrange(pval_nominal, desc(samplesize)) %>%
  filter(row_number() == 1) %>%
  ungroup()

cross_tissue_format = format_data(data.frame(result),
                                        type = "outcome",
                                        chr_col = "chr",
                                        pos_col = "pos",
                                        snp_col = "dbSNP",
                                        beta_col = "beta",
                                        se_col = "se",
                                        eaf_col = "eaf",
                                        effect_allele_col = "alt",
                                        other_allele_col = "ref",
                                        phenotype_col = "tissue",
                                        pval_col = "pval_nominal")
colnames(cross_tissue_format)[colnames(cross_tissue_format) == "outcome"] <- "tissue"
cross_tissue_format$outcome = "F2R cross tissues Gtex"

f2r_gtex_cross_str_exp = cross_tissue_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(. ,threshold = 5e-8, r2 = 0.001) %>% 
  mutate(exposure = "F2R Gtex cross str")

f2r_gtex_cross_wk_exp = cross_tissue_format %>% 
  dplyr::filter(chr.outcome == 5 & pos.outcome<=76835770 & pos.outcome>=76616126) %>% 
  ld_clump_local(. ,threshold = 5e-8, r2 = 0.1) %>% 
  mutate(exposure = "F2R Gtex cross wk")

