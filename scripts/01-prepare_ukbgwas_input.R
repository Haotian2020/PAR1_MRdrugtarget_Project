# Script for phenotype data input into UKB GWAS pipeline

`%notin%` <- Negate(`%in%`)
vars <- c("eid","31-0.0","21022-0.0", # ID, sex, age,
          paste0("41270-0.",seq(0,225)))

df <- fread(paste0(rdsf_path,"data.48733.csv"),
            sep = ",", 
            header = TRUE, 
            select = vars,
            data.table = FALSE)

link <- data.table::fread(paste0(rdsf_personal,"./data/linker_app15825_withexcl.csv"))
colnames(link) <- c("ieu","eid")

df_icd10 = df[,c('eid',paste0("41270-0.",seq(0,225)))]

codes <- reshape2::melt(df_icd10,
                        id.vars = c("eid"),
                        variable.name = "var_code",
                        value.name = "value_code")
codes$var_code <- gsub("-.*","",as.character(codes$var_code))

# FinnGen definition for NS ----------------------------------------------------
# include code start with N04
# exclude code N14_GLOMERULAR defined by FG, N00, N01, N02, N03, N05, N06, N07, N08

cases_rows <- codes[grep("^N04", codes$value_code), ]
# unique 410
exc_rows <- codes[grep("^(N00|N01|N02|N03|N05|N06|N07|N08)", codes$value_code), ]
# unique 3410

case_id = unique(cases_rows$eid)
control_exc_id = unique(exc_rows$eid)[unique(exc_rows$eid)%notin%unique(cases_rows$eid)]

df_icd10_N04 =subset(df,eid%notin%control_exc_id)
df_icd10_N04$ns = ifelse(df_icd10_N04$eid %in% case_id, yes = 1,no = 0)
df_icd10_N04_input = df_icd10_N04[,c("eid","ns")]

df_ns = merge(df_icd10_N04_input,link,by =c("eid"))

# 349 cases
# 459687 control

df_ns = df_ns[,c(3,3,2)]
colnames(df_ns) = c("FID","IID","ns")

# write the input file for ukb pipeline ----------------------------------------

write.table(df_ns,paste0(rdsf_personal,"data/df_ns.txt"),
            sep = " ",
            row.names = F, col.names = T,quote =F)