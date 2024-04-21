library(xml2)
library(XML)
library(dbparser)
library(canvasXpress)

citation("dbparser")

db <- parseDrugBank(paste0(rdsf_personal,"data/drugbank.xml"))

gi <- db$drugs$general_information %>% data.frame()

filtered_gi <- gi[grep("type 2 diabetes", gi$description, ignore.case = TRUE), ]

# load drugs data
drugs <- readRDS(system.file("drugs.RDS", package = "dbparser"))

# load drug groups data
drug_groups <- readRDS(system.file("drug_groups.RDS", package = "dbparser"))

# load drug targets actions data
drug_targets_actions <- readRDS(system.file("targets_actions.RDS", package = "dbparser"))

drugs
drug_groups
drug_targets_actions