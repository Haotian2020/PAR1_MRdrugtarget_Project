# Script for main packages needed in this project

getwd()

# install.packages("remotes")
library(remotes)

# remotes::install_github("MRCIEU/TwoSampleMR",force = T)
library(TwoSampleMR)

# remotes::install_github("mrcieu/ieugwasr",force = T)
library(ieugwasr)

# install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)

library(forestplot)

# other package you might need

library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(tidyverse)
library(lubridate)
library(dplyr)
library(data.table)
library(reshape)
library(stringr)
library(meta)
library(devtools)