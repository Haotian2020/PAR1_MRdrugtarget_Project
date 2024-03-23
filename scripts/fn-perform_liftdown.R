# liftdown function
# 38 to 37

perform_liftdown <- function(data) {
  
  # Load data ------------------------------------------------------------------
  
  df <- data
  
  # Set variable names ---------------------------------------------------------
  
  df <- dplyr::rename(df,
                      "chr_b38" = "Chr",
                      "pos_b38" = "Pos",
                      "rsid_b38" = "MarkerName")
  
  # Restrict to relevant chromosomes -------------------------------------------
  
  df <- df[df$chr_b38 %in% paste0(1:22),] 
  df$chr_b38 <- paste0("chr",df$chr_b38)
  
  # Load chain file ------------------------------------------------------------
  
  ch <- rtracklayer::import.chain(paste0(rdsf_personal,"data/hg38ToHg19.over.chain"))
  # Make liftdown input --------------------------------------------------------
  chr_col <- "chr_b38"
  pos_col <- "pos_b38"
  
  liftdown_input <- GenomicRanges::GRanges(seqnames=df[[chr_col]], 
                                           ranges=IRanges::IRanges(start=df[[pos_col]], end=df[[pos_col]]), 
                                           LIFTDOWNCHRPOS=paste0(df[[chr_col]], ":", df[[pos_col]]))
  # Perform liftdown ----------------------------------------------------------- 
  liftdown_output <- rtracklayer::liftOver(liftdown_input, ch) %>% unlist()
  # Format liftdown output -----------------------------------------------------
  liftdown_output <- liftdown_output %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(LIFTDOWNCHRPOS=LIFTDOWNCHRPOS, LIFTDOWNCHR=seqnames, LIFTDOWNPOS=start)
  
  df$LIFTDOWNCHRPOS <- paste0(df[[chr_col]], ":", df[[pos_col]])
  
  # Combine liftdown input and output ------------------------------------------
  
  df <- merge(df, liftdown_output, by="LIFTDOWNCHRPOS")
  
  # Format data ----------------------------------------------------------------
  
  df <- dplyr::rename(df, "chr_b37" = "LIFTDOWNCHR", "pos_b37" = "LIFTDOWNPOS")
  df$LIFTDOWNCHRPOS <- NULL
  df$chr_b38 <- as.numeric(gsub("chr", "",as.character(df$chr_b38)))
  df$chr_b37 <- as.numeric(gsub("chr", "",as.character(df$chr_b37)))
  
  # Return variant manifest ----------------------------------------------------
  
  return(df)
}






