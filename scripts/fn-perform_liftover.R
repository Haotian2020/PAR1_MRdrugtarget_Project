# liftover function
# 37 to 38

perform_liftover <- function(data) {
  
  # Load data ------------------------------------------------------------------
  
  df <- data
  
  # Set variable names ---------------------------------------------------------
  
  df <- dplyr::rename(df,
                      "chr_b37" = "chrom",
                      "pos_b37" = "pos",
                      "rsid_b37" = "rsid")
  
  # Restrict to relevant chromosomes -------------------------------------------
  
  df <- df[df$chr_b37 %in% paste0(1:22),] 
  df$chr_b37 <- paste0("chr",df$chr_b37)
  
  # Load chain file ------------------------------------------------------------
  
  ch <- rtracklayer::import.chain(paste0(rdsf_personal,"data/hg19ToHg38.over.chain"))
  # Make liftover input --------------------------------------------------------
  chr_col <- "chr_b37"
  pos_col <- "pos_b37"
  
  liftover_input <- GenomicRanges::GRanges(seqnames=df[[chr_col]], 
                                           ranges=IRanges::IRanges(start=df[[pos_col]], end=df[[pos_col]]), 
                                           LIFTOVERCHRPOS=paste0(df[[chr_col]], ":", df[[pos_col]]))
  # Perform liftover ----------------------------------------------------------- 
  liftover_output <- rtracklayer::liftOver(liftover_input, ch) %>% unlist()
  # Format liftover output -----------------------------------------------------
  liftover_output <- liftover_output %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(LIFTOVERCHRPOS=LIFTOVERCHRPOS, LIFTOVERCHR=seqnames, LIFTOVERPOS=start)
  
  df$LIFTOVERCHRPOS <- paste0(df[[chr_col]], ":", df[[pos_col]])
  
  # Combine liftover input and output ------------------------------------------
  
  df <- merge(df, liftover_output, by="LIFTOVERCHRPOS")
  
  # Format data ----------------------------------------------------------------
  
  df <- dplyr::rename(df, "chr_b38" = "LIFTOVERCHR", "pos_b38" = "LIFTOVERPOS")
  df$LIFTOVERCHRPOS <- NULL
  df$chr_b37 <- as.numeric(gsub("chr", "",as.character(df$chr_b37)))
  df$chr_b38 <- as.numeric(gsub("chr", "",as.character(df$chr_b38)))
  
  # Return UKB variant manifest ------------------------------------------------
  
  return(df)
}
