# function to calculate F-stat

F_statistic = function(exposure_data){
  exposure_data$R = get_r_from_bsen(exposure_data$beta.exposure,exposure_data$se.exposure,exposure_data$samplesize.exposure)
  exposure_data$Rsq = exposure_data$R^2
  exposure_data$F_stat = (exposure_data$samplesize.exposure-2)*((exposure_data$Rsq)/(1-exposure_data$Rsq))
  exposure_data
  
  return(exposure_data)
}

