##############################################################
# Get P50 Values from Blood-Gas Datasets
# Author: Wanjun Gu
# Email: wag001@ucsd.edu
# Reference: P50.FOR, Wagner PD, M.D
# University of California, San Diego
# UCSD School of Medicine
# Simonson Lab of Physiological Genomics of Altitude Adaptation
##############################################################

get_p50 = function(sat  = c(90.3,85.2,83.8,84.4,41.5),
                   po2  = c(58.0,50.0,48.0,54.0,24.0),
                   temp = c(37.0,37.0,37.0,37.0,37.0),
                   ph   = c(7.42,7.31,7.37,7.30,7.37),
                   pco2 = c(19.9,23.40,25.20,22.30,45.6),
                   lower_limit = 10,
                   upper_limit = 90,
                   graphics = FALSE){
  qualified_data_index = which(sat >= lower_limit & sat <= upper_limit)
  sat = sat[qualified_data_index]
  po2 = po2[qualified_data_index]
  temp = temp[qualified_data_index]
  ph = ph[qualified_data_index]
  pco2 = pco2[qualified_data_index]
  get_saturation = function(deviation_p50 = 0,
                            po2 = 26, 
                            temperature = 37,
                            ph = 7.4,
                            pco2 = 40){
    constant_list = c(-8532.229, 2121.401, -67.07399, 935960.9, -31346.26, 2396.167, -67.10441)
    pco2_correction_factor = 0.43429 * log(40.0/pco2)
    virtual_po2 = po2 * 10.0^(0.024 * (37.0 - temperature) + 0.4 * (ph - 7.4) + 0.06 * pco2_correction_factor)
    virtual_po2 = 26.8 * virtual_po2 / (26.8 + deviation_p50)
    if(virtual_po2 <= 10){
      saturation = 0.003683 * virtual_po2 + 0.000584 * virtual_po2^2
    }else{
      saturation = (virtual_po2*(virtual_po2*(virtual_po2*(virtual_po2+ constant_list[3])+constant_list[2])+constant_list[1]))/
        (virtual_po2*(virtual_po2*(virtual_po2*(virtual_po2+constant_list[7])+constant_list[6])+constant_list[5])+constant_list[4])
    }
    saturation = saturation * 100
    return(saturation)
  }
  ssq = vector()
  n = 1
  for(DP50 in seq(-16.8,23.2, 0.001)){
    sample_size = min(c(length(sat), length(po2), length(temp), length(ph), length(pco2)))
    virtual_saturation = vector()
    for(i in 1:sample_size){
      virtual_saturation[i] = get_saturation(deviation_p50 = DP50,
                                             po2 = po2[i],
                                             temperature = temp[i],
                                             ph = ph[i],
                                             pco2 = pco2[i])
    }
    ssq[n] = sum((virtual_saturation - sat)^2)
    n = n + 1
  }
  p50 = which.min(ssq)/1000 + 10
  if(graphics){
    P50_guess = seq(-16.8,23.2, 0.001) + 26.8
    plot(P50_guess, ssq,
         type = "l", 
         main = "SSQ Vs. Guessed P50 Values",
         sub = paste0("P50 = ", p50),
         xlab = "Guessed P50",
         ylab = "Corresponding SSQ")
    abline(v = p50)
  }
  return(p50)
}
