#Reference about the Kelman Curve
#Reference: http://www.ventworld.com/resources/oxydisso/dissoc.html

#Reference about the base excess
#Reference: https://media.lanecc.edu/users/driscolln/RT127/Softchalk/Acid_Base_Lesson/Acid_Base_Lesson10.html


#The input of this function should include at least 3 data points
#And at least 1 venous sample to maintain accuracy
#The input should belong to the same subject (Patient)
get_p50 = function(sat = c(42.0,92.7,90.0,90.8,26.1),
                   po2 = c(26,76,71,72,21),
                   temp = c(37,37,37,37,37),
                   ph = c(7.32,7.4,7.3,7.28,7.15),
                   pao2 = rep(mean(71,72,76),5), #What should be filling in this blank? ##
                   pco2 = c(45,3,33,31,64),
                   window_graphic = TRUE){
  
  
  #The virtual po2 will theoredically always be on the standard curve
  to_virtual_po2 = function(actual_po2,
                            temp = 37,
                            pH = 7.4,
                            PaO2 = 84,
                            PaCO2 = 40){
    
    virtual_po2 = actual_po2*10^(0.024*(37-temp)+0.40*(pH-7.40)+0.06*(log10(40)-log10(PaCO2)))
    return(virtual_po2)
    
  }
  
  ideal_curve = function(po2){
    
    #Parameters are empirical
    a1 = -8.5322289 * 10^3
    a2 = 2.1214010 * 10^3
    a3 = -6.7073989 * 10^1
    a4 = 9.3596087 * 10^5
    a5 = -3.1346258 * 10^4
    a6 = 2.3961674 * 10^3
    a7 = -6.7104406 * 10^1
    
    spO2 = 100*(a1*po2+a2*po2^2+a3*po2^3+po2^4)/(a4+a5*po2+a6*po2^2+a7*po2^3+po2^4)
    
    return(spO2)
    
  }
  
  #Print the standard Hemo-Oxy dessociation curve
  print_curve = function(){
    
    po2 = seq(1,100,0.01)
    spo2 = vector(length = length(seq(1,100,0.01)))
    
    n = 1
    for(trial_po2 in seq(1,100,0.01)){
      
      spo2[n] = ideal_curve(trial_po2)
      n = n + 1
      
    }
    
    plot(po2, spo2, 
         type = "l",
         main = "Standard Hemo-Oxy Dissociation Curve",
         xlab = "Partial Pressure of Oxygen",
         ylab = "Percent Hemoglobin Saturation of Oxygen")
    
  }
  
  #Convert actual po2 input to virtual po2
  virtual_po2 = vector(length = length(sat))
  for(i in 1:length(sat)){
    virtual_po2[i] = to_virtual_po2(po2[i],
                                    temp[i],
                                    ph[i],
                                    pao2[i],
                                    pco2[i])
  }
  
  #Using original data points instead
  virtual_po2 = po2
  
  #Add a shifting parameter lambda
  ssqlist = vector(length = length(seq(-30,30,0.001)))
  n = 1
  for(lambda in seq(-30,30,0.001)){
    
    ssq = vector()
    for(i in 1:length(sat)){
      ssq = sum(ssq, (sat[i] - ideal_curve(po2 = virtual_po2[i] + lambda))^2)
    }
    
    ssqlist[n] = ssq
    n = n + 1
    
  }
  
  #Determine the standard curve correction
  standard_curve_correction = -30 + which.min(ssqlist) * 0.001
  
  #Reconstruct the dissociation curve
  spo2_corrected = vector(length = length(seq(1,100,0.001)))
  n = 1
  for(trial_po2 in seq(1,100,0.001)){
    
    spo2_corrected[n] = ideal_curve(trial_po2 + standard_curve_correction)
    n = n + 1
    
  }
  
  if(window_graphic){
    
    #Open Graphic Devices
    x11()
    
    # #Print Standard Curve and place data points on it
    print_curve()
    points(virtual_po2, sat, pch = 16)
    
    #Plot the dissociation curve with the correction parameter
    plot(seq(1,100,0.001), spo2_corrected, 
         type = "l",
         main = "Fitted Hemo-Oxy Dissociation Curve",
         xlab = "Partial Pressure of Oxygen",
         ylab = "Percent Hemoglobin Saturation of Oxygen")
    #Place virtual points onto the plot
    points(virtual_po2, sat, cex = 1, pch = 16)
    abline(50, 0)
    legend(40, 48, legend = "50% Saturated")
    
  }
  
  #Create a lookup dictionary
  po2_lookup = seq(1,100,0.001)
  spo2_lookup = spo2_corrected
  corrected_lookup = cbind(po2_lookup, spo2_lookup)
  p50 = corrected_lookup[corrected_lookup[,2] == spo2_corrected[which.min(abs(50 - spo2_corrected))]][1]
  return(c(p50,standard_curve_correction))
  
  
}

#Data Input
sat = c(90.3,	85.2,	83.8,	84.4)
po2 = c(58,	50,	48,	54)
temp = rep(37,4)
ph = c(7.422,	7.313,	7.367,	7.305)
pao2 = c(58,	50,	48,	54)
pco2 = c(19.9, 23.4,	25.2,	22.3)


graphics.off()
get_p50(sat, po2, temp, ph, pao2, pco2)
get_p50()