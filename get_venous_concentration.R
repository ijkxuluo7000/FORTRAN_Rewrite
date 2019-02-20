##############################################################
# Calculate Venous Concentration using Fick's Equation
# Author: Wanjun Gu
# Email: wag001@ucsd.edu
# Reference: Fick's Equation
# Fick's Equation: CVO2 = CaO2 - (VO2/(10*QT))
# Fick's Equation: CVCO2 = CaCO2 - (VCO2/(10*QT))
# University of California, San Diego
# UCSD School of Medicine
# Simonson Lab of Physiological Genomics of Altitude Adaptation
##############################################################

get_venous_concentration = function(vo2 = 1941.1,
                                    vco2 = 2484.1,
                                    qt = 10.6,
                                    cao2 = 24.62, 
                                    caco2 = 17.61){
  cvo2 = cao2 - (vo2/(10 * qt))
  cVco2 = caco2 + (vco2/(10 * qt))
  venous_concentration = c(cvo2, cVco2)
  names(venous_concentration) = c("CVO2", "CVCO2")
  return(venous_concentration)
}




