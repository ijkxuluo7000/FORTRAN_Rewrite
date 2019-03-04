##############################################################
# Calculate Pressure from Concentration
# Author: Wanjun Gu
# Email: wag001@ucsd.edu
# Reference: PFROMC.FOR, Wagner PD, M.D
# University of California, San Diego
# UCSD School of Medicine
# Simonson Lab of Physiological Genomics of Altitude Adaptation
##############################################################

hb = 20.5
hct = 61.4
temp = 37
aph = 7.16
bph = 7.04
p50 = 24.6
pb = 457.0
pio2 = 85.8

sameo2 = function(ppo, ppco, arto2c, co2ct2, pio2, o2cnt2){
  blood = function(pco2, hb, hcrit, temp, po2, p50, aph, bph){
    
    get_ph = function(pco2, y, aph, bph, apco2 = 30, bpco2 = 60){
      
      if(pco2 < 0.001){
        pco2 = 0.001
      }
      if(aph < 1){
        ph = 7.59 + y - 0.2741 * log(pco2/20.0)
      }else{
        ph = bph + y + (aph - bph) * log(pco2/bpco2)/log(apco2/bpco2)
      }
      return(ph)
      
    }
    satura = function(po2, pco2, phe, p50){
      
      a1 = -8532.229
      a2 = 2121.401
      a3 = -67.07399
      a4 = 935960.9
      a5 = -31346.26
      a6 = 2396.167
      a7 = -67.10441
      b = 0.43429 * log(40.0/pco2)
      x = po2 * 10.0^(0.024 * (37.0 - temp) + 0.4 * (phe - 7.4) + 0.06 * b)
      x = 26.8 * x / p50
      if(x <= 10){
        sat = 0.003683 * x + 0.000584 * x^2
      }else{
        sat=(x*(x*(x*(x+a3)+a2)+a1))/(x*(x*(x*(x+a7)+a6)+a5)+a4)
      }
      satura = 100.0 * sat
      return(satura)
      
    }
    co2con = function(phe, satn, temp, hcrit, pco2){
      
      p = 7.4 - phe
      pk = 6.086 + 0.042 * p + (38.0 - temp) * (0.00472 + 0.00139 * p)
      sol = 0.0307 + 0.00057 * (37.0 - temp) + 0.00002 * (37.0 - temp) * (37.0 - temp)
      dox = 0.59 + 0.2913 * p-0.0844 * p^2
      dr = 0.664 + 0.2275 * p - 0.0938 * p^2
      ddd = dox + (dr-dox) * (1 - satn/100.0)
      cp = sol * pco2 * (1.0 + 10.0^(phe - pk))
      ccc = ddd * cp
      co2con = (hcrit * ccc * 0.01 + (1.0 - hcrit * 0.01) * cp) * 2.22
      return(co2con)
      
    }
    
    ph1 = get_ph(pco2 = pco2, y = 0, aph = aph, bph = bph)
    y = 0.003 * hb * (1 - satura(po2 = po2, pco2 = pco2, phe = ph1, p50 = p50)/100.0)
    ph2 = get_ph(pco2, y = y, aph = aph, bph = bph)
    satrn = satura(po2 = po2, pco2 = pco2, phe = ph2, p50 = p50)
    o2c = 0.0139 * hb * satrn + 0.003 * po2
    co2c = co2con(pco2 = pco2, phe = ph2, satn = satrn, temp = temp, hcrit = hcrit)
    concentration = c(o2c, co2c)
    names(concentration) = c("o2c", "co2c")
    return(concentration)
    
  }
  
  n = 1
  
  e = 0.0
  f = pio2 + 10.0
  g = (e + f)/2.0
  while(abs(o2cnt2-arto2c) > 0.0001){

    #co2ct2 = blood(po2 = g, pco2 = ppco, hb = hb, p50 = p50, aph = aph, bph = bph)[2]
    if(o2cnt2 < arto2c){
      e = g
    }else if(o2cnt2 > arto2c){
      f = g
    }
    
    o2cnt2 = blood(po2 = g, pco2 = ppco, hb = hb,
                   hcrit = hct, p50 = p50, aph = aph,
                   bph = bph, temp = temp)[1]
    
    n = n + 1
    if(n > 10000){
      print("Infinite Loop")
      break
    }
    
  }
  ppo = g
  return(ppo)
}


