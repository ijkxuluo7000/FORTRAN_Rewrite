########################Variables
inspo2 = c(21,21,21,21)
po2 = c(50,53,54,55)
pco2 = c(25.6,26.3,27.2,24.6)
ph = c(7.323,7.309,7.256,7.185,7.356)
sat = c(86.4,94.4,93.0,85.0)
hb = c(19.2,20.0,20.6,20.5)
hct = c(71,71,72,73)
temp = c(36.8,36.7,36.7,36.7)
dat = as.data.frame(cbind(inspo2,po2,pco2,ph,sat,hb,hct,temp))
###############################
#Function Saturation
saturation = function(PO2,TEMP,PH,PCO2,P50){
  a1 = -8532.229;a2 = 2121.401
  a3 = -67.07399;a4 = 935960.9
  a5 = -31346.26;a6 = 2396.167
  a7 = -67.10441
  b = 0.43429 * log(40.0/PCO2,base = 10)
  xx = PO2 * 10^(0.024 * (37 - TEMP) + 0.4 * (PH - 7.4) + 0.06 * b)
  x = 26.8 * (xx/P50)
  sat = (x*(x*(x*(x+a3)+a2)+a1))/(x*(x*(x*(x+a7)+a6)+a5)+a4)
  return(100 * sat)
}

lower_bond = 10
upper_bond = 60
digit = 1
vp50 = seq(lower_bond,upper_bond,10^(-digit))

size = length(dat[,1])
satura_mat = matrix(ncol = size,nrow = length(vp50))

for(i in 1:size){
  for(j in 1:length(vp50)){
    satura_mat[j,i] = saturation(PO2 = dat$po2[i],
                               TEMP = dat$temp[i],
                               PH = dat$ph[i],
                               PCO2 = dat$pco2[i],
                               P50 = vp50[j])
  }
}

ssqvec = vector(length = length(vp50))
for(j in 1:length(vp50)){
  ssq = vector(length = size)
  for(i in 1:size){
    ssq[i] = (satura_mat[,i] - dat$sat[i])^2
  }
  ssq = sum(ssq)
  ssqvec[j] = ssq
}
