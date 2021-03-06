      PROGRAM CALCON
C
C     PROGRAM IN FORTRAN TO COMPUTE BLOOD GAS CONTENTS FROM PARTIAL
C     PRESSURES USING THE KELMAN DISSOCIATION CURVES
C
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      CHARACTER IOK
      WRITE(*,5)
  5   FORMAT(' OUTPUT TO DISK (0) OR PRINTER (1) ?')
      READ(*,*)IPR
      IF(IPR.EQ.0) OPEN(3,FILE=' ',STATUS='NEW')
      IF(IPR.EQ.1) OPEN(3,FILE='PRN')
      APCO2=30.0
      BPCO2=60.0
      NSETS=0
  10  CONTINUE
C
C     LOOP THRU AS MANY TIMES AS YOU LIKE
C
      NSETS=NSETS+1
C
C     ENTER BLOOD GAS DATA: PO2, PCO2, PH  (AT BODY TEMPERATURE)
C
  20  WRITE(*,30)NSETS
  30  FORMAT(' FOR SET NUMBER',I3,/,
     1' ENTER PO2,PCO2 & PH OF SAMPLE, AT BODY TEMP,FORMAT FREE')
      READ(*,*)PO2,PCO2,PHE
C
C     CHECK DATA ENTRY; IF WRONG RE-ENTER
C
      WRITE(*,40)PO2,PCO2,PHE
  40  FORMAT(3F10.2,10X,'OK?  (Y/N)')
      READ(*,50)IOK
  50  FORMAT(A1)
      IF(IOK.EQ.'N'.OR.IOK.EQ.'n') GO TO 20
  60  WRITE(*,70)
  70  FORMAT(' ENTER HB,% HCRIT, TEMP & STANDARD P50, FORMAT FREE')
C
C     MORE DATA ENTRY: HEMOGLOBIN, HEMATOCRIT, TEMPERATURE  & P50
C
C     P50 IS AT STANDARD CONDITIONS,IE., AT PH=7.40 & TEMP=37.0
C     ITS NORMAL VALUE IS 26.8 TORR WHICH MAY BE ASSUMED IF THE REAL
C     VALUE IS NOT MEASURED
C
      READ(*,*)HB,HCRIT,TEMP,P50
C
C     CHECK DATA ENTRY AS BEFORE, RE-ENTER IF NECESSARY
C
      WRITE(*,80)HB,HCRIT,TEMP,P50
  80  FORMAT(4F10.1,10X,'OK?  (Y/N)')
      READ(*,50)IOK
      IF(IOK.EQ.'N'.OR.IOK.EQ.'n') GO TO 60
C
C     SUBROUTINE BUFLIN COMPUTES IN EFFECT TWO POINTS ON THE PH-PCO2
C     BLOOD BUFFER LINE SO THAT ACID/BASE ALTERATIONS AFFECTING THE
C     PH-PCO2 RELATIONSHIPS CAN BE TAKEN INTO ACCOUNT WHEN COMPUTING
C     O2/CO2 CONCENTRATIONS FROM PARTIAL PRESSURES
C
      CALL BUFLIN(PO2,PCO2,PHE)
C
C     SUBROUTINE BLOOD IS THE ROUTINE TO DO THE REQUIRED CALCULATIONS
C     OF BLOOD GAS CONCENTRATIONS
C
      CALL BLOOD(PO2,PCO2,O2CONT,CO2CNT,SATRN)
C
C     WRITE OUT THE INPUT DATA AND THE RESULTS AS WELL
C
      WRITE(*,90)NSETS,PO2,PCO2,PHE,SATRN,O2CONT,CO2CNT
      WRITE(3,90)NSETS,PO2,PCO2,PHE,SATRN,O2CONT,CO2CNT
  90  FORMAT(' ******* THIS IS SET NUMBER ',I3,' (KELMAN) *******'//,
     1' PO2 =',F7.1,4X,'TORR',/,
     2' PCO2=',F7.1,4X,'TORR',/,
     3' PH  =',F8.2,/,
     4' % O2 SATURATION  =',F8.1,/,
     5' O2 CONCENTRATION =',F10.3,4X,'ML PER 100 ML BLOOD',/,
     6' CO2 CONCENTRATION=',F10.3,4X,'ML PER 100 ML BLOOD',/)
      WRITE(*,100)HB,HCRIT,TEMP,P50,APCO2,APH,BPCO2,BPH
      WRITE(3,100)HB,HCRIT,TEMP,P50,APCO2,APH,BPCO2,BPH
 100  FORMAT(' HEMOGLOBIN CONC=',F7.1,4X,'GM PER DECILITER',/,
     1' % HEMATOCRIT   =',F7.1,/,
     2' TEMPERATURE    =',F7.1,4X,'DEGREES CENTIGRADE',/,
     3' STANDARD P50   =',F7.1,4X,'TORR',/,
     4' PH FOR A PCO2 OF',F6.1,2X,'TORR   =',F7.2,/,
     5' PH FOR A PCO2 OF',F6.1,2X,'TORR   =',F7.2,////)
C
C     CYCLE BACK TO THE TOP FOR MORE DATA SETS IF DESIRED
C
      WRITE(*,110)
 110  FORMAT(' WANT TO RUN MORE DATA?  (Y/N)')
      READ(*,50)IOK
      IF(IOK.EQ.'Y'.OR.IOK.EQ.'y') GO TO 10
      END
      SUBROUTINE BLOOD(PO2,PCO2,O2C,CO2C,SATRN)
C
C     CODE TO CALL THE DISSOCIATION CURVE SUBROUTINES
C 
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      PH1=PH(PCO2,0.0)
      Y=0.003*HB*(1.0-SATURA(PO2,PCO2,PH1)/100.0)
      PH2=PH(PCO2,Y)
      SATRN=SATURA(PO2,PCO2,PH2)
      O2C=0.0139*HB*SATRN + 0.003*PO2
      CO2C=CO2CON(PCO2,PH2,SATRN)
      RETURN
      END
      FUNCTION PH(PCO2,Y)
C
C     KELMAN'S SUBROUTINE FOR COMPUTING PH
C
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      IF(PCO2.LT.0.001) PCO2=0.001
      IF(APH-1.0) 10,10,20
  10  PH=7.59 + Y - 0.2741*ALOG(PCO2/20.0)
      GO TO 30
  20  PH=BPH+Y+(APH-BPH)*ALOG(PCO2/BPCO2)/ALOG(APCO2/BPCO2)
  30  RETURN
      END
      FUNCTION CO2CON(PCO2,PHE,SATN)
C
C     KELMAN'S SUBROUTINE FOR THE CO2 DISSOCIATION CURVE
C
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      P=7.4-PHE
      PK=6.086+0.042*P+(38.0-TEMP)*(0.00472+0.00139*P)
      SOL=0.0307 + 0.00057*(37.0-TEMP) + 0.00002*(37.0-TEMP)*(37.0-TEMP)
      DOX=0.59+0.2913*P-0.0844*P*P
      DR=0.664+0.2275*P-0.0938*P*P
      DDD=DOX+(DR-DOX)*(1.-SATN/100.0)
      CP=SOL*PCO2*(1.0+10.0**(PHE-PK))
      CCC=DDD*CP
      CO2CON=(HCRIT*CCC*0.01 + (1.0-HCRIT*0.01)*CP)*2.22
      RETURN
      END
      FUNCTION SATURA(PO2,PCO2,PHE)
C
C     KELMAN'S SUBROUTINE FOR THE O2 DISSOCIATION CURVE
C     MODIFIED BY INCLUSION OF POSSIBILITY OF VARYING THE P50
C
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      A1=-8532.229
      A2=2121.401
      A3=-67.07399
      A4=935960.9
      A5=-31346.26
      A6=2396.167
      A7=-67.10441
      B=0.43429*ALOG(40.0/PCO2)
      X=PO2*10.0**(0.024*(37.0-TEMP)+0.4*(PHE-7.4)+0.06*B)
      X=26.8*X/P50
      IF(X-10.0) 1,2,2
  1   SAT=0.003683*X + 0.000584*X*X
      GO TO 3
  2   SAT=(X*(X*(X*(X+A3)+A2)+A1))/(X*(X*(X*(X+A7)+A6)+A5)+A4)
  3   SATURA=100.0*SAT
      RETURN
      END
      SUBROUTINE BUFLIN(PO2,PCO2,PHE)
C
C     THIS ROUTINE COMPUTES PH CORRESPONDING TO PCO2'S OF 30 AND 60 TORR
C     TO ALLOW FOR ACID/BASE DISTURBANCES IN CALCULATING O2/CO2 CONTENTS
C
      COMMON HB,HCRIT,TEMP,APH,BPH,APCO2,BPCO2,P50
      Y1=0.003*HB*(100.0-SATURA(PO2,PCO2,PHE))/100.0
      PHX=7.59+Y1-0.2741*ALOG(PCO2/20.0)                              
      DELPH=PHE-PHX
      APH=7.59+DELPH-0.2741*ALOG(APCO2/20.0)
      BPH=7.59+DELPH-0.2741*ALOG(BPCO2/20.0)
      RETURN
      END