
C     PROGRAM P50 
C
C     PROGRAM TO CALCULATE P50 FROM MEASURED SATURATIONS & BLOOD GASES 
C     IT DOES THIS TWO WAYS:  
C
C     1.  FITTING TO KELMAN DISSOCIATION CURVE
C
C     2.  FITTING TO THE HILL EQUN -
C
C           SAT = PO2**N/(PO2**N + P50**N)
C
C     AFTER TRANSFORMING THIS TO:
C
C           LOG[SAT/(1-SAT)]  =  N*LOG(PO2)  -  N*LOG(P50)
C
C     WE LINEARLY BEST FIT LOG[SAT/(1-SAT)] AND LOG(PO2) DATA TO THIS, 
C     SOLVING FOR N AND P50 IN THE USUAL LEAST SQUARES WAY.
C     
      CHARACTER ID,IOK
      DIMENSION PO2(200),PCO2(200),PHE(200),SATM(200),SATC(200),
     1HB(200),HCRIT(200),TEMP(200),FIO2(200),SATMUS(200),ID(40),
     2IUSED(200),PO2V(200)
      COMMON/OXY1/HBX,HCRITX,TEMPX,DP50,PO2VIR
      WRITE(*,10)
  10  FORMAT(' SEND OUTPUT TO PRINTER (1) OR DISK (0) ?')
      READ(*,*)IPR
      IF(IPR.EQ.1) OPEN(3,FILE='PRN')
      IF(IPR.EQ.0) OPEN(3,FILE=' ',STATUS='NEW')
      WRITE(*,20)
  20  FORMAT(' ENTER NAME OF SUBJECT & DATE OF STUDY, 40A1',/,
     1' ****************************************'/)
      READ(*,30)(ID4I),I=1,40)
  30  FORMAT(40A1,/)
      WRITE(3,30)(ID(I),I=1,40)
  40  FORMAT(//)
      WRITE(*,50)
  50  FORMAT(' ENTER NUMBER OF BLOOD GAS SETS')
      READ(*,*)NSETS
      WRITE(*,60)
  60  FORMAT(' ARE YOU RUNNING AN OLD FILE (0) OR CREATING A NEW (1)')
      READ(*,*)ION
      IF(ION.EQ.0) OPEN(6,FILE=' ',STATUS='OLD')
      IF(ION.EQ.1) OPEN(7,FILE=' ',STATUS='NEW')
      WRITE(3,70)
  70  FORMAT(3X,'% INSP O2',4X,'PO2',6X,'PCO2',8X,'PH',4X,'MEAS SAT',
     16X,'HB',6X,'HCRIT',5X,'TEMP'/)
      IF(ION.EQ.0) GO TO 160
      WRITE(7,80)NSETS
  80  FORMAT(I5)
      DO 150 I=1,NSETS
      WRITE(*,90)I
  90  FORMAT(' FOR SET NUMBER:',I5,/)
 100  WRITE(*,110)
 110  FORMAT(' ENTER % INSP O2,PO2,PCO2,PH,SAT,HB,HCT%,TEMP',/)
      READ(*,*)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
 120  FORMAT(7F10.2,F9.2)
      WRITE(*,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
      WRITE(*,130)
 130  FORMAT(' OK (Y) OR NOT (N)')
      READ(*,140)IOK
 140  FORMAT(A1)
      IF(IOK.EQ.'N') GO TO 100
      WRITE(3,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
      WRITE(7,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
 150  CONTINUE
      GO TO 180
 160  CONTINUE
      READ(6,310)NSETS
      DO 170 I=1,NSETS
      READ(6,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
      WRITE(*,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
 170  WRITE(3,120)FIO2(I),PO2(I),PCO2(I),PHE(I),SATM(I),HB(I),HCRIT(I),
     1TEMP(I)
 180  CONTINUE
      WRITE(*,190)
 190  FORMAT(' ENTER LOWER BOUND ON P50')
      READ(*,*)P50
      WRITE(*,200)
 200  FORMAT(' ENTER UPPER BOUND ON P50')
      READ(*,*)UPRP50
C
C     DATA FOR SATURATIONS BELOW AND ABOVE DESIRED LIMITS CAN BE EXCLUDED
C
      WRITE(*,210)
 210  FORMAT(' ENTER LOW AND THEN HI % SATRN LIMITS ON DATA')
      READ(*,*)SLO,SHI
      WRITE(3,220)SLO,SHI
 220  FORMAT(//' P50 CALCS EXCLUDE SATRNS LOWER THAN',F5.1,
     1'  AND HIGHER THAN',F5.1)
      RP50=(UPRP50-P50)*10.0+1.0
      NP50=RP50
      P50=P50-0.1
      SSQMIN=10000.0
      IFLAG=0
      WRITE(*,40)
      WRITE(3,40)
      WRITE(*,230)
      WRITE(3,230)
 230  FORMAT(3X,'TRIAL P50',7X,'SUM OF SQUARES'/)
      DO 260 K=1,NP50
      P50=P50+0.1
      DP50=P50-26.8
      SSQ=0.0
      DO 240 I=1,NSETS
      IF(SATM(I).LT.SLO.OR.SATM(I).GT.SHI) GO TO 240
      HBX=HB(I)
      HCRITX=HCRIT(I)
      TEMPX=TEMP(I)
      SATC(I)=SATURA(PO2(I),PCO2(I),PHE(I))
      SUMDIF=SATM(I)-SATC(I)
      SSQ=SSQ+SUMDIF**2
 240  CONTINUE
      WRITE(3,250)P50,SSQ
      WRITE(*,250)P50,SSQ
 250  FORMAT(F10.1,F20.4)
      IF(SSQ.GT.SSQMIN) IFLAG=IFLAG+1
      IF(SSQ.LT.SSQMIN) RP50=P50
      IF(SSQ.LT.SSQMIN) SSQMIN=SSQ
 260  CONTINUE
      DP50=RP50-26.8
      WRITE(3,270)SSQMIN,RP50
 270  FORMAT(/,' MIN SSQ IS',F8.2,'  AT A P50 OF',F6.1,'  TORR,GIVING:'
     1,//,'  SET #',3X,'VIRT PO2',4X,'MEAS SAT',3X,
     1'KELMAN SAT',5X,'DIFF',5X,'LN PO2',3X,'LN(S/(1-S))')
      DO 280 I=1,NSETS
C      IF(SATM(I).LT.SLO.OR.SATM(I).GT.SHI) GO TO 280
      HBX=HB(I)
      HCRITX=HCRIT(I)
      TEMPX=TEMP(I)
      SATC(I)=SATURA(PO2(I),PCO2(I),PHE(I))
      PO2V(I)=PO2VIR
      DIFSAT=SATM(I)-SATC(I)
      AL=ALOG(PO2V(I))
      AS=ALOG(SATM(I)/(100.0-SATM(I)))
      WRITE(3,290)I,PO2V(I),SATM(I),SATC(I),DIFSAT,AL,AS
 280  CONTINUE
 290  FORMAT(I5,3F12.2,F11.2,F10.2,F12.2)
      WRITE(*,40)
      WRITE(3,40)
      NUSED=0
      DO 300 I=1,NSETS
      IF(SATM(I).LT.SLO.OR.SATM(I).GT.SHI) GO TO 300
      NUSED=NUSED+1
      IUSED(NUSED)=I
      HBX=HB(I)
      HCRITX=HCRIT(I)
      TEMPX=TEMP(I)
      SATC(NUSED)=SATURA(PO2(I),PCO2(I),PHE(I))
      SATMUS(NUSED)=SATM(I)
 300  CONTINUE
 310  FORMAT(I5)
      CALL LINREG(1,NUSED,IUSED,SATMUS,SATC,SLOPE,RINTER)
      WRITE(3,40)
      CALL HILL(NSETS,SLO,SHI,PO2V,SATM,IUSED)
      END
      SUBROUTINE HILL(NSETS,SLO,SHI,PO2,SATM,IUSED)
      DIMENSION PO2(200),SATM(200),IUSED(200)
      DIMENSION X(200),Y(200)
      NP=0
      DO 10 I=1,NSETS
      IF(SATM(I).LT.SLO.OR.SATM(I).GT.SHI) GO TO 10
      NP=NP+1
      X(NP)=ALOG(PO2(I))
      Y(NP)=ALOG(0.01*SATM(I)/(1.0-0.01*SATM(I)))
  10  CONTINUE
      WRITE(3,20)
  20  FORMAT(' LINEAR REGRESSION BETWEEN LOG PO2 & LOG(SAT/(1-SAT))'/)
      WRITE(3,30)
  30  FORMAT(' HILL N IS THE SLOPE; P50 IS EXP(-INTERCEPT/SLOPE)'/)
      CALL LINREG(2,NP,IUSED,X,Y,SLOPE,RINTER)
      REXP=-RINTER/SLOPE
   X  P50=EXP(REXP)
      WRITE(3,40)SLOPE,P50
  40  FORMAT(//' HILL COEFFICIENT =',F6.2,'   P50 =',F6.1)
      RETURN
      END
      FUNCTION SATURA(PO2,PCO2,PHE)
      COMMON/OXY1/HB,HCRIT,TEMP,DP50,XX
      A1=-8532.229
      A2=2121.401
      A3=-67.07399
      A4=935960.9
      A5=-31346.26
      A6=2396.167
      A7=-67.10441
      B=0.43429*ALOG(40.0/PCO2)
      XX=PO2*10.0**(0.024*(37.0-TEMP)+0.4*(PHE-7.4)+0.06*B)
      X=26.8*XX/(26.8+DP50)
      IF(X-10.0) 10,20,20
  10  SAT=0.003683*X + 0.000584*X*X
      GO TO 30
  20  SAT=(X*(X*(X*(X+A3)+A2)+A1))/(X*(X*(X*(X+A7)+A6)+A5)+A4)
  30  SATURA=100.0*SAT
      RETURN
      END
      SUBROUTINE LINREG(IC,NP,IUSED,X,Y,SLOPE,RINTER)
      DIMENSION X(200),Y(200),IUSED(200)
      RNP=NP
      RN1=NP-1
      SUMX=0.0
      SUMY=0.0
      SUMXX=0.0
      SUMYY=0.0
      SUMXY=0.0
      NUSED=0
      DO 10 I=1,NP
      SUMX=SUMX + X(I)
      SUMY=SUMY + Y(I)
      SUMXX=SUMXX + X(I)*X(I)
      SUMYY=SUMYY + Y(I)*Y(I)
  10  SUMXY=SUMXY + X(I)*Y(I)
      XMEAN=SUMX/RNP
      YMEAN=SUMY/RNP
      XSD=SQRT((SUMXX-SUMX*SUMX/RNP)/RN1)
      IF(SUMY.EQ.0.0) GO TO 20
      YSD=SQRT((SUMYY-SUMY*SUMY/RNP)/RN1)
      RNUM=SUMXY-SUMX*SUMY/RNP
      RDENOM=SUMXX-SUMX*SUMX/RNP
      SLOPE=RNUM/RDENOM
      RINTER=YMEAN-SLOPE*XMEAN
      RD1=RDENOM
      RD2=SUMYY-SUMY*SUMY/RNP
      CORCOE=RNUM/SQRT(RD1*RD2)
  20  CONTINUE
      IF(IC.EQ.1)
     1WRITE(3,30)XMEAN,YMEAN,XSD,YSD,SLOPE,RINTER,CORCOE
  30  FORMAT('  MEAN VALUE OF MEASURED   SAT (X) =',F10.1,/,
     1'  MEAN VALUE OF CALCULATED SAT (Y) =',F10.1,/,
     1'  STANDARD DEVIATION OF MEASURED   SATURATIONS =',F10.1,/,
     1'  STANDARD DEVIATION OF CALCULATED SATURATIONS =',F10.1,/,
     1'  SLOPE OF LINEAR REGRESSION LINE =',F10.3,/,
     1'  Y INTERCEPT OF THIS LINE =',F10.3,/,
     1'  CORRELATION COEFFICIENT =',F10.3/)
      IF(IC.EQ.1) WRITE(3,40)
  40  FORMAT(6X,'#',6X,'SAMPLE',5X,'MEAS SAT (X)',3X,'CALC SAT (Y)',4X,
     1'BEST FIT Y',5X,'  DIFF'/)
      IF(IC.EQ.2)
     1WRITE(3,50)XMEAN,YMEAN,XSD,YSD,SLOPE,RINTER,CORCOE
  50  FORMAT('  MEAN VALUE OF      LOG PO2 (X)     =',F10.2,/,
     1'  MEAN VALUE OF LOG(SAT/(1-SAT)) (Y) =',F10.2,/,
     1'  STANDARD DEVIATION OF X =',F10.2,/,
     1'  STANDARD DEVIATION OF Y =',F10.2,/,
     1'  SLOPE OF LINEAR REGRESSION LINE =',F10.3,/,
     1'  Y INTERCEPT OF THIS LINE =',F10.3,/,
     1'  CORRELATION COEFFICIENT =',F10.3//)
      IF(IC.EQ.2) WRITE(3,60)
  60  FORMAT(1X,'#',6X,'SAMPLE',5X,'LOG PO2 (X)',3X,'LOG(SAT/(1-SAT))',
     1' (Y)',2X,'BEST FIT Y',3X,'  DIFF'/)
      DO 70 I=1,NP
      YCALC=RINTER + SLOPE*X(I)
      IF(IC.EQ.1) WRITE(3,80)I,IUSED(I),X(I),Y(I),YCALC,Y(I)-YCALC
      IF(IC.EQ.2) WRITE(3,90)I,IUSED(I),X(I),Y(I),YCALC,Y(I)-YCALC
  70  CONTINUE
  80  FORMAT(I7,I10,3F15.2,F14.2)
  90  FORMAT(I2,I10,2F15.2,7X,F13.2,F12.2)
      RETURN
      END