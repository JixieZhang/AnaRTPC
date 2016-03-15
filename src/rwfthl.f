
      SUBROUTINE RWFTHL(NPT,RF,PF,WFI,ZF,WZF,IOPT,
     1                  VV0,EE0,CH2PH,CH2Z,DEL,DELZ)
      implicit none
C
C-----------------------------------------------------------------------
C! Fast helix fit
C
C   A generalization of the TFTHEL routine to allow it to be called
C   from a routine that contains any list of x and y values XF,YF for a
C   set of NPT points to be fitted.
C
C   Input:  NPT    /I     Number of 3-D points to be fit
C           XF     /R     Array of X-values of points to be fit
C           YF     /R     Array of Y-values of points to be fit
C           RF     /R     Array of R-values of points to be fit
C           PF     /R     Array of PHI-values of points to be fit
C           WFI    /R     Array of 1/(sig(rphi))**2 for each point
C           ZF     /R     Array of Z-values of points to be fit
C           WZF    /R     Array of 1/(sig(z))**2 for each point
C           IOPT = 0 -> DISTANCE**2=X**2+Y**2 MINIMISED
C                  1 -> WEIGHTED WITH 1/SIMA(R*PHI)**2
C                  2 -> ERROR MATRIX CALCULATED
C                  3 -> 3-DIMENSIONAL ITERATION
C  OUTPUT:   VV0 = 1/R*CHARGE    POS. IF CLOCKWISE
C                  TAN(LAMBDA)  {=DZ/DS}TAN(ANGLE TO X,Y PLANE)
C                  PHI0         {0,2PI} ANGLE TO X-AXIS at R=D0
C                  D0*SIGN      [CM]    MINIMAL DIST. TO Z-AXIS,
C                                       +VE IF AXIS ENCIRCLED
C                  Z0           [CM]    Z POS AT R=D0
C          EE0 = INVERSE OF ERROR MATRIX IN TRIANG. FORM
C          CH2PH = CHI SQUARED = SUM (PHI DEVIATIONS/ERRORS)**2
C          CH2Z  = CHI SQUARED = SUM (Z DEVIATIONS/ERRORS)**2
C          DEL = ??
c          DELZ= ??
C  NOTE: DEGREES OF FREEDOM = 2*NPT-5
C----------------------------------------------------------------
C     BASED ON  SUBROUTINE CIRCLE
C     REFERENCE:  COMPUTER PHYSICS COMMUNICATIONS VOL 33,P329
C
C   AUTHORS:  N. CHERNOV, G. OSOSKOV & M. POPPE
C   Modified by:  Fred Weber, 8 Jun 1989
C
C-----------------------------------------------------------------
C** +SEQ,SHIMPNON.
C** +SEQ,SHCCTL.
C** +SEQ,SHCONST.
C
C** +SEQ,RWDIM.
C** +SEQ,RWCNTC.
      INTEGER ITMAX
      REAL EPS
      PARAMETER (EPS = 1.0E-16, ITMAX =15)
      INTEGER MXRWPT
      PARAMETER (MXRWPT = 200)
      INTEGER IOWRIT
      PARAMETER(IOWRIT=6)
      REAL PI, TWOPI, PIBY2
      PARAMETER(PI=3.1415927)
      PARAMETER(TWOPI= 6.2831854)
      PARAMETER(PIBY2= 1.57074635)
      REAL   PF(MXRWPT),RF(MXRWPT),PC(MXRWPT),SP2(MXRWPT),
     1       VV0(*),EE0(*),VV1(5),
     2       DEL(MXRWPT),SXY(MXRWPT),ZF(MXRWPT),WZF(MXRWPT),
     3       SS0(MXRWPT),EEE(MXRWPT),DELZ(MXRWPT),
     4       GRAD(5),COV(15),DV(5),GGRAD(5),
     5       DELN(MXRWPT),DELZN(MXRWPT),WFI(MXRWPT)
      DOUBLE PRECISION XF(MXRWPT),YF(MXRWPT),WF(MXRWPT),XMID,YMID
      DOUBLE PRECISION ALF,ALM,A0,A1,A2,A22,BEM,BET,CUR,
     1   DD,DEN,DET,DY,D2,F,FACT,FG,F1,G,GAM,GAM0,GMM,G1,
     2   H,H2,P2,Q2,RM,RN,
     3   XA,XB,XD,XI,XM,XX,XY,X1,X2,DEN2,
     4   YA,YB,YD,YI,YM,YY,   Y1,Y2,WN,SA2B2,DD0,CU2,PHIC
C **
      INTEGER NPT,I,N,IOPT,ITER,NRANK
      REAL CHI2,CH2PH,CH2Z,RR0,ASYM,SST,PH0,CHECK
      REAL AA0,OME,GG0,HH0,FF0,SUMS,SUMSS,SUMZ,SUMSZ,SUMW
      REAL DENOM,DZDS,ZZ0,ETA,DFD,DFO,DPD,DPO,GGG,DZA,DZD,DZO,CHI1
C
      IF (NPT .LE. 2)  THEN
         WRITE
     *   (IOWRIT,'('' RWFTHL: Cannot fit < 3 points! quitting...'')')
         RETURN
      END IF
      IF (NPT .GT. MXRWPT)  THEN
         WRITE
     *   (IOWRIT,'('' RWFTHL: Cannot fit > MXRWPT points! '')')
         RETURN
      END IF
      DO 11 I=1,NPT
         XF(I) = RF(I)*COS(PF(I))
         YF(I) = RF(I)*SIN(PF(I))
         WF(I) = WFI(I)
   11 CONTINUE
      N = NPT
      XM = 0.
      YM = 0.
      DO 1   I=1,15
         EE0(I)=0.0
  1   CONTINUE
      DO 2   I=1,5
         GRAD(I)=0.0
         VV0(I) =0.0
  2   CONTINUE
      CHI2=0.0
      CH2PH = 0.
      CH2Z = 0.
      DO 10 I = 1,N
         SP2(I) = WF(I)*(RF(I)*RF(I))
   10 CONTINUE
      IF(IOPT.EQ.0) THEN
         DO 100 I= 1, N
            WZF(I)= 1.0
            WF(I) = 1.0
            XM    = XM + XF(I)
            YM    = YM + YF(I)
  100    CONTINUE
         RN = 1./FLOAT(N)
      ELSE
         WN=0.0
         DO 101 I= 1, N
            XM = XM + XF(I)*WF(I)
            YM = YM + YF(I)*WF(I)
            WN = WN + WF(I)
  101    CONTINUE
         RN = 1.D0/WN
      END IF
C **
      XM = XM * RN
      YM = YM * RN
      X2 = 0.
      Y2 = 0.
      XY = 0.
      XD = 0.
      YD = 0.
      D2 = 0.
      DO 102 I= 1, N
         XI = XF(I) - XM
         YI = YF(I) - YM
         XX = XI**2
         YY = YI**2
         X2 = X2 + XX*WF(I)
         Y2 = Y2 + YY*WF(I)
         XY = XY + XI*YI*WF(I)
         DD = XX + YY
         XD = XD + XI*DD*WF(I)
         YD = YD + YI*DD*WF(I)
         D2 = D2 + DD**2*WF(I)
  102 CONTINUE
C **
      X2 = X2*RN
      Y2 = Y2*RN
      XY = XY*RN
      D2 = D2*RN
      XD = XD*RN
      YD = YD*RN
      F = 3.D0*X2 + Y2
      G = 3.D0*Y2 + X2
      FG = F*G
      H = XY + XY
      H2 = H**2
      P2 = XD**2
      Q2 = YD**2
      GAM0 = X2 + Y2
         FACT = GAM0**2
      A2 = (FG-H2-D2)/FACT
         FACT = FACT*GAM0
      A1 = (D2*(F+G) - 2.D0*(P2+Q2))/FACT
         FACT = FACT*GAM0
      A0 = (D2*(H2-FG) + 2.D0*(P2*G + Q2*F) - 4.D0*XD*YD*H)/FACT
      A22 = A2 + A2
      YB = 1.0E30
      ITER = 0
      XA = 1.D0
C **                MAIN ITERATION
  103 YA = A0 + XA*(A1 + XA*(A2 + XA*(XA-4.D0)))
      IF (ITER .GE. ITMAX)                      GO TO 105
      DY = A1 + XA*(A22 + XA*(4.D0*XA - 12.D0))
      XB = XA - YA/DY
      IF (ABS(YA).GT.ABS(YB)) XB=0.5D0*(XB+XA)
      IF (ABS(XA-XB) .LT. EPS)                  GO TO 105
      XA = XB
      YB = YA
      ITER = ITER + 1
      GO TO 103
C **
  105 CONTINUE
      GAM = GAM0*XB
      F1 = F - GAM
      G1 = G - GAM
      X1 = XD*G1 - YD*H
      Y1 = YD*F1 - XD*H
      DET = F1*G1 - H2
      DEN2= 1.D0/(X1**2 + Y1**2 + GAM*DET**2)
      IF(DEN2.LE.0.D0)                GO TO 999
      DEN = DSQRT(DEN2)
      CUR = DET*DEN                  + 0.0000000001D0
      ALF = -(XM*DET + X1)*DEN
      BET = -(YM*DET + Y1)*DEN
      RM = XM**2 + YM**2
C      GAM = ((RM-GAM)*DET + 2.D0*(XM*X1 + YM*Y1))*DEN*0.5D0
C
C--------> CALCULATION OF STANDARD CIRCLE PARAMETERS
C          NB: CUR IS ALWAYS POSITIVE
      ASYM = BET*XM-ALF*YM
      SST = 1.0
      IF(ASYM.LT.0.0) SST=-1.0
      RR0 = SST*CUR
      IF((ALF*ALF+BET*BET).LE.0.D0)              GO TO 999
      SA2B2 = 1.D0/DSQRT(ALF*ALF+BET*BET)
      DD0 = (1.D0-1.D0/SA2B2)/CUR
      PHIC = DASIN(ALF*SA2B2)+PIBY2
      IF(BET.GT.0)    PHIC=TWOPI-PHIC
      PH0 = PHIC+PIBY2
      IF(RR0.LE.0)    PH0=PH0-PI
      IF(PH0.GT.TWOPI) PH0=PH0-TWOPI
      IF(PH0.LT.0.0)  PH0=PH0+TWOPI
      VV0(1) = RR0
      VV0(3) = PH0
      VV0(4) = DD0
      CHECK=SST*RR0*DD0
      IF(CHECK.EQ.1.) THEN
        DD0=DD0-.007
        VV0(4)=DD0
      ENDIF
C
C-----> CALCULATE PHI DISTANCES TO MEASURED POINTS
C
      AA0 =SST
      OME =RR0
      GG0 = OME*DD0-AA0
C      HH0 = 0.0000000001
C      IF(ABS(GG0).LT.(1.0/HH0)) HH0=1.0/GG0
      HH0=1.0/GG0
      DO 300 I=1,N
        ASYM   = BET*XF(I)-ALF*YF(I)
        SS0(I) =1.0
        IF(ASYM.LT.0.0) SS0(I)=-1.0
        FF0   = OME*(RF(I)*RF(I)-DD0*DD0)/(2.0*RF(I)*GG0) + DD0/RF(I)
C
        IF(FF0.LT.-1.0) FF0 = -1.0
        IF(FF0.GT.1.0)  FF0 = 1.0
C
        DEL(I)= PH0 + (SS0(I)-AA0)*PIBY2 + SS0(I)*ASIN(FF0) - PF(I)
        IF(DEL(I).GT.PI) DEL(I)=DEL(I)-TWOPI
        IF(DEL(I).LT.-PI)DEL(I)=DEL(I)+TWOPI
  300 CONTINUE
C
C-----> FIT STRAIGHT LINE IN S-Z
C
      DO 350 I=1,N
        EEE(I) = 0.5*VV0(1)
     1          *SQRT(ABS( (RF(I)*RF(I)-VV0(4)*VV0(4))
     2                /(1.0-AA0*VV0(1)*VV0(4))     ) )
C
        IF(EEE(I).GT.0.99990)  THEN
           EEE(I)= 0.99990
C           WRITE(6,901)
C 901       FORMAT(" TRACK CIRCLES TOO MUCH FOR THIS CODE BAD DZDS")
        ENDIF
        IF(EEE(I).LT.-0.99990) THEN
           EEE(I)= -0.99990
C           WRITE(6,901)
        ENDIF
C
        SXY(I)=2.0*ASIN(EEE(I))/OME
  350 CONTINUE
      SUMS  = 0.0
      SUMSS = 0.0
      SUMZ  = 0.0
      SUMSZ = 0.0
      SUMW  = 0.0
      DO 360 I=1,N
        SUMW  = SUMW  +                 WZF(I)
        SUMS  = SUMS  + SXY(I)        * WZF(I)
        SUMSS = SUMSS + SXY(I)*SXY(I) * WZF(I)
        SUMZ  = SUMZ  + ZF(I)         * WZF(I)
        SUMSZ = SUMSZ + ZF(I)*SXY(I)  * WZF(I)
  360 CONTINUE
      DENOM = SUMW*SUMSS - SUMS*SUMS
      IF(ABS(DENOM).LT.1.E-6) DENOM=SIGN(1.E-6,DENOM)
      DZDS  = (SUMW*SUMSZ-SUMS*SUMZ) /DENOM
      ZZ0   = (SUMSS*SUMZ-SUMS*SUMSZ)/DENOM
      VV0(2)= DZDS
      VV0(5)= ZZ0
C
C-----> CALCULATION CHI**2
C
      DO 370 I=1,N
        DELZ(I)= ZZ0+DZDS*SXY(I)-ZF(I)
        CH2PH = CH2PH + SP2(I)*DEL(I)*DEL(I)
        CH2Z = CH2Z + WZF(I)*DELZ(I)*DELZ(I)
        CHI2 = CH2PH + CH2Z
c        write(6,1011) I,ZZ0+DZDS*SXY(I)-ZF(I),SP2(I)*DEL(I)*DEL(I),
c    $       WZF(I)*DELZ(I)*DELZ(I)
 1011   FORMAT("CHI2 CONTRIBUTIONS: ",I3,1x,E10.4,1x,E10.4,1x,E10.4)
  370 CONTINUE
c      write(6,1012) ch2ph, ch2z, chi2
 1012 format("chi2s are ",3(e10.4,1x))
C
      IF(IOPT.LT.2)                          RETURN
C
C-----> CALCULATION OF THE ERROR MATRIX
C
      DO 380 I=1,N
        FF0 = OME*(RF(I)*RF(I)-DD0*DD0)/(2.0*RF(I)*GG0) + DD0/RF(I)
        IF (FF0.GT. 0.99990) FF0 =  0.99990
        IF (FF0.LT.-0.99990) FF0 = -0.99990
        ETA = SS0(I)/SQRT(ABS((1.0+FF0)*(1.0-FF0)))
        DFD = (1.0+HH0*HH0*(1.0-OME*OME*RF(I)*RF(I)))/(2.0*RF(I))
        DFO = -AA0*(RF(I)*RF(I)-DD0*DD0)*HH0*HH0/(2.0*RF(I))
        DPD = ETA*DFD
        DPO = ETA*DFO
C-----> DERIVATIVES OF Z COMPONENT
        GGG = EEE(I)/SQRT(ABS( (1.0+EEE(I))*(1.0-EEE(I))))
        DZA = SXY(I)
        CHECK=RF(I)*RF(I)-VV0(4)*VV0(4)
        IF(CHECK.EQ.0.) CHECK=2.*.007
        DZD = 2.0*( VV0(2)/VV0(1) ) * ABS( GGG )
     1       *(  0.5*AA0*VV0(1)/( 1.0-AA0*VV0(4)*VV0(1) )
     2           -VV0(4)/CHECK    )
        DZO = -VV0(2)*SXY(I)/VV0(1)
     1        +VV0(2)* GGG/( VV0(1)*VV0(1))
     2         *( 2.0+ AA0*VV0(1)*VV0(4)/(1.0-AA0*VV0(1)*VV0(4)) )
C-----> ERROR MARTIX
        EE0(1) = EE0(1) + SP2(I)*  DPO*DPO  + WZF(I) * DZO*DZO
        EE0(2) = EE0(2)                     + WZF(I) * DZA*DZO
        EE0(3) = EE0(3)                     + WZF(I) * DZA*DZA
        EE0(4) = EE0(4) + SP2(I)*  DPO
        EE0(5) = EE0(5)
        EE0(6) = EE0(6) + SP2(I)
        EE0(7) = EE0(7) + SP2(I)*  DPO*DPD  + WZF(I) * DZO*DZD
        EE0(8) = EE0(8)                     + WZF(I) * DZA*DZD
        EE0(9) = EE0(9) + SP2(I)*      DPD
        EE0(10)= EE0(10)+ SP2(I)*  DPD*DPD  + WZF(I) * DZD*DZD
        EE0(11)= EE0(11)                    + WZF(I) * DZO
        EE0(12)= EE0(12)                    + WZF(I) * DZA
        EE0(13)= EE0(13)
        EE0(14)= EE0(14)                    + WZF(I) * DZD
        EE0(15)= EE0(15)                    + WZF(I)
C-----> GRADIENT VECTOR
        GRAD(1)=GRAD(1) - DEL(I) *SP2(I)*DPO - DELZ(I)*WZF(I)*DZO
        GRAD(2)=GRAD(2) -                      DELZ(I)*WZF(I)*DZA
        GRAD(3)=GRAD(3) - DEL(I) *SP2(I)
        GRAD(4)=GRAD(4) - DEL(I) *SP2(I)*DPD - DELZ(I)*WZF(I)*DZD
        GRAD(5)=GRAD(5) -                      DELZ(I)*WZF(I)
  380 CONTINUE
C
      IF(IOPT.LT.3) RETURN
C
C------> NEWTON'S NEXT GUESS
      DO 401 I=1,15
  401   COV(I)=EE0(I)
      CALL RWSMIN(COV,VV1,5,0,NRANK)
      CALL RWSMAV(DV,COV,GRAD,5)
      DO 402 I=1,5
  402   VV1(I)=VV0(I)+DV(I)
C------>  NEW DIFFERENCES IN PHI AND Z
      GG0 = VV1(1)*VV1(4)-AA0
      DO 410 I=1,N
        FF0    = VV1(1)*(RF(I)*RF(I)-VV1(4)**2)
     1          /(2.0*RF(I)*GG0) + VV1(4)/RF(I)
C
      IF(FF0.GT.1.0)  FF0 = 1.0
      IF(FF0.LT.-1.0) FF0 = -1.0
C
        DELN(I) = VV1(3) + (SS0(I)-AA0)*PIBY2+SS0(I)*ASIN(FF0)-PF(I)
        IF(DELN(I).GT.PI) DELN(I)=DELN(I)-TWOPI
        IF(DELN(I).LT.-PI)DELN(I)=DELN(I)+TWOPI
        EEE(I) = 0.5*VV1(1)
     1          *SQRT(ABS( (RF(I)*RF(I)-VV1(4)*VV1(4))
     2                /(1.0-AA0*VV1(1)*VV1(4))     ))
        IF(EEE(I).GT.0.99990)  EEE(I)= 0.99990
        IF(EEE(I).LT.-0.99990) EEE(I)= -0.99990
        SXY(I) = 2.0*ASIN(EEE(I))/VV1(1)
        DELZN(I)= VV1(5)+VV1(2)*SXY(I)-ZF(I)
  410 CONTINUE
C-----> CALCULATION CHI**2
      CHI1 = 0.0
      CH2PH = 0.0
      CH2Z = 0.0
      DO 420 I=1,N
        CHI1   = CHI1  + SP2(I)*DELN(I)*DELN(I)
     1                 + WZF(I)*DELZN(I)*DELZN(I)
        CH2PH = CH2PH + SP2(I)*DELN(I)*DELN(I)
        CH2Z = CH2Z + WZF(I)*DELZN(I)*DELZN(I)
  420 CONTINUE
c      write(6,1012) ch2ph, ch2z, chi2
      IF(CHI1.LT.CHI2) THEN
         DO 430 I=1,5
  430    VV0(I) = VV1(I)
C         CHI2 = CHI1
      END IF
      RETURN
  999 CONTINUE
C------> JUMP HERE IF SOMETHING GOES CRAZY
      CH2PH = 1.0E30
      CH2Z  = 1.0E30
      RETURN
      END
C** +DECK,RWSMAV.
      SUBROUTINE RWSMAV(R,A,V,N)
C -----------------------------------------------------------
C!    R(N) = A(N,N)*V(N)
C!    AUTHOR: MARTIN POPPE
C -----------------------------------------------------------
C
C** +SEQ,SHIMPNON.
      INTEGER N,I,II,K,KK,IND
      REAL R(*),A(*),V(*)
C---> ADDRESSE IN TRIAGULAR MATRIX, ROW II COLUMN KK
      IND(II,KK)=(II*II-II)/2+KK
C---> ADDRESSE IN TRIAGULAR MATRIX, ROW II COLUMN KK
      DO 30 I=1,N
        R(I) = 0.0
        DO 20 K=1,N
          IF(I.GE.K)THEN
            R(I)=R(I)+A(IND(I,K))*V(K)
          ELSE
            R(I)=R(I)+A(IND(K,I))*V(K)
          END IF
  20    CONTINUE
  30  CONTINUE
      RETURN
      END

C** +DECK,RWSMIN.
      SUBROUTINE RWSMIN(V,B,N,M,NRANK)
C     AUTHOR: V.BLOBEL
C
C     SUBROUTINE SMINV
C     ----------------
C!    OBTAIN SOLUTION OF A SYSTEM OF LINEAR EQUATIONS V *  X  =  B  WITH
C!    SYMMETRIC MATRIX V AND INVERSE (FOR M =  1)  OR  MATRIX  INVERSION
C!    ONLY (FOR M = 0)
C
C                   - - - -
C        CALL SMINV(V,B,N,M,NRANK)
C                   - -     -----
C
C           V = SYMMETRIC N-BY-N MATRIX IN SYMMETRIC STORAGE MODE
C               V(1) = V11, V(2) = V12, V(3) = V22, V(4) = V13, . . .
C               REPLACED BY INVERSE MATRIX
C           B = N-VECTOR   (FOR M = 0 USE A DUMMY ARGUMENT)
C               REPLACED BY SOLUTION VECTOR
C           M = SEE ABOVE
C
C
C     METHOD OF SOLUTION IS BY ELIMINATION SELECTING THE  PIVOT  ON  THE
C     DIAGONAL EACH STAGE. THE RANK OF THE MATRIX IS RETURNED IN  NRANK.
C     FOR NRANK NE N, ALL REMAINING  ROWS  AND  COLS  OF  THE  RESULTING
C     MATRIX V AND THE CORRESPONDING ELEMENTS OF  B  ARE  SET  TO  ZERO.
C     SMINV USES A WORK ARRAY OF 2*N WORDS IN COMMON/MATCOM/. FOR N> 200
C     THE USER HAS TO DEFINE COMMON/MATCOM/ WITH 2*N WORDS.
C
C** +SEQ,SHIMPNON.
      INTEGER N,M,NRANK,I,II,NI,K,KK,J,JJ,JL,JK,LK,L,IJ
      REAL V(*),B(*),VKK,C,D,E
      REAL EPS,DR(2,200)
      DATA EPS/1.E-6/
C
C     CONSTRUCT TABLE
C
      DO 10 I=1,N
   10 DR(1,I)=1.0
      NI=N
   14 II=0
      DO 16 I=1,N
      II=II+I
   16 DR(2,I)=ABS(V(II))
C
C     LOOP BEGIN
C
      NRANK=N-NI
      DO 60 I=1,NI
C
C     SEARCH FOR PIVOT AND TEST FOR LINEARITY AND ZERO MATRIX
C
      K=0
      JJ=0
      VKK=0.0
      DO 20 J=1,N
      JJ=JJ+J
      IF(DR(1,J).EQ.0.0) GOTO 20
      IF(ABS(V(JJ)).LE.VKK) GOTO 20
      IF(ABS(V(JJ)).LT.EPS*DR(2,J)) GOTO 20
      VKK=ABS(V(JJ))
      K=J
      KK=JJ
   20 CONTINUE
      IF(K.EQ.0) GOTO 80
C
C     PREPARATION FOR ELIMINATION
C
      NRANK=NRANK+1
      DR(1,K)=0.0
      D=1.0/V(KK)
      V(KK)=-D
      IF(M.EQ.1) B(K)=B(K)*D
      JK=KK-K
      JL=0
C
C          ELIMINATION
C
      DO 50 J=1,N
      IF(J-K) 24,22,26
   22 JK=KK
      JL=JL+J
      GOTO 50
   24 JK=JK+1
      GOTO 28
   26 JK=JK+J-1
   28 E=V(JK)
      V(JK)=D*E
      IF(M.EQ.1) B(J)=B(J)-B(K)*E
      LK=KK-K
      DO 40 L=1,J
      JL=JL+1
      IF(L-K) 34,32,36
   32 LK=KK
      GOTO 40
   34 LK=LK+1
      GOTO 38
   36 LK=LK+L-1
   38 V(JL)=V(JL)-V(LK)*E
   40 CONTINUE
   50 CONTINUE
   60 CONTINUE
C
C          CHANGE SIGN
C
      IJ=0
      DO 70 I=1,N
      DO 70 J=1,I
      IJ=IJ+1
   70 V(IJ)=-V(IJ)
      GOTO 100
C
C          CLEAR REST OF MATRIX
C
   80 IJ=0
      DO 90 I=1,N
      IF(M.EQ.1.AND.DR(1,I).NE.0.0) B(I)=0.0
      DO 90 J=1,I
      IJ=IJ+1
      IF(DR(1,I)+DR(1,J).NE.0.0) V(IJ)=0.0
   90 V(IJ)=-V(IJ)
  100 RETURN
      END


