C---------------------------------------------------------------------
C  program PHSH2.F
C  MAIN PROGRAM FOR CALCULATION OF PHASE SHIFTS USING POTENTIAL
C---------------------------------------------------------------------
      PROGRAM PHSH2
      CHARACTER(LEN=255)	:: INPUT, ARG, DATAFILE, OUTPUT
      CHARACTER(LEN=255)    :: STRING(3)
      INTEGER 				:: I, IPHSH

C  SPECIFY DEFAULTS
      IPHSH = 0 !Rel=0, Cav=1, Wil=2
      DATAFILE = "dataph.d"
      INPUT = "mufftin.d"
      OUTPUT = "phasout"

C  GET PROGRAM ARGUMENTS
      I = 1
      DO
        CALL GETARG(I, ARG)
        IF ((ARG.EQ.'-i').OR.(ARG.EQ.'--input')) THEN
          I = I+1
          CALL GETARG(I, INPUT)
        ENDIF
        IF ((ARG.EQ.'-d').OR.(ARG.EQ.'--datafile')) THEN
          I = I+1
          CALL GETARG(I, DATAFILE)
        ENDIF
        IF ((ARG.EQ.'-c').OR.(ARG.EQ.'--cav')) THEN
          IPHSH = 1
        ENDIF
        IF ((ARG.EQ.'-r').OR.(ARG.EQ.'--rel')) THEN
          IPHSH = 0
        ENDIF
        IF ((ARG.EQ.'-w').OR.(ARG.EQ.'--wil')) THEN
          IPHSH = 2
        ENDIF
        IF ((ARG.EQ.'-o').OR.(ARG.EQ.'--output')) THEN
          I = I+1
          CALL GETARG(I, OUTPUT)
        ENDIF
        IF ((ARG.EQ.'-h').OR.(ARG.EQ.'--help')) THEN
	      WRITE(*,*) 'phsh2 usage:-'
          WRITE(*,*)
          STRING(1) = 'phsh2 -i <file> -d <file> -o <file>'
          STRING(2) = ' [-r|-c|-w] '
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))
		  WRITE(*,*) ''
          WRITE(*,*) 'where:-'
          STRING(1) = '-i or --input <file>'
          STRING(2) = ' specifies mufftin file path'
          STRING(3) = ' (default: "mufftin.d")'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))//TRIM(STRING(3))
          STRING(1) = '-d or --datafile <file>'
          STRING(2) = ' specifies output data file dump'
          STRING(3) = ' (default: "dataph.d")'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))//TRIM(STRING(3))
          STRING(1) = '-o or --output <output> specifies the phase '
          STRING(2) = ' shift output file path (default: "phasout")'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))
          STRING(1) = '-c or --cav'
          STRING(2) = ' use CAVLEED phase shift calculations'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))//TRIM(STRING(3))
          STRING(1) = '-r or --rel use relativistic phase shift'
          STRING(2) = ' calculations'
          WRITE(*,*) TRIM(STRING(1))//TRIM(STRING(2))
    	  STRING(1) = "-w or --wil use Williams' phase shift calculations"
          WRITE(*,*) TRIM(STRING(1))
          STRING(1) = '-h or --help print help and exit'
          WRITE(*,*) TRIM(STRING(1))
          STOP
        ENDIF
        IF (I.GE.IARGC()) THEN
          EXIT
        ENDIF
        I = I+1
      END DO

C  CALL DESIRED PHASE SHIFT CALCULATION PROGRAM
      IF (IPHSH.EQ.0) THEN
        CALL PHSH_REL(INPUT, OUTPUT, DATAFILE)
      ELSE IF (IPHSH.EQ.1) THEN
        CALL PHSH_CAV(INPUT, OUTPUT, DATAFILE)
      ELSE IF (IPHSH.EQ.2) THEN
        CALL PHSH_WIL(INPUT, OUTPUT, DATAFILE)
      ENDIF

      STOP
      END

C---------------------------------------------------------------------
C  subroutine PHSH2CAV
C---------------------------------------------------------------------
C
C
C  POTENTIAL-TO-PHASE-SHIFT CALCULATION(CAVLEED PACKAGE)
C
C  USES LOUCKS GRID (E.G. AS SUPPLIED BY THE MUFFIN-TIN POTENTIAL
C  PROGRAM).  ENERGIES INPUT IN HARTREES.
      SUBROUTINE PHSH_CAV(MUFFTIN_FILE, PHASOUT_FILE, DATAPH_FILE)
      CHARACTER(LEN=*), INTENT(IN)	:: MUFFTIN_FILE
      CHARACTER(LEN=*), INTENT(IN)	:: PHASOUT_FILE, DATAPH_FILE
      DIMENSION V(250),RX(250),PHS(20)
      REAL*8 NAME(2),MTZ,delstore(401,15),estore(401)
C
C First input channels
C
      OPEN (UNIT=4,FILE=MUFFTIN_FILE,STATUS='OLD')
C
C Now output channels
C
      OPEN (UNIT=6,FILE='zph.o',STATUS='UNKNOWN')
      OPEN (UNIT=9,FILE=PHASOUT_FILE,STATUS='UNKNOWN')
      OPEN (UNIT=8,FILE=DATAPH_FILE,STATUS='UNKNOWN')
C
C standard values for phase shifts calculation
C
      write(8,110)
      emin=1.
      emax=12.
      estep=.25
      ianz=(emax-emin)/estep +1.01
      nl=12
      READ(4,103)NR
      DO 2  KKK=1,NR
      READ(4,100)NAME
      READ(4,101)Z,RMT,MTZ
      READ(4,103)NTAB
      MTZ=MTZ/2.
      DO 19 IX=1,NTAB
19    READ(4,219)RX(IX),V(IX)
      WRITE(6,200)NAME,Z,RMT,MTZ
      WRITE(9,181)(NAME(I),I=1,2)
  181 FORMAT('NON-RELATIVISTIC PHASE SHIFTS FOR ',2A4)
      WRITE(9,1030) emin,estep,IANZ,nl
 1030 FORMAT(2F9.4,2(2X,I3))
      E=EMIN
      ncount=0
1     E=2.*E
      ncount=ncount+1
      CALL PS(V,RX,NTAB,RMT,E,PHS,NL)
      E=0.5*E
      WRITE(6,201)E,(PHS(L),L=1,NL)
      WRITE(9,1040)E*27.21,(PHS(L),L=1,NL)
 1040    FORMAT(F9.4,8F8.4)
C  store phase shifts
      do 144 kk=1,nl
144   delstore(ncount,kk)=phs(kk)
      estore(ncount)=E
      E=E+ESTEP
      IF(E.LE.EMAX)GOTO 1
C  write phase shifts as function of energy for plotting
      do 145 kk=1,nl
         write(8,107) kk-1
         do 146 ii=1,ncount
146           write(8,*) estore(II),delstore(II,kk)
145        write(8,*)
2     continue
      STOP
71    FORMAT(1F7.4,/,10F7.4)
100   FORMAT(2A8)
101   FORMAT(3F8.4)
102   FORMAT(I4/(5E14.5))
103   FORMAT(I4)
107   format(3H"L=,i2)
110   format (11HTitleText: ,8HDELTA(E))
200   FORMAT(18H1PHASE SHIFTS FOR ,2A8,2X,15HATOMIC NUMBER?,,F6.1/
     + 19H0MUFFIN-TIN RADIUS?,F8.4,6X,
     + 23H MUFFIN-TIN ZERO LEVEL?,F8.4,9H HARTREES)
201   FORMAT(8H0ENERGY?,F8.4,9H HARTREES/(10F12.5))
202   FORMAT(2A8/2F8.4)
203   FORMAT(10F8.4)
219   FORMAT(2E14.5)
      RETURN
      END
C***********************************************************************
      SUBROUTINE PS(V,RX,NGRID,RAD,E,PHS,NL)
      DIMENSION V(NGRID),RX(NGRID),PHS(NL)
      DIMENSION WF(250),BJ(25),BN(25),XR(10),FR(5)
C  SUBROUTINE TO CALCULATE NL PHASE SHIFTS (L=0,NL-1) FOR AN
C  ATOMIC POTENTIAL TABULATED ON THE LOUCKS RADIAL GRID.
C     V?  ATOMIC POTENTIAL (RYDBERGS)
C     RX?  LOUCKS' EXPONENTIAL GRID  RX(I)=EXP(-8.8+0.05(I-1))
C     NGRID?  NUMBER OF TABULATED POINTS
C     RAD?  LIMIT OF INTEGRATION OF SCHRODINGER EQUATION (A.U.)
C     E?  ENERGY (RYDBERGS)
C     PHS?  PHASE SHIFTS FOR L=0 TO NL-1
C  REFERENCE? LOUCKS T L, (1967), A.P.W. METHOD, BENJAMIN, NY.
C***********************************************************************
      DATA PI,DX,DAC,DB/3.141592653589,0.05,2.083333E-04,
     + 2.083333E-03/
      INDEX(X)=20.*(ALOG(X)+8.8)+2.
C  TABULATION OF SPHERICAL BESSEL FUNCTIONS IN BJ AND BN
      ES=SQRT(E)
      X=ES*RAD
      Z=X
      LL=NL+1
      CALL CALCBF(BJ,BN,LL,X)
C  INTEGRATION OF THE RADIAL SCHRODINGER EQUATION BY THE NUMEROV
C  METHOD (SEE LOUCKS P56 FF).  WF CONTAINS WAVE FUNCTION X RADIUS,
C  ON THE LOUCKS GRID
      X1=EXP(-8.8)
      DO 6 L1=1,NL
      FL=L1-1
      FL2=FL+0.5
      FF=FL2*FL2
      Y1=X1**FL2
      Y2=EXP(DX*FL2)*Y1
      WRITE(6,60)FL,Y1,Y2
60    FORMAT('0L?',F5.1,5X,'Y1,Y2?',2E14.5)
      GAM1=FF+RX(1)*RX(1)*(V(1)-E)
      GAM2=FF+RX(2)*RX(2)*(V(2)-E)
      WF(1)=Y1*SQRT(RX(1))
      WF(2)=Y2*SQRT(RX(2))
      DO 2 IX=3,NGRID
      GAM=FF+RX(IX)*RX(IX)*(V(IX)-E)
      A=1.-DAC*GAM
      B=-2.-DB*GAM2
      C=1.-DAC*GAM1
      YN=-(B*Y2+C*Y1)/A
      WF(IX)=YN*SQRT(RX(IX))
      Y1=Y2
      Y2=YN
      GAM1=GAM2
2     GAM2=GAM
C  LAGRANGIAN INTERPOLATION FOR WAVEFUNCTION AND DERIVATIVE AT
C  RADIUS X.  WFN HOLDS WAVEFUNCTION X RADIUS, AND DWFN DERIVATIVE X
C  RADIUS
      X=RAD
      JR=INDEX(RAD)
      DO 3 J=1,5
      XR(J)=RX(JR-5+J)
      XR(J+5)=XR(J)
3     FR(J)=WF(JR-5+J)
      WFN=0.
      DWFN=0.
      A=(X-XR(1))*(X-XR(2))*(X-XR(3))*(X-XR(4))*(X-XR(5))
      DO 5 I=1,5
      TERM=A/(X-XR(I))/(XR(I)-XR(I+1))/(XR(I)-XR(I+2))
     +      /(XR(I)-XR(I+3))/(XR(I)-XR(I+4))
      SUM=0.
      DO 4 J=1,5
      IF(I.EQ.J)GOTO 4
      SUM=SUM+TERM/(X-XR(J))
4     CONTINUE
      WFN=WFN+TERM*FR(I)
5     DWFN=DWFN+SUM*FR(I)
C  LOGARITHMIC DERIVATIVE
      DLOGA=DWFN/WFN-1./RAD
C  PHASE SHIFTS
      X=ES*RAD
      A=FL*BJ(L1)/X-BJ(L1+1)
      B=FL*BN(L1)/X-BN(L1+1)
      A=ES*A-DLOGA*BJ(L1)
      B=ES*B-DLOGA*BN(L1)
      PHS(L1)=PI/2.
      IF(ABS(B).GT.1.0E-8)PHS(L1)=ATAN(A/B)
      WRITE(6,78)PHS(L1)
78    FORMAT('0PHASE SHIFT?',F10.4)
6     CONTINUE
      RETURN
C
      END
C***********************************************************************
      SUBROUTINE CALCBF(BJ,BN,NL,X)
      DIMENSION BJ(NL),BN(NL)
      IF(ABS(X).LT.1.0E-6)GOTO 7
      BJ(1)=SIN(X)/X
      BN(1)=-COS(X)/X
      IF(NL.EQ.1)RETURN
      BJ(2)=(BJ(1)-COS(X))/X
      BN(2)=(BN(1)-SIN(X))/X
      IF(NL.EQ.2)RETURN
      IF(FLOAT(NL*(NL+1)).GT.X*X)GOTO 2
C  FORWARD RECURRENCE FOR BJ'S
      FL=3.0
      DO 1 L=3,NL
      BJ(L)=FL*BJ(L-1)/X-BJ(L-2)
1     FL=FL+2.
      GOTO 5
C  BACKWARD RECURRENCE FOR BJ'S
2     BJ0=BJ(1)
      BJ1=BJ(2)
      NN=MAX0(10,2*NL)
      A=0.
      B=1.
      FL=FLOAT(2*NN+1)
      DO 3 I=1,NN
      L=NN-I+1
      C=FL*B/X-A
      IF(L.LE.NL)BJ(L)=C
      A=B
      B=C
3     FL=FL-2.
C  NORMALISATION
      B=BJ0/BJ(1)
      IF(ABS(BJ0).LT.0.01)B=BJ1/BJ(2)
      DO 4 L=1,NL
4     BJ(L)=B*BJ(L)
C  FORWARD RECURRENCE FOR BN'S
5     FL=3.
      DO 6 L=3,NL
      BN(L)=FL*BN(L-1)/X-BN(L-2)
6     FL=FL+2.
      RETURN
7     WRITE(6,200)X
      RETURN
200   FORMAT(13H0** ARGUMENT?,E12.4,29H TOO SMALL FOR ROUTINE CALCBF)
      END
C
C---------------------------------------------------------------------
C  subroutine PHSH_WIL
C---------------------------------------------------------------------
C  A.R. WILLIAMS^ PHASE SHIFT PROGRAM (GIVEN A MUFFIN-TIN POTENTIAL)
      SUBROUTINE PHSH_WIL(MUFFTIN_FILE, PHASOUT_FILE, DATAPH_FILE)
      CHARACTER(LEN=*), INTENT(IN)	:: MUFFTIN_FILE
      CHARACTER(LEN=*), INTENT(IN)	:: PHASOUT_FILE, DATAPH_FILE
      REAL E(401),S(401,15),C(401,15),DEL(15),DELOLD(15)
      REAL DELL(9),delstore(8,401,15)
      INTEGER TLP1
      COMMON / CM16 / E1, E2, NE, IX,NEUO
      COMMON / CMRV / R(201), V(201, 15), NR, NL, Z
      COMMON / CM5 / Y(30,4), F(30,4), ILST
      NAMELIST / NL2 / IP,NRR
C
C First input channels
C
      OPEN (UNIT=5,FILE=MUFFTIN_FILE,STATUS='OLD')
C
C Now output channels
C
      OPEN (UNIT=6,FILE='zph.o',STATUS='UNKNOWN')
      OPEN (UNIT=7,FILE=PHASOUT_FILE,STATUS='UNKNOWN')
      OPEN (UNIT=8,FILE=DATAPH_FILE,STATUS='UNKNOWN')
C
      PI=3.1415926535
C  READ IP
C   IP=0: ONLY RADIAL WAVEFUNCTION
C   IP=1: PHASE SHIFTS IN ADDITION
C   IP=2: S AND C
C   IP=3: PRODUCE LOGARITHM OF PHASE SHIFTS
C   NRR= number of inequivalent atoms for which we want phase shifts
      IP=1
      READ(5, NL2)
      WRITE(6, NL2)
      write(8,110)
110   format (11HTitleText: ,8HDELTA(E))
C  INPUT
      DO 2  KKK=1,NRR
      CALL S16
      TX = 2. * R(NR)/ FLOAT(NR - 1)
      DE = (E2 - E1) / FLOAT(MAX0(NE - 1, 1))
      DO 6 I = 1, NE
      E(I) = E1 + FLOAT(I - 1) * DE
C  RADIAL INTEGRATION
      CALL S10(E(I))
      T3 = R(NR) * E(I)
      T4 = R(NR) * T3
      DO 6 LP1 = 1, NL
      L = LP1 - 1
      TLP1 = 2 * L + 1
      T5 = R(NR) ** LP1
      UT = F(TLP1, ILST)/TX + FLOAT(L) * Y(TLP1,ILST) / R(NR)
      T1 = (F44(L,T4) * Y(2*LP1,ILST) + T3 * F44(LP1,T4) * Y(TLP1,ILST))
     1 * T5
      T2 = (F45(L,T4) * UT - T3 * F45(L - 1,T4) * Y(TLP1,ILST)) * R(NR)/
     1 T5
      S(I, LP1) = T1
6     C(I, LP1) = T2
      IS = 2
      I4 = 9
      IF(IP .LT. 1) GO TO 15
C  PRODUCE PHASE SHIFTS
      DO 8 LP=1,NL
8     DELOLD(LP)=0.0
      DO 10 I = 1, NE
      DO 11 LP = 1, NL
11    DEL(LP) = ATAN(-ABS(E(I)) ** (LP - .5) * S(I,LP) / C(I, LP))
C  REMOVE DISCONTINUITIES BY MULTIPLES OF PI
      DO 117 LP=1,NL
      LS=0
111   DELDIF=DEL(LP)-DELOLD(LP)
      IF (ABS(DELDIF).LT.0.7) GO TO 117
      LS=LS+1
      DEL(LP)=DEL(LP)-SIGN(PI,DELDIF)
      IF (LS.LT.5) GO TO 111
      WRITE (6,115) LP
115   FORMAT(36H TOO LARGE CHANGE IN PHASE SHIFT [L=,1I4,
     %20H] SINCE LAST ENERGY ,/,41H DISCONTINUITY BY MULTIPLE OF PI POSS
     %IBLE)
117   DELOLD(LP)=DEL(LP)
      IF (NEUO.EQ.2) E(I)=0.5*E(I)
C  PRINT PHASE SHIFTS
      WRITE(6, 12) E(I), (DEL(LP), LP = 1, NL)
12    FORMAT(1P8E14.7, /, 14X, 1P7E14.7, /)
C  PUNCH PHASE SHIFTS IN FORMAT USED BY LEED PROGRAM
C      WRITE(7,71) E(I),(DEL (LP),LP=1,NL)
71    FORMAT(1F7.4)
72    FORMAT(10F7.4)
C  store phase shifts
      do 144 kk=1,nl
144   delstore(KKK,I,kk)=del(kk)
      IF(IP .LT. 3) GO TO 10
      DO 14 J = 1, 9
      DELL(J) = -4.
      IF(DEL(J) .LT. 1.0E-4) GO TO 14
      DELL(J) = ALOG10(DEL(J))
14    CONTINUE
10    CONTINUE
C  write phase shifts as function of energy for plotting
      do 145 kk=1,nl
         write(8,100) kk-1
         do 146 ii=1,NE
146           write(8,*) E(II),delstore(kkk,II,kk)
145        write(8,*)
100      format(3H"L=,i2)
15    CONTINUE
2     CONTINUE
C      IF(IP .LT. 2) GO TO 2
       WRITE(7,*) 'BE CAREFUL ABOUT THE ORDER OF THE ELEMENTS'
       do 148 ii=1,NE
        WRITE(7,71) E(II)
        do 147 i=1,nrr
          WRITE(7,72)(delstore(i,ii,LP),LP=1,NL)
147     continue
148   continue
      if(IP .LT. 2) goto 22
      DO 9 LP1 = 1, NL
      CALL S41(E, S(1, LP1), NE)
9     CALL S41(E, C(1, LP1), NE)
C      GO TO 2
22    continue
      RETURN
      END
C***********************************************************************
C  SUBROUTINE S16
C  S16 INPUTS DATA
C   CS: CORE SHIFT (POSITION OF ZERO OF ENERGY)
C   Z: ATOMIC NUMBER
C   E1,E2: FIRST AND LAST ENERGIES DESIRED (IN RYDBERGS OR
C           HARTREES, CF. NEUI)
C   NE: NUMBER OF ENERGIES DESIRED
C   NL: NUMBER OF PHASE SHIFTS DESIRED (=LMAX+1)
C   NR: NUMBER OF RADIAL GRID POINTS USED IN CALCULATION
C   IX=0: SIGNAL TO STOP
C   IX=1: SIGNAL TO EXPECT A (NEW) POTENTIAL TO BE READ IN
C   RT: MUFFIN-TIN RADIUS (IN BOHR RADII)
C   NEUI,NEUO
C    IF =1: RYDBERG UNIT USED FOR INPUT (NEUI) AND OUTPUT (NEUO) OF
C           ENERGIES AND POTENTIAL
C    IF =2: HARTREE UNIT (DOUBLE RYDBERG) USED INSTEAD OF RYDBERG
C           UNIT FOR INPUT (NEUI) AND OUTPUT (NEUO)
C   POTYP=1: RADIAL INPUT AS V(R)
C   POTYP=2: RADIAL INPUT AS R*V(R)
C
C***********************************************************************
      SUBROUTINE S16
      COMMON / CM16 / E1, E2, NE, IX,NEUO
      COMMON / CMRV / R, V, NR, NL, Z
      DIMENSION R(201), V(201, 15)
      REAL RS(200),          ZS(200), ZTT(201)
      DIMENSION FMT(18)
      NAMELIST / NL16 / CS,Z,E1,E2,NE,NL,NR,IX,RT,NEUI,NEUO,POTYP
C  SET DEFAULT VALUES OF VARIABLES IN NAMELIST /NL16/
      IX=1
      E1 = 4.
      E2 = 24.0
      NE = 30
      NL = 9
      NR = 101
      NEUI=1
      NEUO=2
      POTYP=2
      READ (5, NL16)
      CS=0.0
      IF(IX .LT. 1) STOP
      IF (NEUI.EQ.1) GO TO 5
      CS=2.0*CS
      E1=2.0*E1
      E2=2.0*E2
5     CONTINUE
      WRITE (6, NL16)
      DRDN2 = (FLOAT(NR - 1))** 2 / RT
C  READ FORMAT USED FOR INPUT OF R VS. V(R) OR R VS. R*V(R)
C  (V IS ASSUMED POSITIVE)
C      READ(5, 8) (FMT(I), I = 1, 18)
8     FORMAT(18A4)
111   FORMAT(2E14.5)
      DO 16 I=1,200
      READ(5,111) RS(I), ZS(I)
C the next lines assume that the input potential and cs are negative
      ZS(I)=-ZS(I)
      IF( RS(I) .LT. 0) GO TO 17
16    CONTINUE
17    NRS = I - 1
      IF (NEUI.EQ.1) GO TO 174
      DO 172 I=1,NRS
172   ZS(I)=2.0*ZS(I)
174   IF (POTYP.EQ.2) GO TO 178
      DO 176 I=1,NRS
176   ZS(I)=(ZS(I)-CS)*RS(I)
      GO TO 21
178   CONTINUE
      DO 20 I = 2, NRS
20    ZS(I) = (ZS(I) / RS(I) - CS) * RS(I)
21    IV = 1
      R(1) = 0.
      ZTT(1) = Z + Z
      DO 1 I = 2, NR
      R(I) = (FLOAT(I - 1)) ** 2 / DRDN2
40    IF(R(I) .LE. RS(IV + 2)) GO TO 50
      IF ( IV + 3 .GE. NRS) GO TO 50
      IV = IV + 1
      GO TO 40
50    ZTT(I) = F12(RS(IV), ZS(IV), R(I), 4)
2     DO 1 LP1 = 1, NL
1     V(I, LP1) = -ZTT(I) / R(I)
      RETURN
      END
C***********************************************************************
C  F12 PERFORMS ITERATIVE INTERPOLATION IN A TABLE OF N VALUES OF
C  X AND Y TO FIND THE VALUE OF Y AT Z
C***********************************************************************
      FUNCTION F12(X, Y, Z, N)
      REAL X(10), Y(10), W(20)
      W(1) = Y(1)
      DO 1 I = 2, N
      W(I) = Y(I)
      U = Z - X(I)
      IP1 = I + 1
      DO 1 J = 2, I
      K = IP1 - J
1     W(K) = W(K + 1) + U * (W(K) - W(K + 1)) / (X(K) - X(I))
      F12 = W(1)
      RETURN
      END
C***********************************************************************
C  S5 -- HAMMING^S METHOD FOR THE INTEGRATION OF SYSTEMS OF FIRST
C  ORDER DIFFERENTIAL EQUATIONS
C***********************************************************************
      SUBROUTINE S5(E)
      REAL EEST(30), VME(15)
      COMMON / CMRV / R(201), V(201, 15), NR, NL, Z
      COMMON / CM5 / Y(30, 4), F(30, 4), IP1
      NJ = 2 * NL
      DO 5 J = 1, NJ
5     EEST(J) = 0.
      DO 4 I = 5, NR
      DO 6 LP1 = 1, NL
6     VME(LP1) = (V(I, LP1) - E) * R(I)
      T1 = 2. / FLOAT(I - 1)
      IP1 = MOD(I - 1, 4) + 1
      IM2 = MOD(IP1, 4) + 1
      IM1 = MOD(IM2, 4) + 1
      IP0 = MOD(IM1, 4) + 1
      DO 1 J = 1, NJ
      F(J, IM2) = Y(J, IP1) + (2. * (F(J, IP0) + F(J, IM2)) -
     1 F(J, IM1)) / 0.75
1     Y(J, IP1) = F(J, IM2) - 0.925619835 * EEST(J)
      DO 2 J = 1, NJ, 2
      JP1 = J + 1
      LP1 = JP1 / 2
      FLP1 = LP1
      F(J, IP1) = (FLP1 * Y(J, IP1) + R(I) * Y(JP1, IP1)) * T1
2     F(JP1, IP1) = (VME(LP1) * Y(J, IP1) - FLP1 * Y(JP1, IP1)) * T1
      DO 3 J = 1, NJ
      Y(J, IP1) = Y(J, IP0) + (Y(J, IP0) - Y(J, IM2) + 3. * (F(J, IP1)
     1 + 2. * F(J, IP0) - F(J, IM1))) / 8.
      EEST(J) = F(J, IM2) - Y(J, IP1)
3     Y(J, IP1) = Y(J, IP1) + .743801653E-1 * EEST(J)
      DO 4 J = 1, NJ, 2
      JP1 = J + 1
      LP1 = JP1 / 2
      FLP1 = LP1
      F(J, IP1) = (FLP1 * Y(J, IP1) + R(I) * Y(JP1, IP1)) * T1
4     F(JP1, IP1) = (VME(LP1) * Y(J, IP1) - FLP1 * Y(JP1, IP1)) * T1
      RETURN
      END
C***********************************************************************
C  S10   POWER SERIES EXPANSION OF THE SOLUTION ABOUT THE ORIGIN
C        AND RADIAL INTEGRATION IN S5
C***********************************************************************
      SUBROUTINE S10(E)
      INTEGER TLP1
      REAL A(10), B(10), TR(4)
      COMMON / CMRV / R(201), V(201, 15), NR, NL, Z
      COMMON / CM5 / Y(30, 4) , F(30, 4) , ILST
      NI = 2 * NL
      TZ = 2. * Z
      A(1) = 1.
      DO 5 I = 1, NI, 2
      LP1 = (I + 1) / 2
      TLP1 = 2 * LP1
      EP = E - V(4, LP1) - TZ/ R(4)
      Y(I, 1) = 0.
      Y(I + 1, 1) = 0.
      A(1) = A(1) / FLOAT(2 * LP1 - 1)
      B(1) = - Z * A(1) / FLOAT(LP1)
      DO 1 J = 2, 4
      TR(J) = R(J) ** LP1
      Y(I, J) = A(1) * TR(J)
1     Y(I + 1, J) = B(1) * TR(J)
      DO 3 K = 1, 9
      A(K + 1) = B(K) / FLOAT(K)
      B(K + 1) = -(EP * A(K) + TZ * A(K + 1)) / FLOAT(TLP1 + K)
      DO 2 J = 2, 4
      TR(J) = TR(J) * R(J)
      Y(I, J) = Y(I, J) + TR(J) * A(K + 1)
2     Y(I + 1, J) = Y(I + 1, J) + TR(J) * B(K + 1)
      IF(ABS(TR(4) * A(K + 1) / Y(I, 4)) .LT. 1.0E-4) GO TO 5
3     CONTINUE
      WRITE (6, 4) E, LP1, R(4), (A(K), K = 1, 10)
4     FORMAT(1PE10.2, I10, 11E10.2)
5     CONTINUE
      DO 6 J = 2, 4
      T1 = 2. / FLOAT(J - 1)
      DO 6 I = 1, NI, 2
      IP1 = I + 1
      LP1 = IP1 / 2
      FLP1 = LP1
      F(I, J) = (FLP1 * Y(I, J) + R(J) * Y(IP1, J)) * T1
6     F(IP1,J) = ((V(J, LP1) - E) * R(J) * Y(I, J) - FLP1 * Y(IP1, J))
     1  * T1
      CALL S5(E)
      RETURN
      END
C***********************************************************************
C  F44  EVALUATES THE SPECIAL VERSION OF THE SPHERICAL BESSEL FUNCT.
C***********************************************************************
      FUNCTION F44(L, X)
      REAL S(20)
      JS = L + L + 1
      IF(ABS(X / FLOAT(JS)) .GT. 10.) GO TO 5
      FI = 1.
      IF(L .LT. 1) GO TO 2
      DO 1 K = 3, JS, 2
1     FI = FI * FLOAT(K)
2     T1 = 1. / FI
      DT = 1.
      T = 1.
      I = 0
      DO 3 K = 1, 100
      I = I + 2
      DT = -DT * X / FLOAT(I * (I + JS))
      T = T + DT
      IF(ABS(DT).LT. 1.E-8) GO TO 4
3     CONTINUE
4     T1 = T1 * T
      F44 = T1
      RETURN
5     T = SQRT(ABS(X))
      IF(X .LT. 0.) GO TO 9
      S(2) = SIN(T) / T
      IF(L .GT. 0) GO TO 6
11    F44 = S(2)
      RETURN
6     S(1) = COS(T)
      GO TO 10
9     S(2) = SINH(T) / T
      IF(L .LT. 1) GO TO 11
      S(1) = COSH(T)
10    IS = L + 2
      DO 7 I = 3, IS
7     S(I) = (S(I - 1) * FLOAT(2 * I - 5) - S(I - 2)) / X
      F44 = S(IS)
      RETURN
      END
C***********************************************************************
C  F45  EVALUATES SPECIAL VERSION OF THE SPHERICAL NEUMANN FUNCTION
C***********************************************************************
      FUNCTION F45(L, X)
      REAL S(20)
      IF(L .LT. 0) GO TO 8
      LP1 = L + 1
      JS = L + L + 1
      IF(ABS(X / FLOAT(JS)) .GT. 10.) GO TO 5
      FI = 1.
      IF(L .LT. 1) GO TO 2
      DO 1 K = 3, JS, 2
1     FI = FI * FLOAT(K)
2     T1 = FI / FLOAT(JS)
      DT = 1.
      T = 1.
      I = 0
      DO 3 K = 1, 100
      I = I + 2
      DT = -DT * X / FLOAT(I * (I - JS))
      T = T + DT
      IF(ABS(DT) .LT. 1.E-8) GO TO 4
3     CONTINUE
4     T1 = T1 * T
      F45 = T1
      RETURN
5     T = SQRT(ABS(X))
      IF(X .LT. 0.) GO TO 9
      S(2) = COS(T)
      IF(L .GT. 0) GO TO 6
11    F45 = S(2)
      RETURN
6     S(1) = -SIN(T) / T
      GO TO 10
9     S(2) = COSH(T)
      IF(L .LT. 1) GO TO 11
      S(1) = -SINH(T) / T
10    IS = L + 2
      DO 7 I = 3, IS
7     S(I) = S(I - 1) * FLOAT(2 * I - 5) - X * S(I - 2)
      F45 = S(IS)
      RETURN
8     F45 = - F44(L+1, X)
      RETURN
      END
C***********************************************************************
C  S41 PLOTS Y AGAINST X
C***********************************************************************
      SUBROUTINE S41(X, Y, N)
      DIMENSION X(100), Y(100), P(97)
      DATA B, C, O, D / 1H , 1H*, 1H0, 1HI /
      Y1 = 0
      Y2 = 0
      DO 1 I = 1, N
      Y1 = AMIN1(Y1, Y(I))
1     Y2 = AMAX1(Y2, Y(I))
      DO 2 I = 1, 97
2     P(I) = B
      T = 96/ (Y2 - Y1)
      J0 = -Y1 * T + 1.5
      P(J0) = 0
      IF(N .GE. 30) WRITE(6, 3) P
3     FORMAT(1H1, 34X, 97A1, //)
      IF(N .LT. 30) WRITE(6, 6) P
6     FORMAT(////, 35X, 97A1, //)
      P(J0) = D
      DO 5 I = 1, N
      J = T * (Y(I) - Y1) + 1.5
      P(J) = C
      WRITE(6, 4) X(I), Y(I), P
4     FORMAT(1H , 1P2E16.6, 2X, 97A1)
      P(J) = B
5     P(J0) = D
      RETURN
      END
C
C---------------------------------------------------------------------
C  subroutine PHSH_REL
C---------------------------------------------------------------------
      SUBROUTINE PHSH_REL(MUFFTIN_FILE, PHASOUT_FILE, DATAPH_FILE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER(LEN=*), INTENT(IN)	:: MUFFTIN_FILE
      CHARACTER(LEN=*), INTENT(IN)	:: PHASOUT_FILE, DATAPH_FILE
      CHARACTER OPT*3,OPT1*3,OPTS*3,ANAME*2,AN*30,BDATA*28
      CHARACTER SUB*3,RECORD*3,TL*1,SL*1,SS1*6,SS2*6,WRD*6
      CHARACTER AMS(5)*4
      REAL JF,JFS,JFD,JFU
      REAL NAME(4)
      DIMENSION JF(250,18),ENERG (250)
C------>PROGRAM ATMSET (INPUT,OUTPUT,PUNCH)
C     COMPILATION                           SEPT. 13, 1973        SWWLIU
      COMMON/ZZZZ/ZP(340),VS,IPT,JRI
      COMMON /Z/ RMAXI
      DIMENSION ADATA(7)
      DATA ZERO/0.0D0/,ONE/1.0D0/,TWO/2.0D0/,ANINE/9.0D0/,HALF/.5D0/
      DATA ZILCH/1.0D-4/,TOL/.005D+0/,DES/.025D+0/
      DATA AMS/'NC= ','L= ',' ES=',' DE=','ID= '/
      DATA TL/'L'/,SL/'S'/,SS1/'NOSPIN'/,SS2/' SPIN '/
      DATA SUB/'SUB'/,RECORD/'NOS'/
c    1 FORMAT (3D12.4,T41,2A3,I3)
    1 FORMAT (3D12.4,4X,I3,4X,D12.4)
    2 FORMAT (5E14.6)
 3    FORMAT    (///,10X,A6,////,10X,18HMUFFIN-TIN RADIUS=,
     1F10.6,4X,14HCONSTANT POT.=,F10.6,//,10X,'ATOMIC DATA SET FOR Z=',
     2I3,' AND L=',I3,4X,'     ',A2,//,14X,5HE(EV),12X,7HL - 0.5,13X,
     37HL + 0.5,12X,10HS-AVERAGED,/)
  4   FORMAT (10X,F10.6,3F20.8)
   8  FORMAT (F9.4,8F8.5)
  9   FORMAT (1X,A1,I3,3D15.7,22X,A3)
   12 FORMAT (5D14.6)
   13 FORMAT (1H1,//,T61,'INPUT DATA',//)
   15 FORMAT (A28,T29,A2,T35,A30,T65,F10.7,T76,A3)
   17 FORMAT (1H1,/,10X,'RELATIVISTIC PHASE SHIFTS FOR ',A30,//,10X,'EX
     1CA =',F10.6,4X,'EXCB =',F10.6,4X,'EXCO =',F10.6//,10X,'LATTICE CON
     2STANT =',F10.6,' ,',F10.6,' ,',F10.6)
   18 FORMAT(F10.4,F9.4,2I5)
      OPEN(UNIT=4,FILE='inpdat',STATUS='unknown')
      OPEN(UNIT=5,FILE=MUFFTIN_FILE,STATUS='OLD')
C      OPEN(UNIT=6,FILE='prot',STATUS='unknown')
      OPEN(UNIT=7,FILE=PHASOUT_FILE,STATUS='unknown')
      OPEN (UNIT=8,FILE=DATAPH_FILE,STATUS='UNKNOWN')
C
C
      PI=4.D0*DATAN(1.D0)
      PI2=0.5D0*PI
      WRITE (4,13)
10    READ(5,217,END=999,ERR=999) (NAME(I),I=1,4)
217   FORMAT(4A4)
c      READ(5,1)ES,DE,UE,OPT,OPT1,LSM                                    C
      READ(5,1)ES,DE,UE,LSM,VC
c nl is the number of plotted phase shifts
      nl=8
      WRITE (4,11) ES,DE,UE,OPT,OPT1,LSM
   11 FORMAT (3D12.4,4X,2A3,I3)
c  10  READ (5,15,END=999,ERR=999) BDATA,ANAME,    AN   ,VC,OPTS
c      WRITE (4,75) BDATA,ANAME,AN,VC,OPTS
   75 FORMAT (A28,A2,4X,A30,F10.7,1X,A3)
      READ (5,16) NZ,ADATA(1),JRI,ALC,BLC,CLC,EXCA,EXCB,EXCO
   16 FORMAT (I4,F10.6,I4,T21,6F10.6)
      WRITE (4,76) NZ,ADATA(1),JRI,ALC,BLC,CLC,EXCA,EXCB,EXCO
   76 FORMAT (I4,F10.6,I4,2X,6F10.6)
      VS=0.
      IF(OPTS.EQ.SUB) VS=VC
      IF ((JRI.LE.0).OR.(JRI.GT.340))  GO TO 999
      READ(5,2) (ZP(J),J=1,JRI)                                         C
      WRITE (4,12) (ZP(J),J=1,JRI)
      RHOZ=-0.90306514D+01
      DELRHO=0.3125D-01
      RM=DEXP(RHOZ)
      XRX=DEXP(DELRHO)
      DO 19 J=1,JRI
      IF(ZP(J).LT.0.0) ZP(J)=-ZP(J)
      IF(J.EQ.JRI) GO TO 21
 19   RM=XRX*RM
 21   CONTINUE
      IF (DE.GT.ZERO)  GO TO 20
      ES=-HALF
      DE=DES
      UE=ONE
  20  N=(UE-ES)/DE+HALF
      N=N+1
      WRITE(07,181)(NAME(I),I=1,4)
C     WRITE(15,181) AN
  181 FORMAT('RELATIVISTIC PHASE SHIFTS FOR ',4A4)
      WRITE(07,18) ES,DE,N,LSM
C     WRITE(15,18) ES,DE,N,LSM
      ES=ES/13.6
      DE=DE/13.6
      UE=UE/13.6
      IF (N.GT.250)  N=250
      L=0
      E=ES
      IPT=2
      IF(OPT.EQ.RECORD) GO TO 23
      WRD=SS2
      GO TO 24
 23   IPT=-2
      WRD=SS1
 24   CONTINUE
C      WRITE (6,17) AN,EXCA,EXCB,EXCO,ALC,BLC,CLC
C      WRITE (6,3) WRD,RM,VC,NZ,L,ANAME
      KAP=-1
      L=1
      DO 30 J=1,N
      DXAZ=0.0D0
      TTR=DLGKAP(E,KAP)/(12.5663706)
      CALL SBFIT(TTR,E,L-1,RMAXI,JFS)
      JF(J,L)=JFS
      E=E*13.6
      ENERG(J)=E
C      WRITE(6,4) E,JFS,JFS,JFS
      E=E/13.6
      E=E+DE
  30  CONTINUE
   40 IF(L.GT.LSM) GO TO 80
C      WRITE (6,17) AN,EXCA,EXCB,EXCO,ALC,BLC,CLC
C      WRITE (6,3) WRD,RM,VC,NZ,L,ANAME
      KAP=-(L+1)
      LIND=L+1
      E=ES
      LVOR=0
      DO 50 J=1,N
      DLU=DLGKAP (E,KAP)/(12.5663706)
      DLD=DLGKAP(E,L)/(12.5663706)
      CALL SBFIT(DLD,E,L,RMAXI,JFD)
      CALL SBFIT(DLU,E,L,RMAXI,JFU)
      LK=0
C      JFDIFF=-(JFD-JFU)*LVOR
      ZFDIFF=-(JFD-JFU)*LVOR
      IF (ZFDIFF .GT. 2.5 )LK=L
      IF (ZFDIFF .LT. -2.5 )LK=L+1
      JFS=L*JFD-KAP*JFU+LVOR*LK*PI
      JFS=JFS/(2*L+1)
      IF (JFS .GT.PI2 ) JFS=JFS-PI2*2.
      IF (JFS .LT. -PI2 ) JFS=JFS+PI2*2.
      JF(J,LIND)=JFS
      IF (LK .EQ. 0) LVOR=SIGN(1.,JFS)
      E=E*13.6
C      WRITE(6,4) E,JFD,JFU,JFS
      E=E/13.6
      E=E+DE
  50  CONTINUE
      L=L+1
      GO TO 40
   80 DO 90 I=1,N
      LSM1=LSM+1
      WRITE(7,8) ENERG(I),(JF(I,L),L=1,LSM1)
C     WRITE(15,8) (JF(I,L),L=1,LSM1)
   90 CONTINUE
      do 145 kk=1,nl
         write(8,100) kk-1
         do 146 ii=1,N
146           write(8,*) ENERG(II),JF(II,kk)
145        write(8,*)
100      format(3H"L=,i2)
      ES=ES*13.6
      DE=DE*13.6
      UE=UE*13.6
      GO TO 10
  999 CONTINUE
      WRITE (4,900)
  900 FORMAT (//,T57,'END OF INPUT DATA')
      STOP
      END
C***********************************************************************
C...  DLGKAP CALCULATES THE LOGRITHMIC DERIVATIVE OF THE LARGE
C...  COMPONENT USING THE PROCEDURE DESCRIBED BY LOUCKS IN APPENDIX 7.
C+++  THE SMALL MULTIPLICATIVE FACTOR IS INCLUDED.
C...  POTENTIAL  IN THE FORM OF 2ZP IS TO BE PASSED IN THE COMMON /ZZZZ/
C...  THE RADIAL FUNCTIONS ARE MADE AVAILABLE IN THE COMMON /RADFUN/
C***  WABER MESH (XNOT=-9.03065 AND DX=1/32) IS USED
C***  JRI IS THE MESH POINT OF THE APW SPHERE RADIUS
C...  E IS THE ENERGY TO BE USED (IN RYDBERGS)
C***** 4 PI R**2 INSERTED NOW.   FOR COMPOUNDS ONLY.
C***********************************************************************
      FUNCTION DLGKAP (E,KAPPA)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION POT(340),U(340),W(340),UP(340),WP(340),SXK(4),SXM(4)
      COMMON /Z/ T                                                      CCCCCCC
C
      COMMON /RADFUN/ U,W
      COMMON/ZZZZ/POT,VCZ,IPT,JRI
      DATA USTART/1.D-25/,ZILCH/1.D-30/,TEST/1.D+6/,XS/-9.03065133D+00/,
     1 DX/3.125D-2/,C/2.740746D+2/,CIN/1.3312581146D-5/,HF/.5D+0/,
     2 TH/.3333333333D+0/,T2/2.D+0/,T7/7.D+0/,T11/11.D+0/,T12/12.D+0/,
     3 T14/14.D+0/,T26/26.D+0/,T32/32.D+0/,ZERO/.1D+0/
   83 FORMAT (10H HARD TEST,D20.8, I5, 4D20.8 )
C**********SET UP FOR RELATIVISTIC OR NO RELATIVISTIC EFFECT************
      IF(IPT.GT.0) GO TO 88
      CIN=0.0D00
 88   CONTINUE
C******************* SET UP STARTING VALUES  ***************************
      DX2=HF*DX
      X20=(.3)*DX
      XMFT=(4.4444444444444D-2)*DX
      TS=DEXP(XS)
      TDX=DEXP(DX)
      HOC=(VCZ*TS+POT(1))/C
      XK=KAPPA
   49 U(1)=USTART
      IF (DABS(HOC/XK).GT. .05) GO TO 6
      P=(XK+DABS(XK))/HOC - HF*HOC/DABS(XK)
      GO TO 7
    6 P=(XK+DSQRT(XK*XK-HOC*HOC))/HOC
 7    TC=(E+VCZ)*TS+POT(1)
      VC=CIN*TC
      W(1)=C*P*USTART
C***************** START RUNGE-KUTTE PROCEDURE  ************************
   11 X = XS
      N=1
   25 IK= 0
      NP1=N+1
      XC=X
      BGC=POT(N)
      WC= W(N)
      UC= U(N)
  20  IK=IK+1
      T=DEXP(XC)
      TC=(E+VCZ)*T+BGC
      VC=CIN*TC
  12  SXK(IK)=DX2*(WC*(VC+T)-XK*UC)
      SXM(IK)=DX2*(XK*WC-TC*UC)
   15 GO TO (16,17,18,19),IK
  16  XC=XC+DX2
      UC=UC+SXK(1)
      WC=WC+SXM(1)
      BGC=HF*(BGC+POT(NP1))
      GO TO 20
   17 UC=UC+SXK(2)-SXK(1)
      WC=WC+SXM(2)-SXM(1)
      GO TO 20
  18  XC=XC+DX2
      UC=UC+T2*SXK(3)-SXK(2)
      WC=WC+T2*SXM(3)-SXM(2)
      BGC=POT(NP1)
      GO TO 20
  19  W(NP1)=W(N)+(SXM(1)+SXM(4)+T2*(SXM(2)+SXM(3)))*TH
      U(NP1)=U(N)+(SXK(1)+SXK(4)+T2*(SXK(2)+SXK(3)))*TH
      UP(NP1)=(VC+T)*W(NP1)-XK*U(NP1)
      WP(NP1)=XK*W(NP1)-TC*U(NP1)
  24  X=X+DX
      N=NP1
      IF (N.LT.6) GO TO 25
C***********************************************************************
C     END OF STARTING INTEGRATION.  BEGIN MILNE PROCEDURE.
C***********************************************************************
      T=DEXP(X)
   26 T=T*TDX
      NP1=N+1
      NM1=N-1
      NM2=N-2
      NM3=N-3
      NM4=N-4
      NM5=N-5
      TC=(E+VCZ)*T+POT(NP1)
      VC=CIN*TC
  27  UNP=U(NM5)+X20*(T11*(UP(N)+UP(NM4))+T26*UP(NM2)                   -
     C                                 -T14*(UP(NM1)+UP(NM3)))
      WNP=W(NM5)+X20*(T11*(WP(N)+WP(NM4))+T26*WP(NM2)                   -
     C                                      -T14*(WP(NM1)+WP(NM3)))
      NIT=0
  33  UP(NP1)=(VC+T)*WNP-XK*UNP
      WP(NP1)=XK*WNP-TC*UNP
      UNP2=U(NM3)+(T7*(UP(NP1)+UP(NM3))+T32*(UP(NM2)+UP(N))+T12*UP(NM1))-
     C*XMFT
      WNP2=W(NM3)+(T7*(WP(NP1)+WP(NM3))+T32*(WP(NM2)+WP(N))+T12*WP(NM1))-
     C*XMFT
C********  COMPARE PREDICTOR WITH CORRECTOR  ***************************
      IF (DABS(TEST*(UNP2 -UNP)).GT.DABS(UNP2)) GO TO 31
      IF (DABS(TEST*(WNP2-WNP)).LE.DABS(WNP2)) GO TO 32
   31 IF (NIT.LT.5) GO TO 81
C   82 WRITE(6,83)   E,N,UNP2,UNP,WNP2,WNP
      GO TO 32
  81  NIT=NIT+1
      WNP=WNP2
      UNP=UNP2
      GO TO 33
  32  W(NP1)=WNP2
      U(NP1)=UNP2
      N=NP1
      IF (N.LT.JRI) GO TO 26
C******************  END OF MILNE PROCEDURE  ***************************
      IF ( DABS(U(JRI)).GT.ZILCH)  GO TO 37
      U(JRI)=DSIGN(ZILCH,U(JRI))
  37  P=(T+VC)/T
      WNP=P*W(JRI)/U(JRI)
      UNP=WNP-(KAPPA+1)/T
          DLGKAP=(12.5663706)*T*T*UNP
   46 RETURN
      END
C***********************************************************************
      SUBROUTINE SBFIT(T,E,L,R,JFS)
      IMPLICIT DOUBLE PRECISION(E,R,T)
      REAL JFS,KAPPA
      SE=SNGL(E)
      SR=SNGL(R)
      ST=SNGL(T)
      KAPPA=SQRT(SE)
      X=KAPPA*SR
      BJ1=SIN(X)/X
      BN1=-COS(X)/X
      BJ2=BJ1/X + BN1
      BN2=BN1/X - BJ1
      IF(L) 3,3,1
    1 LS=1
    2 LS=LS+1
      BJT=(2*LS-1)*BJ2/X - BJ1
      BNT=(2*LS-1)*BN2/X - BN1
      BJ1=BJ2
      BJ2=BJT
      BN1=BN2
      BN2=BNT
      IF(L+1-LS) 3,3,2
    3 DL=ST / (SR*SR)
      DL=DL - L/SR
      AN=DL*BJ1 + KAPPA*BJ2
      AD=DL*BN1 + KAPPA*BN2
      JFS=3.141592654/2.0
      IF(ABS(AD)-1.0E-8) 5,5,4
    4 JFS=ATAN(AN/AD)
    5 RETURN
      END
