      SUBROUTINE   DNSF(N, P, L, ALF, C, Y, CALCA, INC, IINC, IV,
     1                  LIV, LV, V, UIPARM, URPARM, UFPARM)
C
C  ***  SOLVE SEPARABLE NONLINEAR LEAST SQUARES USING
C  ***  FINITE-DIFFERENCE DERIVATIVES.
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER IINC, L, LIV, LV, N, P
C/6
C     INTEGER INC(IINC,P), IV(LIV), UIPARM(1)
C     DOUBLE PRECISION ALF(P), C(L), URPARM(1), V(LV), Y(N)
C/7
      INTEGER INC(IINC,P), IV(LIV), UIPARM(*)
      DOUBLE PRECISION ALF(P), C(L), URPARM(*), V(LV), Y(N)
C/
      EXTERNAL CALCA, UFPARM
C
C  ***  PARAMETERS  ***
C
C      N (IN)  NUMBER OF OBSERVATIONS.
C      P (IN)  NUMBER OF NONLINEAR PARAMETERS TO BE ESTIMATED.
C      L (IN)  NUMBER OF LINEAR PARAMETERS TO BE ESTIMATED.
C    ALF (I/O) NONLINEAR PARAMETERS.
C                 INPUT = INITIAL GUESS,
C                 OUTPUT = BEST ESTIMATE FOUND.
C      C (OUT) LINEAR PARAMETERS (ESTIMATED).
C      Y (IN)  RIGHT-HAND SIDE VECTOR.
C  CALCA (IN)  SUBROUTINE TO COMPUTE A MATRIX.
C    INC (IN)  INCIDENCE MATRIX OF DEPENDENCIES OF COLUMNS OF A ON
C                 COMPONENTS OF ALF -- INC(I,J) = 1 MEANS COLUMN I
C                 OF A DEPENDS ON ALF(J).
C   IINC (IN)  DECLARED LEAD DIMENSION OF INC.  MUST BE AT LEAST L+1.
C     IV (I/O) INTEGER PARAMETER AND SCRATCH VECTOR.
C    LIV (IN)  LENGTH OF IV.  MUST BE AT LEAST
C                 122 + 2*M + 4*P + 2*L + MAX(L+1,6*P), WHERE  M  IS
C                 THE NUMBER OF ONES IN INC.
C     LV (IN)  LENGTH OF V.  MUST BE AT LEAST
C                 105 + 2*N*(L+3) + JLEN + L*(L+3)/2 + P*(2*P + 18),
C                 WHERE  JLEN = (L+P)*(N+L+P+1),  UNLESS NEITHER A
C                 COVARIANCE MATRIX NOR REGRESSION DIAGNOSTICS ARE
C                 REQUESTED, IN WHICH CASE  JLEN = N*P.  IF THE LAST
C                 ROW OF INC CONTAINS ONLY ZEROS, THEN LV CAN BE 4*N
C                 LESS THAN JUST DESCRIBED.
C      V (I/O) FLOATING-POINT PARAMETER AND SCRATCH VECTOR.
C              IF A COVARIANCE ESTIMATE IS REQUESTED, IT IS FOR
C              (ALF,C) -- NONLINEAR PARAMETERS ORDERED FIRST,
C              FOLLOWED BY LINEAR PARAMETERS.
C UIPARM (I/O) INTEGER VECTOR PASSED WITHOUT CHANGE TO CALCA.
C URPARM (I/O) FLOATING-POINT VECTOR PASSED WITHOUT CHANGE TO CALCA.
C UFPARM (I/O) SUBROUTINE PASSED (WITHOUT HAVING BEEN CALLED) TO CALCA.
C
C
C--------------------------  DECLARATIONS  ----------------------------
C
C
C  ***  EXTERNAL SUBROUTINES  ***
C
      EXTERNAL DIVSET, DSM,  DRNSG,DV2AXY,DV7CPY, DV7SCL
C
C DIVSET.... PROVIDES DEFAULT IV AND V VALUES.
C DSM...... DETERMINES EFFICIENT ORDER FOR FINITE DIFFERENCES.
C  DRNSG... CARRIES OUT NL2SOL ALGORITHM.
C DV2AXY.... ADDS A MULTIPLE OF ONE VECTOR TO ANOTHER.
C DV7CPY.... COPIES ONE VECTOR TO ANOTHER.
C DV7SCL... SCALES AND COPIES ONE VECTOR TO ANOTHER.
C
C  ***  LOCAL VARIABLES  ***
C
      LOGICAL PARTJ
      INTEGER A0, A1, AJ, ALP1, BWA1, D0, DA0, DA1, DAJ, GPTR1, GRP1,
     1        GRP2, I, I1, IN0, IN1, IN2, INI, INLEN, IPNTR1, IV1, IWA1,
     2        IWALEN, J1, JN1, JPNTR1, K, L1, LP1, M, M0, NF, NG, NGRP0,
     3        NGRP1, NGRP2, RSAVE0, RSAVE1, RSVLEN, X0I, XSAVE0, XSAVE1
      DOUBLE PRECISION DELTA, DI, H, XI
      DOUBLE PRECISION NEGONE, ONE, ZERO
C
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER AMAT, COVREQ, D, DAMAT, DLTFDJ, GPTR, GRP, IN, IVNEED,
     1        L1SAV, MAXGRP, MODE, MSAVE, NEXTIV, NEXTV, NFCALL, NFGCAL,
     2        PERM, RESTOR, TOOBIG, VNEED, XSAVE
C
C  ***  IV SUBSCRIPT VALUES  ***
C
C/6
C     DATA AMAT/113/, COVREQ/15/, D/27/, DAMAT/114/, DLTFDJ/43/,
C    1     GPTR/117/, GRP/118/, IN/112/, IVNEED/3/, L1SAV/111/,
C    2     MAXGRP/116/, MODE/35/, MSAVE/115/, NEXTIV/46/, NEXTV/47/,
C    3     NFCALL/6/, NFGCAL/7/, PERM/58/, RESTOR/9/, TOOBIG/2/,
C    4     VNEED/4/, XSAVE/119/
C/7
      PARAMETER (AMAT=113, COVREQ=15, D=27, DAMAT=114, DLTFDJ=43,
     1           GPTR=117, GRP=118, IN=112, IVNEED=3, L1SAV=111,
     2           MAXGRP=116, MODE=35, MSAVE=115, NEXTIV=46, NEXTV=47,
     3           NFCALL=6, NFGCAL=7, PERM=58, RESTOR=9, TOOBIG=2,
     4           VNEED=4, XSAVE=119)
C/
      DATA NEGONE/-1.D+0/, ONE/1.D+0/, ZERO/0.D+0/
C
C++++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++
C
      LP1 = L + 1
      IF (IV(1) .EQ. 0) CALL DIVSET(1, IV, LIV, LV, V)
      IF (P .LE. 0 .OR. L .LT. 0 .OR. IINC .LE. L) GO TO 80
      IV1 = IV(1)
      IF (IV1 .EQ. 14) GO TO 120
      IF (IV1 .GT. 2 .AND. IV1 .LT. 12) GO TO 120
      IF (IV1 .EQ. 12) IV(1) = 13
      IF (IV(1) .NE. 13) GO TO 50
C
C  ***  FRESH START ***
C
      IF (IV(PERM) .LE. XSAVE) IV(PERM) = XSAVE + 1
C
C  ***  CHECK INC, COUNT ITS NONZEROS
C
      L1 = 0
      M = 0
      DO 40 I = 1, P
         M0 = M
         IF (L .EQ. 0) GO TO 20
         DO 10 K = 1, L
            IF (INC(K,I) .LT. 0 .OR. INC(K,I) .GT. 1) GO TO 80
            IF (INC(K,I) .EQ. 1) M = M + 1
 10         CONTINUE
 20      IF (INC(LP1,I) .NE. 1) GO TO 30
            M = M + 1
            L1 = 1
 30      IF (M .EQ. M0 .OR. INC(LP1,I) .LT. 0
     1                 .OR. INC(LP1,I) .GT. 1) GO TO 80
 40      CONTINUE
C
C     *** NOW L1 = 1 MEANS A HAS COLUMN L+1 ***
C
C     *** COMPUTE STORAGE REQUIREMENTS ***
C
      IWALEN = MAX0(LP1, 6*P)
      INLEN = 2 * M
      IV(IVNEED) = IV(IVNEED) + INLEN + 3*P + L + IWALEN + 3
      RSVLEN = 2 * L1 * N
      L1 = L + L1
      IV(VNEED) = IV(VNEED) + 2*N*L1 + RSVLEN + P
C
 50   CALL  DRNSG(V, ALF, C, V, IV, IV, L, 1, N, LIV, LV, N, M, P, V, Y)
      IF (IV(1) .NE. 14) GO TO 999
C
C  ***  STORAGE ALLOCATION  ***
C
      IV(IN) = IV(NEXTIV)
      IV(AMAT) = IV(NEXTV)
      IV(DAMAT) = IV(AMAT) + N*L1
      IV(XSAVE) = IV(DAMAT) + N*L1
      IV(NEXTV) = IV(XSAVE) + P + RSVLEN
      IV(L1SAV) = L1
      IV(MSAVE) = M
C
C  ***  DETERMINE HOW MANY GROUPS FOR FINITE DIFFERENCES
C  ***  (SET UP TO CALL DSM)
C
      IN1 = IV(IN)
      JN1 = IN1 + M
      DO 70 K = 1, P
         DO 60 I = 1, LP1
            IF (INC(I,K) .EQ. 0) GO TO 60
               IV(IN1) = I
               IN1 = IN1 + 1
               IV(JN1) = K
               JN1 = JN1 + 1
 60         CONTINUE
 70      CONTINUE
      IN1 = IV(IN)
      JN1 = IN1 + M
      IWA1 = IN1 + INLEN
      NGRP1 = IWA1 + IWALEN
      BWA1 = NGRP1 + P
      IPNTR1 = BWA1 + P
      JPNTR1 = IPNTR1 + L + 2
      CALL DSM(LP1, P, M, IV(IN1), IV(JN1), IV(NGRP1), NG, K, I,
     1         IV(IPNTR1), IV(JPNTR1), IV(IWA1), IWALEN, IV(BWA1))
      IF (I .EQ. 1) GO TO 90
         IV(1) = 69
         GO TO 50
 80   IV(1) = 66
      GO TO 50
C
C  ***  SET UP GRP AND GPTR ARRAYS FOR COMPUTING FINITE DIFFERENCES
C
C  ***  THERE ARE NG GROUPS.  GROUP I CONTAINS ALF(GRP(J)) FOR
C  ***  GPTR(I) .LE. J .LE. GPTR(I+1)-1.
C
 90   IV(MAXGRP) = NG
      IV(GPTR) = IN1 + 2*L1
      GPTR1 = IV(GPTR)
      IV(GRP) = GPTR1 + NG + 1
      IV(NEXTIV) = IV(GRP) + P
      GRP1 = IV(GRP)
      NGRP0 = NGRP1 - 1
      NGRP2 = NGRP0 + P
      DO 110 I = 1, NG
         IV(GPTR1) = GRP1
         GPTR1 = GPTR1 + 1
         DO 100 I1 = NGRP1, NGRP2
            IF (IV(I1) .NE. I) GO TO 100
            IV(GRP1) = I1 - NGRP0
            GRP1 = GRP1 + 1
 100        CONTINUE
 110     CONTINUE
      IV(GPTR1) = GRP1
      IF (IV1 .EQ. 13) GO TO 999
C
C  ***  INITIALIZE POINTERS  ***
C
 120  A1 = IV(AMAT)
      A0 = A1 - N
      DA1 = IV(DAMAT)
      DA0 = DA1 - N
      IN1 = IV(IN)
      IN0 = IN1 - 2
      L1 = IV(L1SAV)
      IN2 = IN1 + 2*L1 - 1
      D0 = IV(D) - 1
      NG = IV(MAXGRP)
      XSAVE1 = IV(XSAVE)
      XSAVE0 = XSAVE1 - 1
      RSAVE1 = XSAVE1 + P
      RSAVE0 = RSAVE1 + N
      ALP1 = A1 + L*N
      DELTA = V(DLTFDJ)
      IV(COVREQ) = -IABS(IV(COVREQ))
C
 130  CALL  DRNSG(V(A1), ALF, C, V(DA1), IV(IN1), IV, L, L1, N, LIV, LV,
     1            N, L1, P, V, Y)
      IF (IV(1)-2) 140, 150, 999
C
C  ***  NEW FUNCTION VALUE (R VALUE) NEEDED  ***
C
 140  NF = IV(NFCALL)
      CALL CALCA(N, P, L, ALF, NF, V(A1), UIPARM, URPARM, UFPARM)
      IF (NF .LE. 0) IV(TOOBIG) = 1
      IF (L1 .LE. L) GO TO 130
      IF (IV(RESTOR) .EQ. 2) CALL DV7CPY(N, V(RSAVE0), V(RSAVE1))
      CALL DV7CPY(N, V(RSAVE1), V(ALP1))
      GO TO 130
C
C  ***  COMPUTE DR = GRADIENT OF R COMPONENTS  ***
C
 150  IF (L1 .GT. L .AND. IV(NFGCAL) .EQ. IV(NFCALL))
     1      CALL DV7CPY(N, V(RSAVE0), V(RSAVE1))
      GPTR1 = IV(GPTR)
      DO 230 K = 1, NG
         CALL DV7CPY(P, V(XSAVE1), ALF)
         GRP1 = IV(GPTR1)
         GRP2 = IV(GPTR1+1) - 1
         GPTR1 = GPTR1 + 1
         DO 160 I1 = GRP1, GRP2
            I = IV(I1)
            XI = ALF(I)
            J1 = D0 + I
            DI = V(J1)
            IF (DI .LE. ZERO) DI = ONE
            H = DELTA * DMAX1(DABS(XI), ONE/DI)
            IF (XI .LT. ZERO) H = -H
            X0I = XSAVE0 + I
            V(X0I) = XI + H
 160        CONTINUE
         CALL CALCA(N, P, L, V(XSAVE1), IV(NFGCAL), V(DA1),
     1              UIPARM, URPARM, UFPARM)
         IF (IV(NFGCAL) .GT. 0) GO TO 170
            IV(TOOBIG) = 1
            GO TO 130
 170     JN1 = IN1
         DO 180 I = IN1, IN2
 180        IV(I) = 0
         PARTJ = IV(MODE) .LE. P
         DO 220 I1 = GRP1, GRP2
            I = IV(I1)
            DO 210 J1 = 1, L1
               IF (INC(J1,I) .EQ. 0) GO TO 210
               INI = IN0 + 2*J1
               IV(INI) = I
               IV(INI+1) = J1
               X0I = XSAVE0 + I
               H = ONE / (V(X0I) - ALF(I))
               DAJ = DA0 + J1*N
               IF (PARTJ) GO TO 190
C                 *** FULL FINITE DIFFERENCE FOR COV. AND REG. DIAG. ***
                  AJ = A0 + J1*N
                  CALL DV2AXY(N, V(DAJ), NEGONE, V(AJ), V(DAJ))
                  GO TO 200
 190           IF (J1 .GT. L)
     1            CALL DV2AXY(N, V(DAJ), NEGONE, V(RSAVE0), V(DAJ))
 200           CALL DV7SCL(N, V(DAJ), H, V(DAJ))
 210           CONTINUE
 220        CONTINUE
         IF (K .GE. NG) GO TO 240
         IV(1) = -2
         CALL  DRNSG(V(A1), ALF, C, V(DA1), IV(IN1), IV, L, L1, N, LIV,
     1               LV, N, L1, P, V, Y)
         IF (-2 .NE. IV(1)) GO TO 999
 230     CONTINUE
 240  IV(1) = 2
      GO TO 130
C
 999  RETURN
C
C  ***  LAST CARD OF   DNSF FOLLOWS  ***
      END
