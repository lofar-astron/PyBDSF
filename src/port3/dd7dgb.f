        SUBROUTINE DD7DGB(B, D, DIG, DST, G, IPIV, KA, L, LV, P, PC,
     1                    NWTST, STEP, TD, TG, V, W, X0)
C
C  ***  COMPUTE DOUBLE-DOGLEG STEP, SUBJECT TO SIMPLE BOUNDS ON X  ***
C
      INTEGER LV, KA, P, PC
      INTEGER IPIV(P)
      DOUBLE PRECISION B(2,P), D(P), DIG(P), DST(P), G(P), L(1),
     1                 NWTST(P), STEP(P), TD(P), TG(P), V(LV), W(P),
     2                 X0(P)
C
C     DIMENSION L(P*(P+1)/2)
C
      DOUBLE PRECISION DD7TPR, DR7MDC, DV2NRM
      EXTERNAL DD7DOG, DD7TPR, I7SHFT, DL7ITV, DL7IVM, DL7TVM,DL7VML,
     1         DQ7RSH, DR7MDC, DV2NRM,DV2AXY,DV7CPY, DV7IPR, DV7SCP,
     2         DV7SHF, DV7VMP
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K, P1, P1M1
      DOUBLE PRECISION DNWTST, GHINVG, GNORM, GNORM0, NRED, PRED, RAD,
     1                 T, T1, T2, TI, X0I, XI
      DOUBLE PRECISION HALF, MEPS2, ONE, TWO, ZERO
C
C  ***  V SUBSCRIPTS  ***
C
      INTEGER DGNORM, DST0, DSTNRM, GRDFAC, GTHG, GTSTEP, NREDUC,
     1        NWTFAC, PREDUC, RADIUS, STPPAR
C
C/6
C     DATA DGNORM/1/, DST0/3/, DSTNRM/2/, GRDFAC/45/, GTHG/44/,
C    1     GTSTEP/4/, NREDUC/6/, NWTFAC/46/, PREDUC/7/, RADIUS/8/,
C    2     STPPAR/5/
C/7
      PARAMETER (DGNORM=1, DST0=3, DSTNRM=2, GRDFAC=45, GTHG=44,
     1           GTSTEP=4, NREDUC=6, NWTFAC=46, PREDUC=7, RADIUS=8,
     2           STPPAR=5)
C/
C/6
C     DATA HALF/0.5D+0/, ONE/1.D+0/, TWO/2.D+0/, ZERO/0.D+0/
C/7
      PARAMETER (HALF=0.5D+0, ONE=1.D+0, TWO=2.D+0, ZERO=0.D+0)
      SAVE MEPS2
C/
      DATA MEPS2/0.D+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      IF (MEPS2 .LE. ZERO) MEPS2 = TWO * DR7MDC(3)
      GNORM0 = V(DGNORM)
      V(DSTNRM) = ZERO
      IF (KA .LT. 0) GO TO 10
         DNWTST = V(DST0)
         NRED = V(NREDUC)
 10   PRED = ZERO
      V(STPPAR) = ZERO
      RAD = V(RADIUS)
      IF (PC .GT. 0) GO TO 20
         DNWTST = ZERO
         CALL DV7SCP(P, STEP, ZERO)
         GO TO 140
C
 20   P1 = PC
      CALL DV7CPY(P, TD, D)
      CALL DV7IPR(P, IPIV, TD)
      CALL DV7SCP(PC, DST, ZERO)
      CALL DV7CPY(P, TG, G)
      CALL DV7IPR(P, IPIV, TG)
C
 30   CALL DL7IVM(P1, NWTST, L, TG)
      GHINVG = DD7TPR(P1, NWTST, NWTST)
      V(NREDUC) = HALF * GHINVG
      CALL DL7ITV(P1, NWTST, L, NWTST)
      CALL DV7VMP(P1, STEP, NWTST, TD, 1)
      V(DST0) = DV2NRM(PC, STEP)
      IF (KA .GE. 0) GO TO 40
         KA = 0
         DNWTST = V(DST0)
         NRED = V(NREDUC)
 40   V(RADIUS) = RAD - V(DSTNRM)
      IF (V(RADIUS) .LE. ZERO) GO TO 100
      CALL DV7VMP(P1, DIG, TG, TD, -1)
      GNORM = DV2NRM(P1, DIG)
      IF (GNORM .LE. ZERO) GO TO 100
      V(DGNORM) = GNORM
      CALL DV7VMP(P1, DIG, DIG, TD, -1)
      CALL DL7TVM(P1, W, L, DIG)
      V(GTHG) = DV2NRM(P1, W)
      KA = KA + 1
      CALL DD7DOG(DIG, LV, P1, NWTST, STEP, V)
C
C     ***  FIND T SUCH THAT X - T*STEP IS STILL FEASIBLE.
C
      T = ONE
      K = 0
      DO 70 I = 1, P1
         J = IPIV(I)
         X0I = X0(J) + DST(I)/TD(I)
         XI = X0I + STEP(I)
         IF (XI .LT. B(1,J)) GO TO 50
         IF (XI .LE. B(2,J)) GO TO 70
              TI = (B(2,J) - X0I) / STEP(I)
              J = I
              GO TO 60
 50      TI = (B(1,J) - X0I) / STEP(I)
         J = -I
 60      IF (T .LE. TI) GO TO 70
              K = J
              T = TI
 70      CONTINUE
C
C  ***  UPDATE DST, TG, AND PRED  ***
C
      CALL DV7VMP(P1, STEP, STEP, TD, 1)
      CALL DV2AXY(P1, DST, T, STEP, DST)
      V(DSTNRM) = DV2NRM(PC, DST)
      T1 = T * V(GRDFAC)
      T2 = T * V(NWTFAC)
      PRED = PRED - T1*GNORM * ((T2 + ONE)*GNORM)
     1                 - T2 * (ONE + HALF*T2)*GHINVG
     2                  - HALF * (V(GTHG)*T1)**2
      IF (K .EQ. 0) GO TO 100
      CALL DL7VML(P1, W, L, W)
      T2 = ONE - T2
      DO 80 I = 1, P1
 80      TG(I) = T2*TG(I) - T1*W(I)
C
C     ***  PERMUTE L, ETC. IF NECESSARY  ***
C
      P1M1 = P1 - 1
      J = IABS(K)
      IF (J .EQ. P1) GO TO 90
         CALL DQ7RSH(J, P1, .FALSE., TG, L, W)
         CALL I7SHFT(P1, J, IPIV)
         CALL DV7SHF(P1, J, TG)
         CALL DV7SHF(P1, J, TD)
         CALL DV7SHF(P1, J, DST)
 90   IF (K .LT. 0) IPIV(P1) = -IPIV(P1)
      P1 = P1M1
      IF (P1 .GT. 0) GO TO 30
C
C     ***  UNSCALE STEP, UPDATE X AND DIHDI  ***
C
 100  CALL DV7SCP(P, STEP, ZERO)
      DO 110 I = 1, PC
         J = IABS(IPIV(I))
         STEP(J) = DST(I) / TD(I)
 110     CONTINUE
C
C  ***  FUDGE STEP TO ENSURE THAT IT FORCES APPROPRIATE COMPONENTS
C  ***  TO THEIR BOUNDS  ***
C
      IF (P1 .GE. PC) GO TO 140
      CALL DV2AXY(P, TD, ONE, STEP, X0)
      K = P1 + 1
      DO 130 I = K, PC
         J = IPIV(I)
         T = MEPS2
         IF (J .GT. 0) GO TO 120
            T = -T
            J = -J
            IPIV(I) = J
 120     T = T * DMAX1(DABS(TD(J)), DABS(X0(J)))
         STEP(J) = STEP(J) + T
 130     CONTINUE
C
 140  V(DGNORM) = GNORM0
      V(NREDUC) = NRED
      V(PREDUC) = PRED
      V(RADIUS) = RAD
      V(DST0) = DNWTST
      V(GTSTEP) = DD7TPR(P, STEP, G)
C
 999  RETURN
C  ***  LAST LINE OF DD7DGB FOLLOWS  ***
      END
