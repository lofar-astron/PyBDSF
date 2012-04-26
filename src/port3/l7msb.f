      SUBROUTINE  L7MSB(B, D, G, IERR, IPIV, IPIV1, IPIV2, KA, LMAT,
     1                  LV, P, P0, PC, QTR, RMAT, STEP, TD, TG, V,
     2                  W, WLM, X, X0)
C
C  ***  COMPUTE HEURISTIC BOUNDED NEWTON STEP  ***
C
      INTEGER IERR, KA, LV, P, P0, PC
      INTEGER IPIV(P), IPIV1(P), IPIV2(P)
      REAL B(2,P), D(P), G(P), LMAT(1), QTR(P), RMAT(1),
     1                 STEP(P,3), TD(P), TG(P), V(LV), W(P), WLM(1),
     2                 X0(P), X(P)
C     DIMENSION LMAT(P*(P+1)/2), RMAT(P*(P+1)/2), WLM(P*(P+5)/2 + 4)
C
      REAL  D7TPR
      EXTERNAL  D7MLP,  D7TPR,  L7MST,  L7TVM,  Q7RSH,  S7BQN,
     1         V2AXY, V7CPY,  V7IPR,  V7SCP,  V7VMP
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K, K0, KB, KINIT, L, NS, P1, P10, P11
      REAL DS0, NRED, PRED, RAD
      REAL ONE, ZERO
C
C  ***  V SUBSCRIPTS  ***
C
      INTEGER DST0, DSTNRM, GTSTEP, NREDUC, PREDUC, RADIUS
C
C/6
C     DATA DST0/3/, DSTNRM/2/, GTSTEP/4/, NREDUC/6/, PREDUC/7/,
C    1     RADIUS/8/
C/7
      PARAMETER (DST0=3, DSTNRM=2, GTSTEP=4, NREDUC=6, PREDUC=7,
     1           RADIUS=8)
C/
      DATA ONE/1.E+0/, ZERO/0.E+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      P1 = PC
      IF (KA .LT. 0) GO TO 10
         NRED = V(NREDUC)
         DS0 = V(DST0)
         GO TO 20
 10   P0 = 0
      KA = -1
C
 20   KINIT = -1
      IF (P0 .EQ. P1) KINIT = KA
      CALL V7CPY(P, X, X0)
      CALL V7CPY(P, TD, D)
C     *** _USE_ STEP(1,3) AS TEMP. COPY OF QTR ***
      CALL V7CPY(P, STEP(1,3), QTR)
      CALL  V7IPR(P, IPIV, TD)
      PRED = ZERO
      RAD = V(RADIUS)
      KB = -1
      V(DSTNRM) = ZERO
      IF (P1 .GT. 0) GO TO 30
         NRED = ZERO
         DS0 = ZERO
         CALL  V7SCP(P, STEP, ZERO)
         GO TO 90
C
 30   CALL  V7VMP(P, TG, G, D, -1)
      CALL  V7IPR(P, IPIV, TG)
      P10 = P1
 40   K = KINIT
      KINIT = -1
      V(RADIUS) = RAD - V(DSTNRM)
      CALL  V7VMP(P1, TG, TG, TD, 1)
      DO 50 I = 1, P1
 50      IPIV1(I) = I
      K0 = MAX0(0, K)
      CALL  L7MST(TD, TG, IERR, IPIV1, K, P1, STEP(1,3), RMAT, STEP,
     1            V, WLM)
      CALL  V7VMP(P1, TG, TG, TD, -1)
      P0 = P1
      IF (KA .GE. 0) GO TO 60
         NRED = V(NREDUC)
         DS0 = V(DST0)
C
 60   KA = K
      V(RADIUS) = RAD
      L = P1 + 5
      IF (K .LE. K0) CALL  D7MLP(P1, LMAT, TD, RMAT, -1)
      IF (K .GT. K0) CALL  D7MLP(P1, LMAT, TD, WLM(L), -1)
      CALL  S7BQN(B, D, STEP(1,2), IPIV, IPIV1, IPIV2, KB, LMAT,
     1            LV, NS, P, P1, STEP, TD, TG, V, W, X, X0)
      PRED = PRED + V(PREDUC)
      IF (NS .EQ. 0) GO TO 80
      P0 = 0
C
C  ***  UPDATE RMAT AND QTR  ***
C
      P11 = P1 + 1
      L = P10 + P11
      DO 70 K = P11, P10
         J = L - K
         I = IPIV2(J)
         IF (I .LT. J) CALL  Q7RSH(I, J, .TRUE., QTR, RMAT, W)
 70      CONTINUE
C
 80   IF (KB .GT. 0) GO TO 90
C
C  ***  UPDATE LOCAL COPY OF QTR  ***
C
      CALL  V7VMP(P10, W, STEP(1,2), TD, -1)
      CALL  L7TVM(P10, W, LMAT, W)
      CALL V2AXY(P10, STEP(1,3), ONE, W, QTR)
      GO TO 40
C
 90   V(DST0) = DS0
      V(NREDUC) = NRED
      V(PREDUC) = PRED
      V(GTSTEP) =  D7TPR(P, G, STEP)
C
 999  RETURN
C  ***  LAST LINE OF  L7MSB FOLLOWS  ***
      END
