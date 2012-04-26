      SUBROUTINE  S7LUP(A, COSMIN, P, SIZE, STEP, U, W, WCHMTD, WSCALE,
     1                  Y)
C
C  ***  UPDATE SYMMETRIC  A  SO THAT  A * STEP = Y  ***
C  ***  (LOWER TRIANGLE OF  A  STORED ROWWISE       ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER P
      REAL A(1), COSMIN, SIZE, STEP(P), U(P), W(P),
     1                 WCHMTD(P), WSCALE, Y(P)
C     DIMENSION A(P*(P+1)/2)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, J, K
      REAL DENMIN, SDOTWM, T, UI, WI
C
C     ***  CONSTANTS  ***
      REAL HALF, ONE, ZERO
C
C  ***  EXTERNAL FUNCTIONS AND SUBROUTINES  ***
C
      REAL  D7TPR,  V2NRM
      EXTERNAL  D7TPR,  S7LVM,  V2NRM
C
C/6
C     DATA HALF/0.5E+0/, ONE/1.E+0/, ZERO/0.E+0/
C/7
      PARAMETER (HALF=0.5E+0, ONE=1.E+0, ZERO=0.E+0)
C/
C
C-----------------------------------------------------------------------
C
      SDOTWM =  D7TPR(P, STEP, WCHMTD)
      DENMIN = COSMIN *  V2NRM(P,STEP) *  V2NRM(P,WCHMTD)
      WSCALE = ONE
      IF (DENMIN .NE. ZERO) WSCALE = AMIN1(ONE,  ABS(SDOTWM/DENMIN))
      T = ZERO
      IF (SDOTWM .NE. ZERO) T = WSCALE / SDOTWM
      DO 10 I = 1, P
 10      W(I) = T * WCHMTD(I)
      CALL  S7LVM(P, U, A, STEP)
      T = HALF * (SIZE *  D7TPR(P, STEP, U)  -   D7TPR(P, STEP, Y))
      DO 20 I = 1, P
 20      U(I) = T*W(I) + Y(I) - SIZE*U(I)
C
C  ***  SET  A = A + U*(W**T) + W*(U**T)  ***
C
      K = 1
      DO 40 I = 1, P
         UI = U(I)
         WI = W(I)
         DO 30 J = 1, I
              A(K) = SIZE*A(K) + UI*W(J) + WI*U(J)
              K = K + 1
 30           CONTINUE
 40      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF  S7LUP FOLLOWS  ***
      END
