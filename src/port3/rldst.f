      REAL FUNCTION  RLDST(P, D, X, X0)
C
C  ***  COMPUTE AND RETURN RELATIVE DIFFERENCE BETWEEN X AND X0  ***
C  ***  NL2SOL VERSION 2.2  ***
C
      INTEGER P
      REAL D(P), X(P), X0(P)
C
      INTEGER I
      REAL EMAX, T, XMAX, ZERO
C/6
C     DATA ZERO/0.E+0/
C/7
      PARAMETER (ZERO=0.E+0)
C/
C
C  ***  BODY  ***
C
      EMAX = ZERO
      XMAX = ZERO
      DO 10 I = 1, P
         T =  ABS(D(I) * (X(I) - X0(I)))
         IF (EMAX .LT. T) EMAX = T
         T = D(I) * ( ABS(X(I)) +  ABS(X0(I)))
         IF (XMAX .LT. T) XMAX = T
 10      CONTINUE
       RLDST = ZERO
      IF (XMAX .GT. ZERO)  RLDST = EMAX / XMAX
 999  RETURN
C  ***  LAST CARD OF  RLDST FOLLOWS  ***
      END
