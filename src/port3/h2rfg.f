      REAL FUNCTION  H2RFG(A, B, X, Y, Z)
C
C  ***  DETERMINE X, Y, Z SO  I + (1,Z)**T * (X,Y)  IS A 2X2
C  ***  HOUSEHOLDER REFLECTION SENDING (A,B)**T INTO (C,0)**T,
C  ***  WHERE  C = -SIGN(A)*SQRT(A**2 + B**2)  IS THE VALUE  H2RFG
C  ***  RETURNS.
C
      REAL A, B, X, Y, Z
C
      REAL A1, B1, C, T
C/+
      REAL  SQRT
C/
      REAL ZERO
      DATA ZERO/0.E+0/
C
C  ***  BODY  ***
C
      IF (B .NE. ZERO) GO TO 10
         X = ZERO
         Y = ZERO
         Z = ZERO
          H2RFG = A
         GO TO 999
 10   T =  ABS(A) +  ABS(B)
      A1 = A / T
      B1 = B / T
      C =  SQRT(A1**2 + B1**2)
      IF (A1 .GT. ZERO) C = -C
      A1 = A1 - C
      Z = B1 / A1
      X = A1 / C
      Y = B1 / C
       H2RFG = T * C
 999  RETURN
C  ***  LAST LINE OF  H2RFG FOLLOWS  ***
      END
