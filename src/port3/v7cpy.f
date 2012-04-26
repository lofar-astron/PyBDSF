      SUBROUTINE V7CPY(P, Y, X)
C
C  ***  SET Y = X, WHERE X AND Y ARE P-VECTORS  ***
C
      INTEGER P
      REAL X(P), Y(P)
C
      INTEGER I
C
      DO 10 I = 1, P
 10      Y(I) = X(I)
      RETURN
      END
