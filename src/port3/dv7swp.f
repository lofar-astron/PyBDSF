      SUBROUTINE DV7SWP(N, X, Y)
C
C  ***  INTERCHANGE N-VECTORS X AND Y.  ***
C
      INTEGER N
      DOUBLE PRECISION X(N), Y(N)
C
      INTEGER I
      DOUBLE PRECISION T
C
      DO 10 I = 1, N
         T = X(I)
         X(I) = Y(I)
         Y(I) = T
 10      CONTINUE
 999  RETURN
C  ***  LAST CARD OF DV7SWP FOLLOWS  ***
      END
