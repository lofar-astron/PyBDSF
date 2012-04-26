      SUBROUTINE V7SWP(N, X, Y)
C
C  ***  INTERCHANGE N-VECTORS X AND Y.  ***
C
      INTEGER N
      REAL X(N), Y(N)
C
      INTEGER I
      REAL T
C
      DO 10 I = 1, N
         T = X(I)
         X(I) = Y(I)
         Y(I) = T
 10      CONTINUE
 999  RETURN
C  ***  LAST CARD OF V7SWP FOLLOWS  ***
      END
