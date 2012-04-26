      SUBROUTINE  L7IVM(N, X, L, Y)
C
C  ***  SOLVE  L*X = Y, WHERE  L  IS AN  N X N  LOWER TRIANGULAR
C  ***  MATRIX STORED COMPACTLY BY ROWS.  X AND Y MAY OCCUPY THE SAME
C  ***  STORAGE.  ***
C
      INTEGER N
      REAL X(N), L(1), Y(N)
      REAL  D7TPR
      EXTERNAL  D7TPR
      INTEGER I, J, K
      REAL T, ZERO
C/6
C     DATA ZERO/0.E+0/
C/7
      PARAMETER (ZERO=0.E+0)
C/
C
      DO 10 K = 1, N
         IF (Y(K) .NE. ZERO) GO TO 20
         X(K) = ZERO
 10      CONTINUE
      GO TO 999
 20   J = K*(K+1)/2
      X(K) = Y(K) / L(J)
      IF (K .GE. N) GO TO 999
      K = K + 1
      DO 30 I = K, N
         T =  D7TPR(I-1, L(J+1), X)
         J = J + I
         X(I) = (Y(I) - T)/L(J)
 30      CONTINUE
 999  RETURN
C  ***  LAST CARD OF  L7IVM FOLLOWS  ***
      END
