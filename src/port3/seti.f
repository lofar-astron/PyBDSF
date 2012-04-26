      SUBROUTINE SETI(N,V,B)
C
C     SETI SETS THE N INTEGER ITEMS IN B TO V
C
      INTEGER B(1),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
