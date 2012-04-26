      SUBROUTINE SETR(N,V,B)
C
C     SETR SETS THE N REAL ITEMS IN B TO V
C
      REAL B(1),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
