      SUBROUTINE SETC(N,V,B)
C
C     SETC SETS THE N COMPLEX ITEMS IN B TO V
C
C/R
C     REAL B(2,N), V(2), V1, V2
C     V1 = V(1)
C     V2 = V(2)
C/C
      COMPLEX B(1),V
C/
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
C/R
C       B(1,I) = V1
C10     B(2,I) = V2
C/C
 10     B(I) = V
C/
C
      RETURN
C
      END
