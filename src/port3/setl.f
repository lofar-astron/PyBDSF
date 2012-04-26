      SUBROUTINE SETL(N,V,B)
C
C     SETL SETS THE N LOGICAL ITEMS IN B TO V
C
      LOGICAL B(1),V
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
 10     B(I) = V
C
      RETURN
C
      END
