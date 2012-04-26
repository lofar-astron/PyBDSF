      SUBROUTINE MOVEFC(N,A,B)
C
C     MOVEFC MOVES N COMPLEX ITEMS FROM A TO B
C     USING A FORWARDS DO LOOP
C
C/R
C     REAL A(2,N), B(2,N)
C/C
      COMPLEX A(1),B(1)
C/
C
      IF(N .LE. 0) RETURN
C
      DO 10 I = 1, N
C/R
C       B(1,I) = A(1,I)
C10     B(2,I) = A(2,I)
C/C
 10     B(I) = A(I)
C/
C
      RETURN
C
      END
