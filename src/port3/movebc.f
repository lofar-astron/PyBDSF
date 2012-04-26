      SUBROUTINE MOVEBC(N,A,B)
C
C     MOVEBC MOVES N COMPLEX ITEMS FROM A TO B
C     USING A BACKWARDS DO LOOP
C
C/R
C     REAL A(2,N), B(2,N)
C/C
      COMPLEX A(1),B(1)
C/
C
      I = N
C
 10   IF(I .LE. 0) RETURN
C/R
C       B(2,I) = A(2,I)
C       B(1,I) = A(1,I)
C/C
        B(I) = A(I)
C/
        I = I - 1
        GO TO 10
C
      END
