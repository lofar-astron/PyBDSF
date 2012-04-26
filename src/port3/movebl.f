      SUBROUTINE MOVEBL(N,A,B)
C
C     MOVEBL MOVES N LOGICAL ITEMS FROM A TO B
C     USING A BACKWARDS DO LOOP
C
      LOGICAL A(1),B(1)
C
      I = N
C
 10   IF(I .LE. 0) RETURN
        B(I) = A(I)
        I = I - 1
        GO TO 10
C
      END
