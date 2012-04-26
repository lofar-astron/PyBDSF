      SUBROUTINE MOVEBR(N,A,B)
C
C     MOVEBR MOVES N REAL ITEMS FROM A TO B
C     USING A BACKWARDS DO LOOP
C
      REAL A(1),B(1)
C
      I = N
C
 10   IF(I .LE. 0) RETURN
        B(I) = A(I)
        I = I - 1
        GO TO 10
C
      END
