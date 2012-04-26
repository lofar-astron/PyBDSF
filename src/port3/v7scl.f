      SUBROUTINE  V7SCL(N, X, A, Y)
C
C  ***  SET X(I) = A*Y(I), I = 1(1)N  ***
C
      INTEGER N
      REAL A, X(N), Y(N)
C
      INTEGER I
C
      DO 10 I = 1, N
 10       X(I) = A * Y(I)
 999    RETURN
C  ***  LAST LINE OF  V7SCL FOLLOWS  ***
      END
