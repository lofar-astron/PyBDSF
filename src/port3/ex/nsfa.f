C$TEST NSFA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NSFA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM NSF1
C
C***********************************************************************
C  EXAMPLE PROGRAM FOR NSF1 TO FIT
C  N DATA POINTS (T,Y) TO CURVE
C  C(1)*EXP(T*X)  +  C(2)
C
       INTEGER IWRITE
       REAL C(2), T(8), Y(8), TT(8), YY(8)
       DOUBLE PRECISION S
       EXTERNAL GETAY
       COMMON /DATBLK/TT,YY
       DATA T(1) /12.0/, T(2) /20.0/ ,T(3) /28.0/, T(4) /48.0/,
     1 T(5)/120.0/, T(6) /240.0/, T(7) /900.0/, T(8) /2400.0/
       DATA Y(1) /0.2342/, Y(2) /0.2244/ , Y(3) /0.2204/,
     1 Y(4) /0.2149/, Y(5) /0.2063/, Y(6) /0.1983/,
     2 Y(7) /0.1842/, Y(8)/0.1761/
C
C  SET UP OUTPUT UNIT
C
       IWRITE = I1MACH(2)
C
C  MOVE T AND Y VECTORS TO COMMON
C
       DO 2 I=1,8
         TT(I) = T(I)
         YY(I) = Y(I)
  2    CONTINUE
C
        N = 8
        L =  2
       X1 = -10.0
       X2 =   0.001
C
C  DO THE FIT
C
       CALL NSF1(N, L, X, X1, X2, 1.E-6, C)
       WRITE(IWRITE, 4) X, C(1), C(2)
  4      FORMAT(5H X = , E20.10/8H C(1) = ,E20.10/8H C(2) = , E20.10)
C
       WRITE(IWRITE, 5)
  5      FORMAT(//,19X,1HT,14X,6HREAL Y,14X,5HEST.Y,15X,5HERROR,/)
       DO 100 I=1,N
         YEST = C(1)*EXP(T(I)*X)+C(2)
         YERR = ABS(Y(I)-YEST)
         WRITE(IWRITE, 6) T(I), Y(I), YEST, YERR
 100     S = S + YERR*YERR
   6     FORMAT (4E20.10)
       WRITE(IWRITE, 7) S
   7     FORMAT(//,24HSUM OF ERRORS SQUARED = ,D20.10)
       STOP
       END
       SUBROUTINE GETAY(N,L,X,A,Y)
       INTEGER N,L
       REAL A(N,L),X,Y(N)
       REAL T(8),YY(8)
       COMMON /DATBLK/T,YY
       DO 100 I=1,N
          A(I,1)=EXP(X*T(I))
          A(I,2)=1.0
          Y(I)=YY(I)
 100   CONTINUE
       RETURN
       END
