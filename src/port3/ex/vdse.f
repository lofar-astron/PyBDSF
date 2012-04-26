C$TEST VDSE
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDSE
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM VDSS3
C
C***********************************************************************
C  X(1)*X(2) + EXP(X(2)) + X(3)**2
C
       REAL X(3), A(10,10,10), F, FPRIME(3), PT1, PT2, PT3
       REAL H1, H2, H3, SOL
       REAL DSOL1, DSOL2, DSOL3
       INTEGER IWRITE, N1, N2, N3, NA1, NA2
       INTEGER I,J,K
       REAL U(3),L(3),DIST(3)
C
       IWRITE = I1MACH(2)
       N1 = 10
       N2 = 10
       N3 = 10
       X(1) = 2.3
       X(2) = 2.3
       X(3) = 2.3
       H1 = 1./ FLOAT(N1-1)
       H2 = 1./ FLOAT(N2-1)
       H3 = 1./ FLOAT(N3-1)
       L(1) = 1.0
       U(1) = 3.0
       L(2) = 2.0
       U(2) = 4.0
       L(3) = 1.5
       U(3) = 3.5
       DIST(1) = U(1) - L(1)
       DIST(2) = U(2) - L(2)
       DIST(3) = U(3) - L(3)
C SET UP THE MESH AND DATA VALUES
       DO 100 I=1,N1
              PT1 = L(1) + DIST(1)*FLOAT(I-1)*H1
              DO 100 J=1,N2
                     PT2 = L(2) + DIST(2)*FLOAT(J-1)*H2
                     DO 100 K=1,N3
                            PT3 = L(3) + DIST(3)*FLOAT(K-1)*H3
                            A(I,J,K) = PT1*PT2 + EXP(PT2) + PT3**2
 100   CONTINUE
       NA1 = 10
       NA2 = 10
       CALL VDSS3 (X,N1,N2,N3,U,L,A,NA1,NA2,F,FPRIME)
C CHECK THE SOLUTION
       SOL = X(1)*X(2) + EXP(X(2)) + X(3)**2
       DSOL1 = X(2)
       DSOL2 = X(1) + EXP(X(2))
       DSOL3 = 2.0*X(3)
       WRITE (IWRITE,101)
 101   FORMAT (45H                     ACTUAL          COMPUTED//)
       WRITE (IWRITE,102) SOL,F
 102   FORMAT (17H          F(X) = ,2E16.8)
       WRITE (IWRITE,103) DSOL1,FPRIME(1)
 103   FORMAT (17H     PARTIAL X = ,2E16.8)
       WRITE (IWRITE,104) DSOL2,FPRIME(2)
 104   FORMAT (17H     PARTIAL Y = ,2E16.8)
       WRITE (IWRITE,105) DSOL3,FPRIME(3)
 105   FORMAT (17H     PARTIAL Z = ,2E16.8)
       STOP
       END
