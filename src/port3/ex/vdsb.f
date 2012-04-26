C$TEST VDSB
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDSB
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM VDSS2
C
C***********************************************************************
C  SIN(X(1)) + X(2)**2 + 3
C
       REAL X(2), A(10,10), F, FPRIME(2), PT1, PT2, H1, H2, SOL
       REAL DSOL1, DSOL2
       INTEGER IWRITE, N1, N2, NA1, I, J
       REAL U(2),L(2),DIST(2)
C
       IWRITE = I1MACH(2)
       N1 = 10
       N2 = 10
       X(1) = 1.7
       X(2) = 1.7
       H1 = 1./ FLOAT(N1-1)
       H2 = 1./ FLOAT(N2-1)
       L(1) = 1.0
       U(1) = 5.0
       L(2) = 1.0
       U(2) = 5.0
       DIST(1) = U(1) - L(1)
       DIST(2) = U(2) - L(2)
C SET UP THE MESH AND DATA VALUES
       DO 100 I=1,N1
              PT1 = L(1) + DIST(1)*FLOAT(I-1)*H1
              DO 100 J=1,N2
                     PT2 = L(2) + DIST(2)*FLOAT(J-1)*H2
                     A(I,J) = SIN(PT1) + PT2**2 + 3.0
 100   CONTINUE
       NA1 = 10
       CALL VDSS2 (X,N1,N2,U,L,A,NA1,F,FPRIME)
C CHECK THE SOLUTION
       SOL = SIN(X(1)) + X(2)**2 + 3.0
       DSOL1 = COS(X(1))
       DSOL2 = 2.0*X(2)
       WRITE (IWRITE,101)
 101   FORMAT (45H                     ACTUAL          COMPUTED//)
       WRITE (IWRITE,102) SOL,F
 102   FORMAT (17H          F(X) = ,2E16.8)
       WRITE (IWRITE,103) DSOL1,FPRIME(1)
 103   FORMAT (17H     PARTIAL X = ,2E16.8)
       WRITE (IWRITE,104) DSOL2,FPRIME(2)
 104   FORMAT (17H     PARTIAL Y = ,2E16.8)
       STOP
       END
