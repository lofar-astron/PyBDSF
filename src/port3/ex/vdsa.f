C$TEST VDSA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE VDSA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM VDSS1
C
C***********************************************************************
C  SIN(X) + X**2 + 3
C
       REAL X, A(10), F, FPRIME, PT, H, SOL
       INTEGER IWRITE, N, I
       REAL U, L, DIST, DSOL
C
       IWRITE = I1MACH(2)
       N = 10
       X = 1.3
       H = 1./ FLOAT(N-1)
       L = 1.0
       U = 5.0
       DIST = U - L
C SET UP THE MESH AND DATA VALUES
       DO 100 I=1,N
              PT = L + DIST*FLOAT(I-1)*H
              A(I) = SIN(PT) + PT**2 + 3.0
 100   CONTINUE
       CALL VDSS1 (X,N,U,L,A,F,FPRIME)
C CHECK THE SOLUTION
       SOL = SIN(X) + X**2 + 3.0
       DSOL = COS(X) + 2.0*X
       WRITE (IWRITE,101)
 101   FORMAT (45H                     ACTUAL          COMPUTED//)
       WRITE (IWRITE,102) SOL,F
 102   FORMAT (17H          F(X) = ,2E16.8)
       WRITE (IWRITE,103) DSOL,FPRIME
 103   FORMAT (17H    DERIVATIVE = ,2E16.8)
       STOP
       END
