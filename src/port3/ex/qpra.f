C$TEST QPRA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE QPRA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM IQP
C
C***********************************************************************
C TEST PROGRAM FOR IQP
        REAL X(10), Q(10,10), A(10,10), BL(10), BU(10)
        REAL C(10), B(10)
        REAL SUM(10), FUNCT
        INTEGER N, I, J, IPRINT, MAXITR, IQ, M, IA
        INTEGER IEQ
        DOUBLE PRECISION DSTAK(2000)
        COMMON /CSTAK/DSTAK
C
        IWRITE = I1MACH(2)
        CALL ISTKIN(2000,4)
        N = 4
        M = 3
C SET UP INITIAL GUESS AND QUADRATIC FUNCTION
        DO 1 I=1,N
            X(I) = I + 1.
            C(I) = 8. - I
            DO 2 J=1,N
               Q(I,J) = FLOAT(IABS(I-J))
 2          CONTINUE
            Q(I,I) = 1.69
C SET UP GENERAL CONSTRAINTS
            DO 16 J=1,M
               A(J,I) = 0.
 16         CONTINUE
 1      CONTINUE
        DO 3 I=1,M
            B(I) = -1. - (I - 1.) * .05
            A(I,I) = -1.
            A(I,I+1) = 1.
 3      CONTINUE
        IQ = 10
        IA = 10
        IEQ = 1
C SET UP SIMPLE CONSTRAINTS
        DO 4 I=1,N
            BL(I) = -I - (I - 1.) * .1
            BU(I) = I
 4      CONTINUE
C GET MACHINE INFINITY FROM PORT
        BU(1) = R1MACH(2)
        IPRINT = 1
        MAXITR = 3*N
C CALL THE QUADRATIC PROGRAMMING PACKAGE
        CALL IQP(N, X, Q, IQ, C, M, A, IA, B, BL, BU, IPRINT,
     1          MAXITR, IEQ)
C COMPUTE FINAL FUNCTION VALUE
        DO 6 J=1,N
            SUM(J) = X(J) * Q(J,J)
 6      CONTINUE
        DO 7 I=2,N
            DO 9 J=1,I-1
                SUM(I) = SUM(I) + X(J)*Q(J,I)
                SUM(J) = SUM(J) + X(I)*Q(J,I)
 9          CONTINUE
 7      CONTINUE
        FUNCT = 0.
        DO 10 I=1,N
            FUNCT = SUM(I) * X(I)/2. + FUNCT + C(I) * X(I)
 10     CONTINUE
        WRITE (IWRITE,1000)
 1000    FORMAT (16H FINAL SOLUTION:)
        DO 11 I=1,N
           WRITE (IWRITE, 1001) I, X(I)
 1001       FORMAT (I5,D14.4)
 11     CONTINUE
        WRITE (IWRITE,1002)
 1002    FORMAT (22H FINAL FUNCTION VALUE:)
        WRITE (IWRITE,1003) FUNCT
 1003    FORMAT (D14.4)
        STOP
        END
