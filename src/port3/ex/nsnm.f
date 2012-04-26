C$TEST NSNM
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NSNM
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SMNSX
C
C***********************************************************************
C *** SMNSX EXAMPLE PROGRAM ***
C
C *** MINIMIZE F(X) = 0.1*S(X)**4 + SUM(I = 1(1)3) (I * (X(I) - 10)**2),
C *** WHERE S(X) = SUM(I = 1(1)3) X(I),
C *** STARTING FROM     X = (2, 30, 9).
C
      INTEGER I, J, IWRITE, P
      REAL FX, S(3,4), STEP, TOL, X(3)
      EXTERNAL I1MACH, MNSX, QF, SMNSX
      INTEGER I1MACH
      REAL MNSX, SMNSX
C
C *** USE COMMON TO FIND NUMBER OF TIMES F(X) IS EVALUATED...
C
      INTEGER NF
      COMMON /SXCOMN/ NF
C
      DATA P/3/
C
C  ***  BODY  ***
C
C
C *** FIRST SOLVE THE PROBLEM USING SMNSX...
C
      X(1) = 2.E0
      X(2) = 3.E1
      X(3) = 9.E0
C
      NF = 0
C     *** STEP AND TOL ARE USED AS BOTH INPUT AND OUTPUT PARAMETERS,
C     *** SO WE MUST NOT PASS CONSTANTS FOR THEM.
      STEP = 1.E0
      TOL = 1.E-10
C
      FX = SMNSX(QF, P, STEP, TOL, X)
C
C *** PRINT OUT THE SOLUTION (ON THE STANDARD OUTPUT UNIT) ***
C
      IWRITE = I1MACH(2)
      WRITE(IWRITE,10) FX, TOL, STEP, X, NF
 10   FORMAT(21H SMNSX RETURNS F(X) =, E13.6,7H, TOL =, E10.3/
     1       11H AND STEP =, E10.3/7H AT X =, 3E14.6/6H AFTER, I5,
     2       21H FUNCTION EVALUATIONS)
C
C *** SOLVE THE PROBLEM AGAIN, THIS TIME USING MNSX...
C
      X(1) = 2.0E0
      X(2) = 30.0E0
      X(3) = 9.0E0
C
C
C *** CREATE INITIAL SIMPLEX...
C
      DO 30 J = 1, 4
         DO 20 I = 1, 3
            S(I,J) = X(I) - 0.5E0
 20         CONTINUE
         IF (J .LE. 3) S(J,J) = X(J) + 0.5E0
 30      CONTINUE
C
      NF = 0
      TOL = 1.E-10
C
      FX = MNSX(QF, 1000, P, P, S, TOL, X)
C
C *** PRINT OUT THE SOLUTION ***
C
      WRITE(IWRITE,40) FX, TOL, X, NF
 40   FORMAT(/20H MNSX RETURNS F(X) =, E13.6,10H AND TOL =, E10.3/
     1       7H AT X =,3E14.6/6H AFTER, I5, 21H FUNCTION EVALUATIONS)
 999  STOP
      END
      REAL FUNCTION QF(P, X)
C
C *** THIS ROUTINE COMPUTES THE OBJECTIVE FUNCTION, F(X)
C
      INTEGER P
      REAL X(P)
C
      INTEGER NF
      COMMON /SXCOMN/ NF
C
      INTEGER I
      REAL PT1, TEN, ZERO
C
      DATA PT1 /0.1E0/, TEN/1.E1/, ZERO/0.E0/
C
C
      NF = NF + 1
      QF = ZERO
      DO 10 I = 1, P
 10      QF = QF + X(I)
      QF = PT1 * QF**4
      DO 20 I = 1, P
 20      QF = QF + I*(X(I) - TEN)**2
 999  RETURN
      END
