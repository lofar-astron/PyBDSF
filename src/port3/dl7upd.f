      SUBROUTINE DL7UPD(BETA, GAMMA, L, LAMBDA, LPLUS, N, W, Z)
C
C  ***  COMPUTE LPLUS = SECANT UPDATE OF L  ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER N
      DOUBLE PRECISION BETA(N), GAMMA(N), L(1), LAMBDA(N), LPLUS(1),
     1                 W(N), Z(N)
C     DIMENSION L(N*(N+1)/2), LPLUS(N*(N+1)/2)
C
C--------------------------  PARAMETER USAGE  --------------------------
C
C   BETA = SCRATCH VECTOR.
C  GAMMA = SCRATCH VECTOR.
C      L (INPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE.
C LAMBDA = SCRATCH VECTOR.
C  LPLUS (OUTPUT) LOWER TRIANGULAR MATRIX, STORED ROWWISE, WHICH MAY
C             OCCUPY THE SAME STORAGE AS  L.
C      N (INPUT) LENGTH OF VECTOR PARAMETERS AND ORDER OF MATRICES.
C      W (INPUT, DESTROYED ON OUTPUT) RIGHT SINGULAR VECTOR OF RANK 1
C             CORRECTION TO  L.
C      Z (INPUT, DESTROYED ON OUTPUT) LEFT SINGULAR VECTOR OF RANK 1
C             CORRECTION TO  L.
C
C-------------------------------  NOTES  -------------------------------
C
C  ***  APPLICATION AND USAGE RESTRICTIONS  ***
C
C        THIS ROUTINE UPDATES THE CHOLESKY FACTOR  L  OF A SYMMETRIC
C     POSITIVE DEFINITE MATRIX TO WHICH A SECANT UPDATE IS BEING
C     APPLIED -- IT COMPUTES A CHOLESKY FACTOR  LPLUS  OF
C     L * (I + Z*W**T) * (I + W*Z**T) * L**T.  IT IS ASSUMED THAT  W
C     AND  Z  HAVE BEEN CHOSEN SO THAT THE UPDATED MATRIX IS STRICTLY
C     POSITIVE DEFINITE.
C
C  ***  ALGORITHM NOTES  ***
C
C        THIS CODE USES RECURRENCE 3 OF REF. 1 (WITH D(J) = 1 FOR ALL J)
C     TO COMPUTE  LPLUS  OF THE FORM  L * (I + Z*W**T) * Q,  WHERE  Q
C     IS AN ORTHOGONAL MATRIX THAT MAKES THE RESULT LOWER TRIANGULAR.
C        LPLUS MAY HAVE SOME NEGATIVE DIAGONAL ELEMENTS.
C
C  ***  REFERENCES  ***
C
C 1.  GOLDFARB, D. (1976), FACTORIZED VARIABLE METRIC METHODS FOR UNCON-
C             STRAINED OPTIMIZATION, MATH. COMPUT. 30, PP. 796-811.
C
C  ***  GENERAL  ***
C
C     CODED BY DAVID M. GAY (FALL 1979).
C     THIS SUBROUTINE WAS WRITTEN IN CONNECTION WITH RESEARCH SUPPORTED
C     BY THE NATIONAL SCIENCE FOUNDATION UNDER GRANTS MCS-7600324 AND
C     MCS-7906671.
C
C------------------------  EXTERNAL QUANTITIES  ------------------------
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      DOUBLE PRECISION DSQRT
C/
C--------------------------  LOCAL VARIABLES  --------------------------
C
      INTEGER I, IJ, J, JJ, JP1, K, NM1, NP1
      DOUBLE PRECISION A, B, BJ, ETA, GJ, LJ, LIJ, LJJ, NU, S, THETA,
     1                 WJ, ZJ
      DOUBLE PRECISION ONE, ZERO
C
C  ***  DATA INITIALIZATIONS  ***
C
C/6
C     DATA ONE/1.D+0/, ZERO/0.D+0/
C/7
      PARAMETER (ONE=1.D+0, ZERO=0.D+0)
C/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      NU = ONE
      ETA = ZERO
      IF (N .LE. 1) GO TO 30
      NM1 = N - 1
C
C  ***  TEMPORARILY STORE S(J) = SUM OVER K = J+1 TO N OF W(K)**2 IN
C  ***  LAMBDA(J).
C
      S = ZERO
      DO 10 I = 1, NM1
         J = N - I
         S = S + W(J+1)**2
         LAMBDA(J) = S
 10      CONTINUE
C
C  ***  COMPUTE LAMBDA, GAMMA, AND BETA BY GOLDFARB*S RECURRENCE 3.
C
      DO 20 J = 1, NM1
         WJ = W(J)
         A = NU*Z(J) - ETA*WJ
         THETA = ONE + A*WJ
         S = A*LAMBDA(J)
         LJ = DSQRT(THETA**2 + A*S)
         IF (THETA .GT. ZERO) LJ = -LJ
         LAMBDA(J) = LJ
         B = THETA*WJ + S
         GAMMA(J) = B * NU / LJ
         BETA(J) = (A - B*ETA) / LJ
         NU = -NU / LJ
         ETA = -(ETA + (A**2)/(THETA - LJ)) / LJ
 20      CONTINUE
 30   LAMBDA(N) = ONE + (NU*Z(N) - ETA*W(N))*W(N)
C
C  ***  UPDATE L, GRADUALLY OVERWRITING  W  AND  Z  WITH  L*W  AND  L*Z.
C
      NP1 = N + 1
      JJ = N * (N + 1) / 2
      DO 60 K = 1, N
         J = NP1 - K
         LJ = LAMBDA(J)
         LJJ = L(JJ)
         LPLUS(JJ) = LJ * LJJ
         WJ = W(J)
         W(J) = LJJ * WJ
         ZJ = Z(J)
         Z(J) = LJJ * ZJ
         IF (K .EQ. 1) GO TO 50
         BJ = BETA(J)
         GJ = GAMMA(J)
         IJ = JJ + J
         JP1 = J + 1
         DO 40 I = JP1, N
              LIJ = L(IJ)
              LPLUS(IJ) = LJ*LIJ + BJ*W(I) + GJ*Z(I)
              W(I) = W(I) + LIJ*WJ
              Z(I) = Z(I) + LIJ*ZJ
              IJ = IJ + I
 40           CONTINUE
 50      JJ = JJ - J
 60      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DL7UPD FOLLOWS  ***
      END
