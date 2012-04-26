      SUBROUTINE Q7RFH(IERR, IPIVOT, N, NN, NOPIVK, P, Q, R, RLEN, W)
C
C  ***  COMPUTE QR FACTORIZATION VIA HOUSEHOLDER TRANSFORMATIONS
C  ***  WITH COLUMN PIVOTING  ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER IERR, N, NN, NOPIVK, P, RLEN
      INTEGER IPIVOT(P)
      REAL Q(NN,P), R(RLEN), W(P)
C     DIMENSION R(P*(P+1)/2)
C
C----------------------------  DESCRIPTION  ----------------------------
C
C    THIS ROUTINE COMPUTES A QR FACTORIZATION (VIA HOUSEHOLDER TRANS-
C FORMATIONS) OF THE MATRIX  A  THAT ON INPUT IS STORED IN Q.
C IF  NOPIVK  ALLOWS IT, THIS ROUTINE DOES COLUMN PIVOTING -- IF
C K .GT. NOPIVK,  THEN ORIGINAL COLUMN  K  IS ELIGIBLE FOR PIVOTING.
C THE  Q  AND  R  RETURNED ARE SUCH THAT COLUMN  I  OF  Q*R  EQUALS
C COLUMN  IPIVOT(I)  OF THE ORIGINAL MATRIX  A.  THE UPPER TRIANGULAR
C MATRIX  R  IS STORED COMPACTLY BY COLUMNS, I.E., THE OUTPUT VECTOR  R
C CONTAINS  R(1,1), R(1,2), R(2,2), R(1,3), R(2,3), ..., R(P,P) (IN
C THAT ORDER).  IF ALL GOES WELL, THEN THIS ROUTINE SETS  IERR = 0.
C BUT IF (PERMUTED) COLUMN  K  OF  A  IS LINEARLY DEPENDENT ON
C (PERMUTED) COLUMNS 1,2,...,K-1, THEN  IERR  IS SET TO  K AND THE R
C MATRIX RETURNED HAS  R(I,J) = 0  FOR  I .GE. K  AND  J .GE. K.
C    THE ORIGINAL MATRIX  A  IS AN N BY P MATRIX.  NN  IS THE LEAD
C DIMENSION OF THE ARRAY  Q  AND MUST SATISFY  NN .GE. N.  NO
C PARAMETER CHECKING IS DONE.
C    PIVOTING IS DONE AS THOUGH ALL COLUMNS OF Q WERE FIRST
C SCALED TO HAVE THE SAME NORM.  IF COLUMN K IS ELIGIBLE FOR
C PIVOTING AND ITS (SCALED) NORM**2 LOSS IS MORE THAN THE
C MINIMUM SUCH LOSS (OVER COLUMNS K THRU P), THEN COLUMN K IS
C SWAPPED WITH THE COLUMN OF LEAST NORM**2 LOSS.
C
C        CODED BY DAVID M. GAY (FALL 1979, SPRING 1984).
C
C--------------------------  LOCAL VARIABLES  --------------------------
C
      INTEGER I, II, J, K, KK, KM1, KP1, NK1
      REAL AK, QKK, S, SINGTL, T, T1, WK
      REAL  D7TPR,  R7MDC,  V2NRM
      EXTERNAL  D7TPR,  R7MDC, V2AXY,  V7SCL,  V7SCP, V7SWP,  V2NRM
C/+
      REAL  SQRT
C/
      REAL BIG, BIGRT, MEPS10, ONE, TEN, TINY, TINYRT,
     1                 WTOL, ZERO
C/6
C     DATA ONE/1.0E+0/, TEN/1.E+1/, WTOL/0.75E+0/, ZERO/0.0E+0/
C/7
      PARAMETER (ONE=1.0E+0, TEN=1.E+1, WTOL=0.75E+0, ZERO=0.0E+0)
      SAVE BIGRT, MEPS10, TINY, TINYRT
C/
      DATA BIGRT/0.0E+0/, MEPS10/0.0E+0/, TINY/0.E+0/, TINYRT/0.E+0/
C
C+++++++++++++++++++++++++++++++  BODY  ++++++++++++++++++++++++++++++++
C
      IERR = 0
      IF (MEPS10 .GT. ZERO) GO TO 10
          BIGRT =  R7MDC(5)
          MEPS10 = TEN *  R7MDC(3)
          TINYRT =  R7MDC(2)
          TINY =  R7MDC(1)
          BIG =  R7MDC(6)
          IF (TINY*BIG .LT. ONE) TINY = ONE / BIG
 10   SINGTL = FLOAT(MAX0(N,P)) * MEPS10
C
C  ***  INITIALIZE W, IPIVOT, AND DIAG(R)  ***
C
      J = 0
      DO 40 I = 1, P
         IPIVOT(I) = I
         T =  V2NRM(N, Q(1,I))
         IF (T .GT. ZERO) GO TO 20
              W(I) = ONE
              GO TO 30
 20      W(I) = ZERO
 30      J = J + I
         R(J) = T
 40      CONTINUE
C
C  ***  MAIN LOOP  ***
C
      KK = 0
      NK1 = N + 1
      DO 130 K = 1, P
         IF (NK1 .LE. 1) GO TO 999
         NK1 = NK1 - 1
         KK = KK + K
         KP1 = K + 1
         IF (K .LE. NOPIVK) GO TO 60
         IF (K .GE. P) GO TO 60
C
C        ***  FIND COLUMN WITH MINIMUM WEIGHT LOSS  ***
C
              T = W(K)
              IF (T .LE. ZERO) GO TO 60
              J = K
              DO 50 I = KP1, P
                   IF (W(I) .GE. T) GO TO 50
                        T = W(I)
                        J = I
 50                CONTINUE
              IF (J .EQ. K) GO TO 60
C
C             ***  INTERCHANGE COLUMNS K AND J  ***
C
                   I = IPIVOT(K)
                   IPIVOT(K) = IPIVOT(J)
                   IPIVOT(J) = I
                   W(J) = W(K)
                   W(K) = T
                   I = J*(J+1)/2
                   T1 = R(I)
                   R(I) = R(KK)
                   R(KK) = T1
                   CALL V7SWP(N, Q(1,K), Q(1,J))
                   IF (K .LE. 1) GO TO 60
                        I = I - J + 1
                        J = KK - K + 1
                        CALL V7SWP(K-1, R(I), R(J))
C
C        ***  COLUMN K OF Q SHOULD BE NEARLY ORTHOGONAL TO THE PREVIOUS
C        ***  COLUMNS.  NORMALIZE IT, TEST FOR SINGULARITY, AND DECIDE
C        ***  WHETHER TO REORTHOGONALIZE IT.
C
 60      AK = R(KK)
         IF (AK .LE. ZERO) GO TO 140
         WK = W(K)
C
C        *** SET T TO THE NORM OF (Q(K,K),...,Q(N,K))
C        *** AND CHECK FOR SINGULARITY.
C
         IF (WK .LT. WTOL) GO TO 70
            T =  V2NRM(NK1, Q(K,K))
            IF (T / AK .LE. SINGTL) GO TO 140
            GO TO 80
 70      T =  SQRT(ONE - WK)
         IF (T .LE. SINGTL) GO TO 140
         T = T * AK
C
C        *** DETERMINE HOUSEHOLDER TRANSFORMATION ***
C
 80      QKK = Q(K,K)
         IF (T .LE. TINYRT) GO TO 90
         IF (T .GE. BIGRT) GO TO 90
            IF (QKK .LT. ZERO) T = -T
            QKK = QKK + T
            S =  SQRT(T * QKK)
            GO TO 110
 90       S =  SQRT(T)
          IF (QKK .LT. ZERO) GO TO 100
             QKK = QKK + T
             S = S *  SQRT(QKK)
             GO TO 110
 100      T = -T
          QKK = QKK + T
          S = S *  SQRT(-QKK)
 110      Q(K,K) = QKK
C
C         ***  SCALE (Q(K,K),...,Q(N,K)) TO HAVE NORM SQRT(2)  ***
C
          IF (S .LE. TINY) GO TO 140
          CALL  V7SCL(NK1, Q(K,K), ONE/S, Q(K,K))
C
          R(KK) = -T
C
C        ***  COMPUTE R(K,I) FOR I = K+1,...,P AND UPDATE Q  ***
C
         IF (K .GE. P) GO TO 999
         J = KK + K
         II = KK
         DO 120 I = KP1, P
              II = II + I
              CALL V2AXY(NK1, Q(K,I), - D7TPR(NK1,Q(K,K),Q(K,I)),
     1                   Q(K,K), Q(K,I))
              T = Q(K,I)
              R(J) = T
              J = J + I
              T1 = R(II)
              IF (T1 .GT. ZERO)  W(I) = W(I) + (T/T1)**2
 120          CONTINUE
 130     CONTINUE
C
C  ***  SINGULAR Q  ***
C
 140  IERR = K
      KM1 = K - 1
      J = KK
      DO 150 I = K, P
         CALL  V7SCP(I-KM1, R(J), ZERO)
         J = J + I
 150     CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF Q7RFH FOLLOWS  ***
      END
