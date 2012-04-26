      SUBROUTINE Q7RAD(N, NN, P, QTR, QTRSET, RMAT, W, Y)
C
C  ***  ADD ROWS W TO QR FACTORIZATION WITH R MATRIX RMAT AND
C  ***  Q**T * RESIDUAL = QTR.  Y = NEW COMPONENTS OF RESIDUAL
C  ***  CORRESPONDING TO W.  QTR, Y REFERENCED ONLY IF QTRSET = .TRUE.
C
      LOGICAL QTRSET
      INTEGER N, NN, P
      REAL QTR(P), RMAT(1), W(NN,P), Y(N)
C     DIMENSION RMAT(P*(P+1)/2)
C/+
      REAL  SQRT
C/
      REAL  D7TPR,  R7MDC,  V2NRM
      EXTERNAL  D7TPR,  R7MDC, V2AXY,  V7SCL,  V2NRM
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER I, II, IJ, IP1, J, K, NK
      REAL ARI, QRI, RI, S, T, WI
      REAL BIG, BIGRT, ONE, TINY, TINYRT, ZERO
C/7
      SAVE BIGRT, TINY, TINYRT
C/
      DATA BIG/-1.E+0/, BIGRT/-1.E+0/, ONE/1.E+0/, TINY/0.E+0/,
     1     TINYRT/0.E+0/, ZERO/0.E+0/
C
C------------------------------ BODY -----------------------------------
C
      IF (TINY .GT. ZERO) GO TO 10
         TINY =  R7MDC(1)
         BIG =  R7MDC(6)
         IF (TINY*BIG .LT. ONE) TINY = ONE / BIG
 10   K = 1
      NK = N
      II = 0
      DO 180 I = 1, P
         II = II + I
         IP1 = I + 1
         IJ = II + I
         IF (NK .LE. 1) T =  ABS(W(K,I))
         IF (NK .GT. 1) T =  V2NRM(NK, W(K,I))
         IF (T .LT. TINY) GOTO  180
         RI = RMAT(II)
         IF (RI .NE. ZERO) GO TO 100
            IF (NK .GT. 1) GO TO 30
               IJ = II
               DO 20 J = I, P
                  RMAT(IJ) = W(K,J)
                  IJ = IJ + J
 20               CONTINUE
               IF (QTRSET) QTR(I) = Y(K)
               W(K,I) = ZERO
               GO TO 999
 30         WI = W(K,I)
            IF (BIGRT .GT. ZERO) GO TO 40
               BIGRT =  R7MDC(5)
               TINYRT =  R7MDC(2)
 40         IF (T .LE. TINYRT) GO TO 50
            IF (T .GE. BIGRT) GO TO 50
               IF (WI .LT. ZERO) T = -T
               WI = WI + T
               S =  SQRT(T * WI)
               GO TO 70
 50         S =  SQRT(T)
            IF (WI .LT. ZERO) GO TO 60
               WI = WI + T
               S = S *  SQRT(WI)
               GO TO 70
 60         T = -T
            WI = WI + T
            S = S *  SQRT(-WI)
 70         W(K,I) = WI
            CALL  V7SCL(NK, W(K,I), ONE/S, W(K,I))
            RMAT(II) = -T
            IF (.NOT. QTRSET) GO TO 80
            CALL V2AXY(NK, Y(K), - D7TPR(NK,Y(K),W(K,I)), W(K,I), Y(K))
            QTR(I) = Y(K)
 80         IF (IP1 .GT. P) GO TO 999
            DO 90 J = IP1, P
               CALL V2AXY(NK, W(K,J), - D7TPR(NK,W(K,J),W(K,I)),
     1                    W(K,I), W(K,J))
               RMAT(IJ) = W(K,J)
               IJ = IJ + J
 90            CONTINUE
            IF (NK .LE. 1) GO TO 999
            K = K + 1
            NK = NK - 1
            GO TO 180
C
 100     ARI =  ABS(RI)
         IF (ARI .GT. T) GO TO 110
            T = T *  SQRT(ONE + (ARI/T)**2)
            GO TO 120
 110     T = ARI *  SQRT(ONE + (T/ARI)**2)
 120     IF (RI .LT. ZERO) T = -T
         RI = RI + T
         RMAT(II) = -T
         S = -RI / T
         IF (NK .LE. 1) GO TO 150
         CALL  V7SCL(NK, W(K,I), ONE/RI, W(K,I))
         IF (.NOT. QTRSET) GO TO 130
            QRI = QTR(I)
            T = S * ( QRI  +   D7TPR(NK, Y(K), W(K,I)) )
            QTR(I) = QRI + T
 130     IF (IP1 .GT. P) GO TO 999
         IF (QTRSET) CALL V2AXY(NK, Y(K), T, W(K,I), Y(K))
         DO 140 J = IP1, P
            RI = RMAT(IJ)
            T = S * ( RI  +   D7TPR(NK, W(K,J), W(K,I)) )
            CALL V2AXY(NK, W(K,J), T, W(K,I), W(K,J))
            RMAT(IJ) = RI + T
            IJ = IJ + J
 140        CONTINUE
         GO TO 180
C
 150     WI = W(K,I) / RI
         W(K,I) = WI
         IF (.NOT. QTRSET) GO TO 160
            QRI = QTR(I)
            T = S * ( QRI + Y(K)*WI )
            QTR(I) = QRI + T
 160     IF (IP1 .GT. P) GO TO 999
         IF (QTRSET) Y(K) = T*WI + Y(K)
         DO 170 J = IP1, P
            RI = RMAT(IJ)
            T = S * (RI + W(K,J)*WI)
            W(K,J) = W(K,J) + T*WI
            RMAT(IJ) = RI + T
            IJ = IJ + J
 170        CONTINUE
 180     CONTINUE
C
 999  RETURN
C  ***  LAST LINE OF Q7RAD FOLLOWS  ***
      END
