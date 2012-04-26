C$TEST PRSJ
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PRSJ
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SPMCE
C
C***********************************************************************
      INTEGER JCOL(10000),IROW(626), IREAD, IWRITE, M
      INTEGER I1MACH, I, J, K, L, ISIZE, MP1
      INTEGER MRP(625), MCP(625), IL(625)
      REAL AROW(10000), A1, Z(625)
      COMMON /CSTAK/ D
      DOUBLE PRECISION D(3000)
      CALL ISTKIN(3000,4)
      IREAD=I1MACH(1)
      IWRITE=I1MACH(2)
  10  READ(IREAD,11)M,A1
  11  FORMAT(I3,E15.5)
      IF (M .EQ. 0) STOP
      MP1=M+1
      N= MP1*MP1-1
C
C SET UP MATRIX OF SIZE N
C
      IROW(1) = 1
      K=1
      L=0
      DO 70 I=1,MP1
         DO 60 J = 1, MP1
            L=L+1
            AROW(K) = -2.0*A1 - FLOAT(I+J-2)
            IF (J. EQ. MP1)AROW(K)=AROW(K) + A1
            JCOL(K) = L
            K=K+1
            IF (J .EQ. 1) GO TO 20
               AROW(K) = A1
               JCOL(K) = L-1
               K=K+1
  20        IF (J .EQ. MP1 .OR. J.EQ.M. AND .I.EQ.MP1) GO TO 30
               AROW(K) = J
               JCOL(K) = L+1
               K=K+1
  30        IF (I.EQ.1) GO TO 40
               AROW(K) = A1
               JCOL(K) = L-MP1
               K=K+1
  40        IF (I.EQ.MP1.OR.J.EQ.MP1. AND. I.EQ.M) GO TO 50
               AROW(K) = I
               JCOL(K) = L+MP1
               K=K+1
  50        IROW(L+1)=K
  60     CONTINUE
  70  CONTINUE
C
C REORDER ROWS OF THE MATRIX
C
      CALL SPMOR(N,IROW,JCOL,MRP,MCP)
C
C SOLVE THE SYSTEM
C
      CALL SPMCE(N,MRP,MCP,AROW,IROW,JCOL,10000,IL,ISIZE,COND,Z)
C
C PRINT RESULTS
C
      WRITE(IWRITE,71)N,IROW(N+1)
  71  FORMAT(/19HNO. OF EQUATIONS = ,I3,20H NO. OF NONZEROES = ,I5)
      WRITE(IWRITE,72)A1,ISIZE
  72  FORMAT(6H A1 = ,E15.5,9H ISIZE = , I5)
      WRITE(IWRITE,73)COND
  73  FORMAT(16H CONDITION NO = ,E15.5)
      GO TO 10
      END
C
C DATA FOR THE EXAMPLE IN THE PORT SHEET...  (REMOVE THE C
C IN COLUMN 1 BEFORE FEEDING THIS DATA TO THE PROGRAM ABOVE.)
C$DATA
C10   2.0
C20   2.0
C22   3.0
C 0   0.0
