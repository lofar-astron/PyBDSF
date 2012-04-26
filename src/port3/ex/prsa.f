C$TEST PRSA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PRSA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SPMLE
C
C***********************************************************************
      INTEGER JCOL(5000),ISTAK(18000),IROW(626), IREAD, IWRITE, M
      INTEGER I1MACH, I, J, K, L, ISIZE, MP1
      REAL AROW(5000), B(625), SASUM, A1
      COMMON /CSTAK/ ISTAK
      CALL ISTKIN(18000,2)
      IREAD=I1MACH(1)
      IWRITE=I1MACH(2)
  10  READ(IREAD,11)M,A1
  11  FORMAT(I3,E15.5)
      IF (M .EQ. 0) STOP
      MP1=M+1
      N= MP1*MP1
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
            JCOL(K) = L
            K=K+1
            IF (J .EQ. 1) GO TO 20
               AROW(K) = A1
               JCOL(K) = L - 1
               K=K+1
  20        IF (J .EQ. MP1) GO TO 30
               AROW(K) = J
               JCOL(K) = L+1
               K=K+1
  30        IF (I.EQ.1) GO TO 40
               AROW(K) = A1
               JCOL(K) = L - MP1
               K=K+1
  40        IF (I.EQ.MP1) GO TO 50
               AROW(K) = I
               JCOL(K) = L+MP1
               K=K+1
  50        IROW(L+1)=K
  60     CONTINUE
  70  CONTINUE
C
C SET UP RIGHT HAND SIDE AND LAST ROW OF THE MATRIX
C
      L=IROW(N)
      DO 80 I=1,N
         AROW(L)=1.0
         JCOL(L)=I
         L=L+1
         B(I)=0.0
  80  CONTINUE
      IROW(N+1)=L
      B(N)=1.0
C
C SOLVE THE SYSTEM
C
      CALL SPMLE(N,.TRUE.,IROW,JCOL,AROW,ISIZE,B,625,1)
C
C PRINT RESULTS
C
      WRITE(IWRITE,81)N,L
  81  FORMAT(/19HNO. OF EQUATIONS = ,I3,19HNO. OF NONZEROES = ,I5)
      WRITE(IWRITE,82)ISIZE
  82  FORMAT(9H ISIZE = , I5)
      WRITE(IWRITE,83)B(N),SASUM(M,B(M),M)
  83  FORMAT(6H L1 = ,E15.5,6H L2 = ,E15.5)
      GO TO 10
      END
