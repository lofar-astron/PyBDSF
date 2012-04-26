C$TEST PRST
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PRST
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SPFLU
C
C***********************************************************************
       INTEGER K, N, I1MACH, IWRITE, MAXUL, I, NEW, IREAD, IERR
       INTEGER MRP(101), MCP(101), IWORK(5000), ISIZE, NERROR
       DOUBLE PRECISION UL(5000), THRESH, EPS, GROWTH
       DOUBLE PRECISION B(101), ERROR, X
       EXTERNAL TOY
       COMMON /TOYS/ X, N, K
       MAXUL = 4000
       IREAD = I1MACH(1)
       IWRITE = I1MACH(2)
C SET THE RECOVERY MODE
       CALL ENTSRC(NEW, 1)
  10   READ(IREAD,11)K
  11   FORMAT(I2)
       IF (K .EQ. 0) STOP
       N = K*K + 1
C
       READ(IREAD,12)X, THRESH, EPS
  12   FORMAT(3D10.2)
       WRITE(IWRITE,13)K, N, X, THRESH, EPS
  13   FORMAT(3H K=,I3,3H N=,I3,3H X=,D10.2,8H THRESH=,D10.2,
     1 5H EPS=,D10.2)
C SET UP PERMUTATION VECTORS TO INDICATE NO PRIOR PIVOTING
       DO 20 I=1,N
          MRP(I) = I
          MCP(I) = I
  20   CONTINUE
       CALL DSPFLU(N, MRP, MCP, TOY, IWORK, UL, MAXUL, THRESH, EPS,
     1 ISIZE, GROWTH)
       IF (NERROR(IERR) .EQ. 0) GO TO 30
C
C TEST FOR SINGULARITY
C
          CALL ERROFF
          WRITE(IWRITE,21)
  21      FORMAT(16H SINGULAR MATRIX)
          GO TO 10
C
  30   WRITE(IWRITE,31) ISIZE, GROWTH
  31   FORMAT(7H ISIZE=,I5,8H GROWTH=,1PD25.15)
       CALL GETB(N, K, B, X)
C
C GENERATE THE RIGHT HAND SIDE AND SOLVE THE SYSTEM
C
       CALL DSPFSL(N, MRP, MCP, IWORK, UL, B, N, 1)
       ERROR = 0.0D0
C
C COMPUTE THE ERROR IN THE SOLUTION
C
       DO 40 I = 1, N
          ERROR = DMAX1(ERROR, DABS(B(I)-1.D0))
  40   CONTINUE
       WRITE(IWRITE,41)ERROR
  41   FORMAT(19H ERROR IN SOLUTION=,1PD25.15)
       GO TO 10
       END
       SUBROUTINE TOY(I, ROW, JCOL, NUM)
       INTEGER I, NUM, JCOL(101)
       INTEGER N, K, J, MODK
       DOUBLE PRECISION ROW(101)
       DOUBLE PRECISION X
       COMMON /TOYS/ X, N, K
       IF (I .LT. N) GO TO 20
C LAST ROW
          DO 10 J=1,N
             ROW(J) = 1.D0
             JCOL(J) = J
  10      CONTINUE
          NUM = N
          RETURN
  20   JCOL(1) = I
       JCOL(2) = N
       ROW(1) = 2.D0
       ROW(2) = 1.D0
       MODK = MOD(I, K)
       JCOL(3) = I-1
       ROW(3) = -1.D0
       JCOL(4) = I+1
       ROW(4) = -1.D0
       NUM = 4
       IF (MODK .GT. 1) GO TO 30
          ROW(1) = 1.D0 + X
          IF (MODK .EQ. 1) JCOL(3) = I+1
          NUM = 3
  30   IF (I .LE. K) RETURN
       IF ((I-1)/K .EQ. 1) GO TO 40
       NUM = NUM + 1
       JCOL(NUM) = I-K
       ROW(NUM) = 1.D0
  40   IF (I .GE. N-K) RETURN
       NUM = NUM + 1
       JCOL(NUM) = I+K
       ROW(NUM) = 2.D0
       RETURN
       END
