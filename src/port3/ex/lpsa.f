C$TEST LPSA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPSS
C
C***********************************************************************
       INTEGER N, K, I, IWRITE, I1MACH, MU
       REAL G(2,100), B(200)
       REAL X, COND, ERR, AMAX1
C CONSTRUCT MATRIX AND RIGHT-HAND SIDE SO TRUE SOLUTION IS
C COMPOSED ENTIRELY OF ONES
       N=100
       X=1
       DO 30  K=1,3
          DO 10 I=1,N
             G(1,I)=2.0
             G(2,I)=-1.0
             B(I)=0.0
  10      CONTINUE
          G(1,1)=1.0+X
          G(1,N)=1.0+X
          B(1)=X
          B(N)=X
C SOLVE THE SYSTEM
          MU=2
          CALL BPSS(N,MU,G,2,B,N,1,COND)
          IWRITE=I1MACH(2)
          WRITE(IWRITE,11)X
  11      FORMAT(/5H X IS,F15.7)
          WRITE(IWRITE,12)COND
  12      FORMAT(20H CONDITION NUMBER IS,1PE15.7)
C COMPUTE THE ERROR
          ERR=0.0
          DO 20 I=1,N
             ERR=AMAX1(ERR,ABS(B(I)-1.0))
  20      CONTINUE
          WRITE(IWRITE,21)ERR
  21      FORMAT(22H FOR BPSS THE ERROR IS,F16.8)
          X=X/100.
  30   CONTINUE
       STOP
       END
