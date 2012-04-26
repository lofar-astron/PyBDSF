C$TEST LPSK
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSK
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPFS
C
C***********************************************************************
       INTEGER N, ML, IG, NM1, K, I, IWRITE, I1MACH, IT, IEND, ITER
       REAL G(2,100), B(200), R(200)
       REAL X, ERR, AMAX1, RNORM, BNORM, R1MACH, ABS
       DOUBLE PRECISION DBLE
C CONSTRUCT MATRIX AND RIGHT HAND SIDE SO TRUE SOLUTION IS
C COMPOSED ENTIRELY OF 1S
       N=100
       X=1
       ML=2
       IG=2
       NM1=N-1
       DO 90 K=1,3
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
          CALL BPLE(N,ML,G,IG,B,N,1)
          IWRITE=I1MACH(2)
          WRITE(IWRITE,11)X
  11      FORMAT(/5H X IS,F16.8)
C COMPUTE THE ERROR
          ERR=0.0
          DO 20 I=1,N
             ERR=AMAX1(ERR,ABS(B(I)-1.0))
  20      CONTINUE
          WRITE(IWRITE,21)ERR
  21      FORMAT(22H FOR BPLE THE ERROR IS,F16.8)
          IEND=I1MACH(11)*IFIX(R1MACH(5)/ALOG10(2.0)+1.0)
C FIND THE NORM OF THE SOLUTION
          BNORM=0.0
          DO 30 I=1,N
              BNORM=AMAX1(BNORM,ABS(B(I)))
  30      CONTINUE
C REFINE THE SOLUTION
          DO 60 ITER=1,IEND
             IT=ITER
C COMPUTE THE RESIDUAL R=B-AX, IN DOUBLE PRECISION
             DO 40 I=2,NM1
                R(I)=DBLE(B(I-1))+DBLE(B(I+1))-2.D0*DBLE(B(I))
  40         CONTINUE
             R(1)=X+B(2)-DBLE(1.0+X)*DBLE(B(1))
             R(N)=X+B(N-1)-DBLE(1.+X)*DBLE(B(N))
C SOLVE A(DELTAX)=R
             CALL BPFS(N,ML,G,IG,R,N,1)
             CALL BPBS(N,ML,G,IG,R,N,1)
C DETERMINE NORM OF CORRECTION AND ADD IN CORRECTION
             RNORM=0.0
             DO 50 I=1,N
                B(I)=B(I)+R(I)
                RNORM=RNORM+ABS(R(I))
  50         CONTINUE
             IF(RNORM.LT.R1MACH(4)*BNORM) GO TO 70
  60      CONTINUE
          WRITE(IWRITE,61)
  61      FORMAT(18H REFINEMENT FAILED)
C COMPUTE NEW ERROR
  70      ERR=0.0
          DO 80 I=1,N
             ERR=AMAX1(ERR,ABS(B(I)-1.0))
  80      CONTINUE
          WRITE(IWRITE,81)IT,ERR
  81      FORMAT(24H ERROR AFTER REFINEMENT ,I4,3H IS,E14.7)
          X=X/100.0
  90   CONTINUE
       STOP
       END
