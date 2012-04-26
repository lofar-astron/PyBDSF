C$TEST LPSM
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSM
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPML
C
C***********************************************************************
         INTEGER IG, N, MU, I, IWRITE, I1MACH
         REAL G(3,20), X(20), B(20)
         REAL UNI, ERR, COND, SASUM, ABS
         IG=3
         N=10
         MU=3
C
C CONSTRUCT MATRIX A AND PACK IT INTO G
C
          DO 10 I=1,N
             G(1,I)=4.0
             G(2,I)=-1.0
             G(3,I)=-1.0
  10      CONTINUE
C
C CONSTRUCT A RANDOM VECTOR
C
          DO 20 I=1,N
             X(I)=UNI(0)
  20      CONTINUE
C
C CONSTRUCT B=AX
C
          CALL BPML(N,MU,G,IG,X,B)
C
C SOLVE THE SYSTEM AX=B
C
          CALL BPSS(N,MU,G,IG,B,N,1,COND)
C
C PRINT OUT THE TRUE SOLUTION AND THE COMPUTED SOLUTION
C
          IWRITE=I1MACH(2)
          WRITE(IWRITE,21)
  21      FORMAT(34H TRUE SOLUTION   COMPUTED SOLUTION)
          WRITE(IWRITE,22)(X(I),B(I),I=1,N)
  22      FORMAT(1H ,2E16.8)
          ERR=0.0
          DO 30 I=1,N
              ERR=ERR+ABS(B(I)-X(I))
  30      CONTINUE
          ERR=ERR/SASUM(N,X,1)
          WRITE(IWRITE,31)ERR
  31      FORMAT(19H RELATIVE ERROR IS ,1PE15.7)
          WRITE(IWRITE,32)COND
  32      FORMAT(20H CONDITION NUMBER IS,1PE15.7)
          STOP
          END
