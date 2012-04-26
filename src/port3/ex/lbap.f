C$TEST LBAP
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LBAP
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BAML
C
C***********************************************************************
         INTEGER IG, M, ML, N, I, IWRITE, I1MACH
         REAL G(5,20), X(20), B(20), UNI, ERR, SASUM, ABS, COND
         IG=5
         M=5
         N=10
         ML=3
C
C CONSTRUCT THE A MATRIX AND PACK IT INTO G
C
          DO 10 I=1,N
             G(1,I)=2.0
             G(2,I)=1.0
             G(3,I)=0.0
             G(4,I)=1.0
             G(5,I)=2.0
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
          CALL BAML(N,ML,M,G,IG,X,B)
C
C SOLVE THE SYSTEM AX=B
C
          CALL BASS(N,ML,M,G,IG,B,N,1,COND)
C
C PRINT OUT THE TRUE SOLUTION AND THE COMPUTED SOLUTION
C
          IWRITE=I1MACH(2)
          WRITE(IWRITE,21)
  21      FORMAT(34H TRUE SOLUTION   COMPUTED SOLUTION)
          WRITE(IWRITE,22)(X(I),B(I),I=1,N)
  22      FORMAT(1H ,2E17.8)
C
C COMPUTE THE RELATIVE ERROR
C
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
