C$TEST LBAK
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LBAK
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BANM
C
C***********************************************************************
      INTEGER IG,  ML,  M,  N,  I,  J, IWRITE, I1MACH
      REAL G(13,  80),  START,  BANM,  TRNORM
      IG=13
      N=80
      DO 30 ML=2,6
C
C CONSTRUCT THE MATRIX A(I,J)=I+J AND PACK IT INTO G
C
         M=2*ML-1
         START=-FLOAT(M-ML)
         DO 20 I=1,N
            G(1,I)=START+FLOAT(2*I)
            DO 10 J=2,M
               G(J,I)=G(J-1,I)+1.0
 10         CONTINUE
 20      CONTINUE
C
C PRINT OUT THE NORM CALCULATED FROM BANM AND THE TRUE NORM
C
         TRNORM=M*(N-ML+1)*2
         IWRITE=I1MACH(2)
         WRITE(IWRITE,21)ML
 21      FORMAT(/6H ML IS ,I4)
         WRITE(IWRITE,22)TRNORM,BANM(N,ML,M,G,IG)
 22      FORMAT(15H THE TRUE NORM=,E15.5,15H COMPUTED NORM=,E15.5)
 30   CONTINUE
      STOP
      END
