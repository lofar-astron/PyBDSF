C$TEST LPSF
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSF
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPLE
C
C***********************************************************************
           INTEGER IG, N, MU, MLM1, I, KBLOK, KK, J
           INTEGER IWRITE, I1MACH
           REAL G(11,100), B(100), X(100)
           REAL ERR, AMAX1
           IG=11
           N=100
           MU=11
C
C SET UP MATRIX FOR ELLIPTIC PDE IN 2 DIMENSIONS
C
           MLM1=MU-1
           I=0
           DO 30 KBLOK=1,MLM1
              DO 20 KK=1,MLM1
                 I=I+1
                 G(1,I)=4.0
                 G(2,I)=-1.0
                 DO 10 J=3,MLM1
                    G(J,I)=0.0
   10            CONTINUE
                 G(MU,I)=-1.0
   20         CONTINUE
              G(2,I)=0.0
   30      CONTINUE
C
C SET UP RIGHT HAND SIDE SO SOLUTION IS ALL 1'S
C
           DO 40 I=1,N
              X(I)=1.0
  40       CONTINUE
           CALL BPML(N,MU,G,IG,X,B)
C
C SOLVE THE SYSTEM
C
           CALL BPLE(N,MU,G,IG,B,100,1)
C
C COMPUTE THE ERROR
C
           ERR=0.0
           DO 50 I=1,N
              ERR=AMAX1(ERR,ABS(B(I)-1.0))
  50       CONTINUE
           IWRITE=I1MACH(2)
           WRITE(IWRITE,51)ERR
  51       FORMAT(31H ERROR IN SOLUTION FROM BPLE IS,F15.8)
           STOP
           END
