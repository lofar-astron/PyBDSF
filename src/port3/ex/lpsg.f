C$TEST LPSG
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSG
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPDC
C
C***********************************************************************
       INTEGER IG, MLM1, IWRITE, I1MACH, K, N, MU
       INTEGER NBLOK, KBLOK, KK, I, J, IT, ILAPSZ, IT2
       REAL G(17, 100), G2(17, 100), G3(17, 100)
       REAL COND
       IG=17
       MLM1=4
       IWRITE=I1MACH(2)
       DO 70  K=1,3
          DO 60 N=48,96,48
             MU=MLM1+1
             I=0
             NBLOK=N/MLM1
C
C SET UP THREE MATRICES FOR ELLIPTIC PDE IN 2 DIMENSION
C
             DO 30 KBLOK=1,NBLOK
                DO 20 KK=1,MLM1
                   I=I+1
                   G(1,I)=4.0
                   G(2,I)=-1.0
                   G(MU,I)=-1.0
                   DO 10 J=3,MLM1
                      G(J,I)=0.0
  10               CONTINUE
  20            CONTINUE
                G(2,I)=0.0
  30         CONTINUE
             DO 50 I=1,N
                DO 40 J=1,MU
                   G2(J,I)=G(J,I)
                   G3(J,I)=G(J,I)
  40            CONTINUE
  50         CONTINUE
             WRITE(IWRITE,51)N,MU
  51         FORMAT(/6H N IS ,I4,30H ,NUMBER OF UPPER DIAGONALS IS,I3)
C TIME DECOMPOSITION BY BPLD
             IT=ILAPSZ(0)
             CALL BPLD(N,MU,G,IG,0.0)
             IT=ILAPSZ(0)-IT
             WRITE(IWRITE,52)IT
  52         FORMAT(14H TIME FOR BPLD,I7)
C TIME DECOMPOSITION BY BPDC
             IT2=ILAPSZ(0)
             CALL BPDC(N,MU,G2,IG)
             IT2=ILAPSZ(0)-IT2
             WRITE(IWRITE,53)IT2
  53         FORMAT(14H TIME FOR BPDC,I7)
C TIME DECOMPOSITION BY BPCE
             IT3=ILAPSZ(0)
             CALL BPCE(N,MU,G3,IG,COND)
             IT3=ILAPSZ(0)-IT3
             WRITE(IWRITE,54)IT3
  54         FORMAT(14H TIME FOR BPCE,I7)
  60     CONTINUE
         MLM1=MLM1*2
  70  CONTINUE
      STOP
      END
      INTEGER FUNCTION ILAPSZ(N)
      INTEGER N
      ILAPSZ = 0
      RETURN
      END
