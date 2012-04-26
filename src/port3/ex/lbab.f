C$TEST LBAB
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LBAB
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BACE
C
C***********************************************************************
       INTEGER IG, IGL, N, ML, M, I, J, MU, IWRITE, I1MACH
       INTEGER INTER(80)
       REAL G(13, 80), B(80), X(80), GL(6, 80)
       REAL START, FLOAT, AINNO, COND, CONDNO, ABS, AINNOI
       IG=13
       IGL=6
       N=80
       IWRITE=I1MACH(2)
       DO 60 ML=2,6
C
C CONSTRUCT THE MATRIX A(I,J)=I+J AND PACK IT INTO G
            M=2*ML - 1
            START=-FLOAT(M-ML)
            DO 20 I=1,N
               G(1,I)=START+FLOAT(2*I)
               DO 10 J=2,M
                  G(J,I)=G(J-1,I)+1.
  10           CONTINUE
  20        CONTINUE
C
C DETERMINE AN ESTIMATE OF THE CONDITION NUMBER
C AND COMPUTE THE LU DECOMPOSITION
C
            CALL BACE(N,ML,M,G,IG,GL,IGL,INTER,MU,COND)
C
C DETERMINE THE NORM OF THE INVERSE MATRIX BY
C SOLVING FOR ONE COLUMN OF THE INVERSE MATRIX
C AT A TIME
C
            AINNO=0.0
            DO 50 I=1,N
C
C FIND THE ITH COLUMN OF THE INVERSE MATRIX BY
C SETTING THE RIGHT HAND SIDE TO THE ITH COLUMN
C OF THE IDENTITY MATRIX
C
                DO 30 J=1,N
                   B(J)=0.0
  30            CONTINUE
                B(I)=1.0
                CALL BAFS(N,ML,GL,IGL,INTER,B,80,1)
                CALL BABS(N,G,IG,B,80,1,MU)
C FIND THE NORM OF THE ITH COLUMN
                AINNOI=0.0
                DO 40 J=1,N
                   AINNOI=AINNOI+ABS(B(J))
  40            CONTINUE
                IF(AINNOI.GT.AINNO)AINNO=AINNOI
  50         CONTINUE
             WRITE(IWRITE,51)ML
  51         FORMAT(/6H ML IS ,I4)
             WRITE(IWRITE,52)COND
  52         FORMAT(22H CONDITION ESTIMATE IS,1PE15.7)
             CONDNO=AINNO*FLOAT(M*(N-ML+1)*2)
             WRITE(IWRITE,53)CONDNO
  53         FORMAT(22H TRUE CONDITION NO. IS,1PE15.7)
  60      CONTINUE
          STOP
          END
