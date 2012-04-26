C$TEST LBAA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LBAA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BASS
C
C***********************************************************************
       INTEGER N, IG, ML, M, I, J, IWRITE, I1MACH
       REAL G(13,80), B(80,2), X(80)
       REAL START, FLOAT, ERR, ERR2, ABS, COND
       IG=13
       N=80
       DO 60 ML=2,6
C
C CONSTRUCT THE MATRIX A(I,J)=I+J AND PACK IT INTO G
C
            M=2*ML-1
            START=-FLOAT(M-ML)
            DO 20 I=1,N
               G(1,I)=START+FLOAT(2*I)
               IF(M.EQ.1) GO TO 20
               DO 10 J=2,M
                  G(J,I)=G(J-1,I)+1.
  10           CONTINUE
  20        CONTINUE
C CONSTRUCT FIRST RIGHT-HAND SIDE SO SOLUTION IS ALL 1S
            DO 30 I=1,N
  30           X(I)=1
            CALL BAML(N,ML,M,G,IG,X,B)
C CONSTRUCT THE SECOND COLUMN SO X(I)=I
            DO 40 I=1,N
  40           X(I)=I
            CALL BAML(N,ML,M,G,IG,X,B(1,2))
C SOLVE THE SYSTEM
            CALL BASS(N,ML,M,G,IG,B,80,2,COND)
C COMPUTE THE ERRORS IN THE SOLUTION
            ERR=0.0
            ERR2=0.0
            DO 50 I=1,N
               ERR=ERR+ABS(B(I,1)-1.0)
               ERR2=ERR2+ABS(B(I,2)-FLOAT(I))
  50        CONTINUE
            ERR=ERR/FLOAT(N)
            ERR2=ERR2/FLOAT(N*(N+1))*2.0
            IWRITE=I1MACH(2)
            WRITE(IWRITE,51)ML,COND
  51        FORMAT(/9H WHEN ML=,I4,21H THE CONDITION NO. IS,1PE15.7)
            WRITE(IWRITE,52)ERR
  52        FORMAT(38H REL. ERROR IN THE FIRST SOLUTION IS  ,1PE15.7)
            WRITE(IWRITE,53)ERR2
  53        FORMAT(38H REL. ERROR IN THE SECOND SOLUTION IS ,1PE15.7)
  60     CONTINUE
  70  CONTINUE
      STOP
      END
