C$TEST LPSB
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LPSB
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM BPCE
C
C***********************************************************************
        INTEGER N, MU, IG, K, I, IWRITE, I1MACH, J
        REAL G(2,100), B(200)
        REAL X, COND, AINVNO, AINORM, ABS
        N=100
        X=1.0
        MU=2
        IG=2
        DO 50 K=1,3
C CONSTRUCT MATRIX
            DO 10 I=1,N
               G(1,I)=2.0
               G(2,I)=-1.0
  10        CONTINUE
            G(1,1)=1.0+X
            G(1,N)=1.0+X
C GET ESTIMATE OF CONDITION NUMBER FROM BPCE
            CALL BPCE(N,MU,G,IG,COND)
            IWRITE=I1MACH(2)
            WRITE(IWRITE,11)X
  11        FORMAT(/10H WHEN X IS,E14.6)
            WRITE(IWRITE,12)COND
  12        FORMAT(25H CONDITION ESTIMATE IS   ,E15.8)
C SINCE CONDITION NUMBER IS NORM(A)*NORM(INVERSE(A)),
C FIND THE NORM OF EACH COLUMN OF INVERSE(A). GENERATE
C THE COLUMNS ONE AT A TIME AND REUSE SPACE
            AINVNO=0.0
            DO 40 I=1,N
C GENERATE ITH COLUMN OF IDENTITY MATRIX AS RIGHT HAND SIDE
                DO 20 J=1,N
                   B(J)=0.0
  20            CONTINUE
                B(I)=1.0
C SOLVE AX=B TO GET ITH COLUMN OF A(INVERSE)
                CALL BPFS(N,MU,G,IG,B,N,1)
                CALL BPBS(N,MU,G,IG,B,N,1)
C FIND NORM OF COLUMN
                AINORM=0.0
                DO 30 J=1,N
                   AINORM=AINORM+ABS(B(J))
  30            CONTINUE
                IF(AINVNO.LT.AINORM)AINVNO=AINORM
  40        CONTINUE
            COND=4.0*AINVNO
            WRITE(IWRITE,41)COND
  41        FORMAT(25H TRUE CONDITION NUMBER IS,E15.8)
            X=X/100.0
  50     CONTINUE
         STOP
         END
