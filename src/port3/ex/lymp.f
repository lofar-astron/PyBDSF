C$TEST LYMP
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LYMP
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SYML
C
C***********************************************************************
        INTEGER N, L, I, J, IWRITE, I1MACH
        REAL C(55), X(10), B(10)
        REAL UNI, ERR, SASUM, ABS
        N=10
C
C CONSTRUCT THE MATRIX A(I,J)=ABS(J-I) AND PACK INTO C
C
        L=0
        DO 20 I=1,N
           DO 10 J=I,N
              L=L+1
              C(L)=J-I
  10       CONTINUE
  20    CONTINUE
C
C CONSTRUCT A RANDOM VECTOR X
C
        DO 30 I=1,N
           X(I)=UNI(0)
  30    CONTINUE
C
C FIND THE VECTOR B=AX
C
       CALL SYML(N,C,X,B)
C
C SOLVE THE SYSTEM AX=B
C
       CALL SYLE(N,C,B,N,1)
C
C PRINT THE COMPUTED AND TRUE SOLUTION
C
       IWRITE=I1MACH(2)
       WRITE(IWRITE,31)
  31   FORMAT(34H TRUE SOLUTION   COMPUTED SOLUTION)
       WRITE(IWRITE,32)(X(I),B(I),I=1,N)
  32   FORMAT(1H ,2E17.8)
C
C COMPUTE THE RELATIVE ERROR
C
       ERR=0.0
       DO 40 I=1,N
          ERR=ERR+ABS(B(I)-X(I))
  40   CONTINUE
       ERR=ERR/SASUM(N,X,1)
       WRITE(IWRITE,41)ERR
  41   FORMAT(19H RELATIVE ERROR IS ,1PE15.7)
       STOP
       END
