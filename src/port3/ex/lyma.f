C$TEST LYMA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LYMA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SYSS
C
C***********************************************************************
         INTEGER N, L, I, IWRITE, I1MACH
         REAL C(5000), B(100)
         REAL SUM, FLOAT, ABS, ERR, COND
         DO 40  N=10,90,40
C
C CREATE THE MATRIX A(I,J)=ABS(I-J), PACK IT INTO
C THE VECTOR C AND FORM THE RIGHT-HAND SIDE SO THE
C SOLUTION HAS ALL ONES.
C
           L=1
           SUM=(N*(N-1))/2
           DO 20 I=1,N
              DO 10 J=I,N
                 C(L)=J-I
                 L=L+1
  10          CONTINUE
              B(I)=SUM
              SUM=SUM+FLOAT(I-(N-I))
  20       CONTINUE
C
C SOLVE THE SYSTEM AND GET THE CONDITION NUMBER OF THE MATRIX
           CALL SYSS(N,C,B,100,1,COND)
C
C COMPUTE THE ERROR IN THE SOLUTION
           ERR=0.0
           DO 30 I=1,N
  30           ERR=ERR+ABS(B(I)-1.0)
           ERR=ERR/FLOAT(N)
           IWRITE=I1MACH(2)
           WRITE(IWRITE,31)N
  31       FORMAT(/8H FOR N= ,I5)
           WRITE(IWRITE,32)COND
  32       FORMAT(23H CONDITION ESTIMATE IS 1PE15.7)
           WRITE(IWRITE,33)ERR
  33       FORMAT(30H RELATIVE ERROR IN SOLUTION IS,1PE15.7)
  40     CONTINUE
         STOP
         END
