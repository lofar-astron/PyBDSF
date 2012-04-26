C$TEST LRPA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LRPA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM LINPR
C
C***********************************************************************
       REAL X(4),B(3),C(4),A(3,4),SIMP(8)
       INTEGER ISIMP(8)
       N=4
       IA=3
       M=3
       IE=1
C
C SET UP GENERAL CONSTRAINTS
C
       DO 10 J=1,N
          A(1,J)=FLOAT(J)
          A(2,J)=0.0
          A(3,J)=0.0
  10   CONTINUE
       A(2,1)=1.0
       A(2,2)=1.0
       A(3,2)=-1.0
       A(3,4)=-1.0
       B(1)=5
       B(2)=1.0
       B(3)=-5.0
C
C SET UP SIMPLE CONSTRAINTS
C
      IS=8
      DO 20 I=1,N
         SIMP(I)=FLOAT(-I)
         ISIMP(I)=I
         SIMP(I+N)=10.0
         ISIMP(I+N)=-I
  20   CONTINUE
C
C SET UP COST VECTOR AND INITIAL GUESS
C
      DO 30 I=1,N
         C(I)=FLOAT(I+1)
         X(I)=1.0
  30  CONTINUE
C
C CALL LINEAR PROGRAMMING PACKAGE
C
      CALL LINPR(A,M,N,IA,B,C,X,15,CTX,IS,SIMP,ISIMP,IE)
      WRITE(6,21)(X(I),I=1,N)
 21   FORMAT(11H SOLUTION: ,4E15.6)
      WRITE(6,22)CTX
 22   FORMAT(17H FUNCTION VALUE: ,E15.5)
      STOP
      END
