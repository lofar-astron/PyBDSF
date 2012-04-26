C$TEST LRPE
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LRPE
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM FEAS
C
C***********************************************************************
       REAL X(4),B(5),A(5,4),SIMP(8)
       INTEGER ISIMP(8)
       DATA B(1)/5.0/,B(2)/9.0/,B(3)/9.0/,B(4)/1.0/,B(5)/-5.0/
       N=4
       IA=5
       M=5
       IE=2
       IWRITE=I1MACH(2)
C
C SET UP GENERAL CONSTRAINTS
C
       DO 10 J=1,N
          A(1,J)=FLOAT(J)
          A(2,J)=FLOAT(J+1)
          A(3,J)=FLOAT(J*J)
          A(4,J)=0.0
          A(5,J)=0.0
  10   CONTINUE
       A(4,1)=1.0
       A(4,2)=1.0
       A(5,2)=-1.0
       A(5,4)=-1.0
C
C SET UP SIMPLE CONSTRAINTS
C
      IS=8
      DO 20 I=1,N
         SIMP(I)=FLOAT(-I)
         ISIMP(I)=I
         SIMP(I+N)=FLOAT(I+2)
         ISIMP(I+N)=-I
  20   CONTINUE
C
C SET UP INITIAL GUESS
C
      DO 30 I=1,N
         X(I)=1.0
  30  CONTINUE
C
C CALL FEASIBLE POINT ALGORITHM
C
      CALL FEAS(A,M,N,IA,B,X,15,IS,SIMP,ISIMP,IE)
      WRITE(IWRITE,31)(X(I),I=1,N)
 31   FORMAT(11H SOLUTION: ,4E15.6)
C
C CHECK ANSWER
C
      DO 40 I=1,M
         S = SDOT(N, A(I,1), IA, X, 1) -B(I)
         WRITE(IWRITE,41)I, S
 40   CONTINUE
 41   FORMAT(28H THE RESIDUAL AT CONSTRAINT ,I4,4H IS ,E15.5)
      STOP
      END
