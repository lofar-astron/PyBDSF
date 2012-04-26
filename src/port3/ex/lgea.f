C$TEST LGEA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGEA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GESS
C
C***********************************************************************
       INTEGER N, IREAD, I1MACH, I, NB, IWRITE, J
       REAL A(5,5), B(5,2), COND
       N=5
       IREAD=I1MACH(1)
C
       DO 10 I=1,N
           READ(IREAD,1) (A(I,J),J=1,N)
    1      FORMAT(1X,5F10.0)
   10  CONTINUE
C
       NB=2
       DO 20 I=1,N
           READ(IREAD,11) (B(I,J),J=1,NB)
   11      FORMAT(1X,2F10.3)
   20  CONTINUE
C
C SOLVE AX = B  BY CALLING GESS
C
       CALL GESS(N,A,N,B,N,NB,COND)
       IWRITE=I1MACH(2)
       WRITE(IWRITE,21) COND
   21  FORMAT(52H AN ESTIMATE OF THE CONDITION NUMBER OF THE MATRIX =,
     1          E14.7)
C
       WRITE(IWRITE,22)
   22  FORMAT(27H THE COMPUTED SOLUTION X IS,//)
       DO 30 I=1,N
           WRITE(IWRITE,23) (B(I,J),J=1,NB)
   23      FORMAT(1H ,5F20.7)
   30  CONTINUE
C
       STOP
       END
C
C DATA FOR THE EXAMPLE IN THE PORT SHEET...  (REMOVE THE C
C IN COLUMN 1 BEFORE FEEDING THIS DATA TO THE PROGRAM ABOVE.)
C$DATA
C        1.       -2.        3.        7.       -9.
C       -2.        8.       -6.        9.       50.
C       11.       -6.       18.      -15.      -18.
C        7.        2.      -15.      273.      173.
C       -9.       50.      -18.        6.     1667.
C       30.    29.419
C     -191.  -190.994
C      133.   133.072
C     -986.  -985.775
C    -6496. -6495.553
