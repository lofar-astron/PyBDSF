C$TEST PRS1
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PRS1
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SPMML
C
C***********************************************************************
      INTEGER I, IWRITE, I1MACH, NEQ, J, NX, L
      INTEGER IROW(101), JCOL(500)
      REAL A(500)
      REAL X(100), B(100), ERR, SASUM
      REAL RSTACK(1800)
      COMMON/CSTAK/ RSTACK
      CALL ISTKIN(1800,3)
      NX = 10
C
C CONSTRUCT THE MATRIX
C
      L=1
      NEQ=1
      DO 20 I=1, NX
         DO 10 J= 1,NX
            IROW(NEQ)=L
            JCOL(L)=NEQ
            A(L)=-4.0
            L=L+1
            JCOL(L)=NEQ-1
            A(L)=1.0
            IF (J.GT.1)L=L+1
            JCOL(L)=NEQ+1
            A(L)=1.0
            IF (J.LT.NX)L=L+1
            JCOL(L)=NEQ-NX
            A(L)=1.0
            IF (I.GT.1)L=L+1
            JCOL(L)=NEQ+NX
            A(L)=1.0
            IF(I.LT.NX)L=L+1
            NEQ=NEQ+1
  10     CONTINUE
  20  CONTINUE
      IROW(NEQ)=L
      NEQ=NEQ-1
C
C CONSTRUCT A RANDOM VECTOR FOR X
C
      DO 30 I=1,NEQ
         X(I)=UNI(0)
  30  CONTINUE
C
C FIND THE VECTOR B=AX
C
      CALL SPMML(NEQ,IROW,JCOL,A,X,B)
C
C SOLVE THE SYSTEM AX=B
C
      CALL SPMLE(NEQ,.TRUE.,IROW,JCOL,A,ISIZE,B,NEQ,1)
C
C FIND THE NORM OF THE ERROR OF THE SOLUTION
C
      ERR=0.0
      IWRITE = I1MACH(2)
      DO 40 I=1,NEQ
         ERR=ERR + ABS(B(I)-X(I))
  40  CONTINUE
      ERR=ERR/SASUM(NEQ,X,1)
      WRITE(IWRITE,41)ERR
  41  FORMAT(19H RELATIVE ERROR IS ,1PE15.7)
      STOP
      END
