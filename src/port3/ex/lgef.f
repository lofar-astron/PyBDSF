C$TEST LGEF
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGEF
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GELE
C
C***********************************************************************
        INTEGER IA, IB, I1MACH, N, I, J, IT, ILAPSZ, IWRITE
        REAL A(100, 100), AA(100, 100), B(100), BB(100)
        REAL SUM, ERR, COND, ABS, TIME, TIMES, AMAX1
        IA=100
        IB =100
C
C GENERATE THE MATRIX AND RIGHT-HAND SIDE
C
        DO 40 N=10,90,40
           DO 20 I=1,N
              SUM=0.0
              DO 10 J=1,N
                 A(I,J)=IABS(I-J)
                 IF (I.GE.J) A(I,J)=A(I,J) + 1.0
                 AA(I,J)=A(I,J)
                 SUM=SUM + AA(I,J)
  10          CONTINUE
              B(I)=SUM
              BB(I)=SUM
  20       CONTINUE
C
C CALL  GELE AND TIME IT
           IT =ILAPSZ(0)
           CALL GELE(N,A,IA,B,IB,1)
           TIME=FLOAT(ILAPSZ(0)-IT)/64.0
C
C COMPUTE THE MAXIMUM ERROR
C
           ERR=0.0
           DO 30 I=1,N
              ERR=AMAX1(ERR, ABS(B(I)-1.0))
  30       CONTINUE
C
C CALL GESS
C
           IT =ILAPSZ(0)
           CALL GESS(N,AA,IA,BB,IB,1,COND)
           TIMES=FLOAT(ILAPSZ(0)-IT)/64.0
           IWRITE=I1MACH(2)
           WRITE(IWRITE,31)N,COND
  31       FORMAT(8H FOR N= ,I4,20H CONDITION NUMBER = ,E15.7)
           WRITE(IWRITE,32)ERR
  32       FORMAT(30H MAXIMUM ERROR IN SOLUTION IS ,F15.7)
           WRITE(IWRITE,33)TIME
  33       FORMAT(34H TIME IN MILLISECONDS FOR GELE IS ,F10.2)
           WRITE(IWRITE,34)TIMES
  34       FORMAT(34H TIME IN MILLISECONDS FOR GESS IS ,F10.2)
  40    CONTINUE
        STOP
        END
      INTEGER FUNCTION ILAPSZ(N)
      INTEGER N
      ILAPSZ = 0
      RETURN
      END
