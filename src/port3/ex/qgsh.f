C$TEST QGSH
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE QGSH
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GQEX
C
C***********************************************************************
      REAL X(5),W(5),CALC,TRUE,PI,ERR
C
      CALL  GQEX(5,X,W)
      IWRITE=I1MACH(2)
      WRITE(IWRITE,30)
      DO 10 J=1,5
   10    WRITE(IWRITE,40) J, X(J),W(J)
      CALC = 0.E0
      DO 20 J=1,5
   20    CALC = CALC+W(J)*X(J)/(1.0 - EXP(-X(J)))
      PI   = 2.E0*ATAN2(1.E0,0.E0)
      TRUE = PI**2/6.E0
      ERR  = TRUE - CALC
      WRITE(IWRITE,50) TRUE,CALC,ERR
      STOP
   30 FORMAT(///14H TEST OF  GQEX//30H0ABSCISSAS AND WEIGHTS FOR N=5)
   40 FORMAT(I4,0P2E16.7)
   50 FORMAT(15H0SAMPLE PROBLEM/6H TRUE=,1PE16.7/
     X   6H CALC=,1PE16.7/6H ERR =,1PE11.2)
      END
