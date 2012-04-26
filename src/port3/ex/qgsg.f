C$TEST QGSG
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE QGSG
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GQ1
C
C***********************************************************************
      REAL X(5),W(5),CALC,TRUE,ERR
C
      CALL  GQ1(5,X,W)
      IWRITE=I1MACH(2)
      WRITE(IWRITE,30)
      DO 10 J=1,5
   10    WRITE(IWRITE,40) J, X(J),W(J)
      CALC = 0.E0
      DO 20 J=1,5
   20    CALC = CALC+W(J)*(1.0/(2.0+X(J)))
      TRUE = ALOG(3.E0)
      ERR  = TRUE-CALC
      WRITE(IWRITE,50) TRUE,CALC,ERR
      STOP
   30 FORMAT(///13H TEST OF  GQ1//30H0ABSCISSAS AND WEIGHTS FOR N=5)
   40 FORMAT(I4,0P2E16.7)
   50 FORMAT(15H0SAMPLE PROBLEM/6H TRUE=,1PE16.7/
     X   6H CALC=,1PE16.7/6H ERR =,1PE11.2)
      END
