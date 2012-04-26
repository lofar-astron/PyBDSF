C$TEST CSPG
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE CSPG
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM CSPIN
C
C***********************************************************************
C
       INTEGER IWRITE,I1MACH,I
       REAL X(9),Y(9),YY(9),XX(9)
C
C COMPUTED THE POINTS AT WHICH THE SPLINE IS TO BE FITTED
C
       DO 10 J=1,9
       X(J)=FLOAT(J-1)/8.
       Y(J)=X(J)**3
   10  CONTINUE
C
C SET THE POINTS AT WHICH INTERPOLATION IS TO BE DONE
C
       XX(1)=.3
       XX(2)=.6
       XX(3)=.9
C
C PERFORM THE INTERPOLATION
C
       CALL CSPIN(X,Y,9,XX,YY,3)
C
C SET THE OUTPUT UNIT
C
       IWRITE=I1MACH(2)
C
       WRITE (IWRITE,9998)
 9998     FORMAT(2X,2HXX,5X,11HINTERPOLATE//)
C
       WRITE (IWRITE,9999) (XX(J), YY(J), J=1,3)
 9999     FORMAT(F6.3,F12.6)
C
       STOP
       END
