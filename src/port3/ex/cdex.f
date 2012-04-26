C$TEST CDEX
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE CDEX
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM CDEXP
C
C***********************************************************************
      DOUBLE PRECISION A(2),EXPON(2)
      IWRITE = I1MACH(2)
C
      A(1) = 3.D0
      A(2) = -1.D0
      CALL CDEXP(A,EXPON)
C
      WRITE(IWRITE,9999) A, EXPON
 9999 FORMAT (18H THE EXPONENTIAL (,1PD10.4,2H, ,1PD11.4,8H) IS    //
     1           4H   (,2PD25.18,2H, ,2PD26.18,1H))
C
      STOP
      END
