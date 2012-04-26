C$TEST RPAD
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE RPAD
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM DRPOLY
C
C***********************************************************************
       INTEGER IWRITE,I1MACH,K
       DOUBLE PRECISION COEFF(6), ZR(5), ZI(5)
C
      COEFF(1) = 8.D0
      COEFF(2) = -84.D0
      COEFF(3) = 9.D0
      COEFF(4) = - 589.D0
      COEFF(5) = 331.D0
      COEFF(6) = -2915.D0
C
      CALL DRPOLY( 5, COEFF, ZR, ZI )
C
      IWRITE = I1MACH(2)
      WRITE(IWRITE,99) (ZR(K),ZI(K),K = 1,5)
  99  FORMAT(1H0,1P2E27.18)
      STOP
      END
