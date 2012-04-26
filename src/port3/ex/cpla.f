C$TEST CPLA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE CPLA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM CPOLY
C
C***********************************************************************
      INTEGER IWRITE,I1MACH,K
      REAL CR(4), CI(4), ZR(3), ZI(3)
C
      CR(1) = 2.0
      CI(1) = 0.0
C
      CR(2) = -8.0
      CI(2) = 13.0
C
      CR(3) = 3.0
      CI(3) = 74.0
C
      CR(4) = 135.0
      CI(4) = 105.0
C
      CALL CPOLY(3, CR, CI, ZR, ZI)
C
      IWRITE = I1MACH(2)
      WRITE(IWRITE,99) (ZR(K),ZI(K),K = 1,3)
  99  FORMAT(1H ,1P2E15.7)
C
      STOP
      END
