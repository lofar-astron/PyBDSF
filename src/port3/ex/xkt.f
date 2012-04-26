C$TEST XKTH
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE XKT
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM XKTH
C
C***********************************************************************
C
      COMMON/CSTAK/DSTAK(500)
      INTEGER IWRITE, I1MACH, N, K
      REAL X(10), XK, XKTH
      REAL RSTAK(1000)
      DOUBLE PRECISION DSTAK
C
      EQUIVALENCE (DSTAK(1),RSTAK(1))
C
C  SET OUTPUT UNIT TO IWRITE .
      IWRITE = I1MACH(2)
C
      N = 8
      X(1) = 3.
      X(2) = 2.
      X(3) = 9.
      X(4) = 7.
      X(5) = 8.
      X(6) = 8.
      X(7) = 5.
      X(8) = 8.
C
      WRITE (IWRITE,98)
 98    FORMAT(1H0,14H  K       XKTH//)
C
      DO 10 K=1,8
      XK = XKTH(N,K,X)
      WRITE (IWRITE,99) K, XK
 99    FORMAT(1H ,I3,F10.1)
 10   CONTINUE
      STOP
      END
