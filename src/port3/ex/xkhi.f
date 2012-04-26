C$TEST XKHI
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE XKHI
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM IXKTH
C
C***********************************************************************
C
      COMMON/CSTAK/DSTAK(500)
      INTEGER IWRITE, N, K
      INTEGER X(10), XK
      INTEGER ISTAK(1000)
      DOUBLE PRECISION DSTAK
C
      EQUIVALENCE (DSTAK(1),ISTAK(1))
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
 98    FORMAT(1H0,15H  K       IXKTH,//)
C
      DO 10 K=1,8
      XK = IXKTH(N,K,X)
      WRITE (IWRITE,99) K, XK
 99    FORMAT(1H ,I3,I10)
 10   CONTINUE
      STOP
      END
