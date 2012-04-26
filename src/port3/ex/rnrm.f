C$TEST RNRM
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE RNRM
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM RNORM
C
C***********************************************************************
C  RNORM - FIRST 10 RANDOM DEVIATES
C
      REAL X
      IWRITE = I1MACH(2)
C
      DO 10 K=1,10
      X = RNORM(0)
      WRITE (IWRITE,99) X
 99     FORMAT(1H ,F11.8)
 10   CONTINUE
C
      STOP
      END
