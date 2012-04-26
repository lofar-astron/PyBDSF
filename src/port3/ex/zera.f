C$TEST ZERA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE ZERA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM ZERO
C
C***********************************************************************
      EXTERNAL F
      INTEGER IWRITE,I1MACH
      REAL A,B,F,T,X,ZERO
C
      IWRITE = I1MACH(2)
      A = 1.0
      B = 3.0
      T = 1.0E-7
      X=ZERO(F,A,B,T)
C
      WRITE (IWRITE,9999) X
 9999 FORMAT (17H THE ROOT IS X = ,1PE15.8)
C
      STOP
      END
C
      REAL FUNCTION F(X)
      REAL X
      F=X*X - 4.
      RETURN
      END
