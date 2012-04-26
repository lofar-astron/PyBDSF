C$TEST ZONB
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE ZONB
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM ZONEJ
C
C***********************************************************************
      EXTERNAL ROSEN,MYJAC
      INTEGER IWRITE,I1MACH
      REAL X(2), FNORM
      IWRITE  =  I1MACH(2)
C
      X(1)  =  -1.2
      X(2)  =  +1.0
C
      CALL ZONEJ( ROSEN, MYJAC, 2, X, 1.E-2, 100, FNORM )
C
      WRITE ( IWRITE, 9999 ) X(1), X(2), FNORM
 9999 FORMAT ( 1P3E15.6 )
      STOP
      END
      SUBROUTINE ROSEN ( N, X, F )
      INTEGER N
      REAL X(2), F(2)
      F(1)  =  10.0 * ( X(2) - X(1)**2 )
      F(2)  =  1.0 - X(1)
      RETURN
      END
      SUBROUTINE MYJAC(ROSEN, N, X, F, DFDX, JUSED)
      EXTERNAL  ROSEN
      INTEGER N,JUSED
      REAL X(2), F(2), DFDX(2,2)
C
C  JACOBIAN OF ROSEN AT X
C
      DFDX(1,1)  =  -20.0*X(1)
      DFDX(1,2)  =   10.0
      DFDX(2,1)  =  -1.0
      DFDX(2,2)  =   0.0
      JUSED = 1
      RETURN
      END
