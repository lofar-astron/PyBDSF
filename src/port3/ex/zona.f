C$TEST ZONA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE ZONA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM ZONE
C
C***********************************************************************
      EXTERNAL ROSEN
      INTEGER IWRITE, I1MACH
      REAL X(2), FNORM
      IWRITE = I1MACH(2)
C
      X(1) = -1.2
      X(2) = +1.0
C
      CALL ZONE( ROSEN, 2, X, 1.E-2, 100, FNORM )
C
      WRITE ( IWRITE, 9999 ) X(1), X(2), FNORM
 9999 FORMAT ( 1P3E15.6 )
      STOP
      END
      SUBROUTINE ROSEN ( N, X, F )
      INTEGER N
      REAL X(2), F(2)
      F(1) = 10.0* ( X(2) - X(1)**2 )
      F(2) = 1.0 - X(1)
      RETURN
      END
