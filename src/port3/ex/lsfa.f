C$TEST LSFA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LSFA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM L2SFF
C
C***********************************************************************
      EXTERNAL F
      INTEGER K,IWRITE,I1MACH,NT
      REAL EESFF, T(100), A(100), ERROR
C
C MAKE THE MESH
C
      K = 4
      CALL UMB (0.0E0,3.14E0,21,K,T,NT)
C
C DO THE FITTING
C
      CALL L2SFF (F, K, T, NT, A)
C
C GET THE ERROR
C
      ERROR = EESFF (K, T, NT, A, F)
C
      IWRITE = I1MACH(2)
      WRITE (IWRITE, 1000) ERROR
 1000 FORMAT (9H ERROR = ,E10.2)
C
      STOP
C
      END
      SUBROUTINE F(X, NX, FX, WX)
C
      INTEGER I,NX
      REAL X(NX), FX(NX), WX(NX)
C
      DO 1000 I = 1,NX
      FX(I) = SIN(X(I))
 1000 CONTINUE
C
      RETURN
C
      END
