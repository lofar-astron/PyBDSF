      SUBROUTINE DD7DUP(D, HDIAG, IV, LIV, LV, N, V)
C
C  ***  UPDATE SCALE VECTOR D FOR  DMNH  ***
C
C  ***  PARAMETER DECLARATIONS  ***
C
      INTEGER LIV, LV, N
      INTEGER IV(LIV)
      DOUBLE PRECISION D(N), HDIAG(N), V(LV)
C
C  ***  LOCAL VARIABLES  ***
C
      INTEGER DTOLI, D0I, I
      DOUBLE PRECISION T, VDFAC
C
C  ***  INTRINSIC FUNCTIONS  ***
C/+
      DOUBLE PRECISION DSQRT
C/
C  ***  SUBSCRIPTS FOR IV AND V  ***
C
      INTEGER DFAC, DTOL, DTYPE, NITER
C/6
C     DATA DFAC/41/, DTOL/59/, DTYPE/16/, NITER/31/
C/7
      PARAMETER (DFAC=41, DTOL=59, DTYPE=16, NITER=31)
C/
C
C-------------------------------  BODY  --------------------------------
C
      I = IV(DTYPE)
      IF (I .EQ. 1) GO TO 10
         IF (IV(NITER) .GT. 0) GO TO 999
C
 10   DTOLI = IV(DTOL)
      D0I = DTOLI + N
      VDFAC = V(DFAC)
      DO 20 I = 1, N
         T = DMAX1(DSQRT(DABS(HDIAG(I))), VDFAC*D(I))
         IF (T .LT. V(DTOLI)) T = DMAX1(V(DTOLI), V(D0I))
         D(I) = T
         DTOLI = DTOLI + 1
         D0I = D0I + 1
 20      CONTINUE
C
 999  RETURN
C  ***  LAST CARD OF DD7DUP FOLLOWS  ***
      END
