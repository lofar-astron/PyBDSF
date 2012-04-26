C$TEST DPT9
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DPT9
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DPOST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(2000)
      COMMON /PARAM/ C
      DOUBLE PRECISION C
      EXTERNAL HANDLE, DPOSTD, BC, AF
      INTEGER NDX, NXC, NXX, I, K, IS(1000)
      INTEGER NU, NV, NX, I1MACH
      REAL ERRPAR(2), RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DEEBSF, EWE(1000), ERR, U(100), V(1), X(100)
      DOUBLE PRECISION ERRR, DT, XC(100), UC(100), WS(500), XX(1000)
      DOUBLE PRECISION D1MACH, TSTOP
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST DPOST ON AUTOMATIC, STATIC MESH REFINEMENT.
C      U SUB T = U SUB XX + C * U SUB X      ON (0,1)
C THE SOLUTION IS
C      U(X,T) = EXP(-C*X).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(2000, 4)
      C = 50
      NU = 1
      NV = 0
      ERRPAR(1) = 1E-1
      ERRPAR(2) = 1E-1
      K = 4
      NDX = 8
      CALL DUMB(0D0, 1D0, NDX, K, XC, NXC)
C INITIAL CONDITIONS FOR UC.
      CALL SETD(NXC-K, 0D0, UC)
C INFINITY.
      ERR = D1MACH(2)
   1  IF (ERR .LE. 1D-2) GOTO  6
C HALVE THE CRUDE X.
         CALL DLUMB(XC, NXC, 3, K, X, NX)
C FITTING POINTS FOR REFINEMENT.
         CALL DLUMD(X, NX, K, XX, NXX)
C UC ON XX.
         CALL DSPLNE(K, XC, NXC, UC, XX, NXX, EWE)
C FIT U TO UC ON MESH.
         CALL DDL2SF(XX, EWE, NXX, K, X, NX, U)
         TSTOP = 1D0/D1MACH(4)
         DT = 1D-6
         I = NX-2*(K-1)
         TEMP = I1MACH(2)
         WRITE (TEMP,  2) I
   2     FORMAT (18H SOLVING FOR NDX =, I3)
         CALL DPOST(U, NU, K, X, NX, V, NV, 0D0, TSTOP, DT, AF, BC, 
     1      DPOSTD, ERRPAR, HANDLE)
C GET RUN-TIME STATISTICS.
         CALL DPOSTX
C ERROR ESTIMATE FOR UC.
         ERR = DEEBSF(K, XC, NXC, UC, X, NX, U)
C ERROR ESTIMATE FOR U.
         ERRR = ERR/16D0
         TEMP = I1MACH(2)
         WRITE (TEMP,  3) ERR, ERRR
   3     FORMAT (21H ERROR ESTIMATES UC =, 1PE10.2, 9H  AND U =, 1P
     1      E10.2)
         NXC = NX
         DO  4 I = 1, NX
            XC(I) = X(I)
   4        CONTINUE
         TEMP = NX-K
         DO  5 I = 1, TEMP
            UC(I) = U(I)
   5        CONTINUE
         GOTO  1
   6  STOP 
      END
      SUBROUTINE AF(T, X, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NX
      INTEGER NV
      DOUBLE PRECISION T, X(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(
     1   NX, NU)
      DOUBLE PRECISION V(1), VT(1), A(NX, NU), AU(NX, NU, NU), AUX(NX, 
     1   NU, NU), AUT(NX, NU, NU)
      DOUBLE PRECISION AUTX(NX, NU, NU), AV(1), AVT(1), F(NX, NU), FU(
     1   NX, NU, NU), FUX(NX, NU, NU)
      DOUBLE PRECISION FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(1), FVT(1)
      COMMON /PARAM/ C
      DOUBLE PRECISION C
      INTEGER I
      DO  1 I = 1, NX
         A(I, 1) = UX(I, 1)+C*U(I, 1)
         AUX(I, 1, 1) = 1
         AU(I, 1, 1) = C
         F(I, 1) = UT(I, 1)
         FUT(I, 1, 1) = 1
   1     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU
      INTEGER NV
      DOUBLE PRECISION T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      DOUBLE PRECISION UTX(NU, 2), V(1), VT(1), B(NU, 2), BU(NU, NU, 2),
     1   BUX(NU, NU, 2)
      DOUBLE PRECISION BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(1), BVT(1)
      COMMON /PARAM/ C
      DOUBLE PRECISION C
      DOUBLE PRECISION DEXP
      B(1, 1) = U(1, 1)-1D0
      B(1, 2) = U(1, 2)-DEXP(-C)
      BU(1, 1, 1) = 1
      BU(1, 1, 2) = 1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV, K, X, 
     1   NX, DT, TSTOP)
      INTEGER NXMK, NU, NX
      INTEGER NV, K
      DOUBLE PRECISION T0, U0(NXMK, NU), V0(1), T, U(NXMK, NU), V(1)
      DOUBLE PRECISION X(NX), DT, TSTOP
      COMMON /TIME/ TT
      DOUBLE PRECISION TT
      EXTERNAL UOFX
      INTEGER I1MACH
      DOUBLE PRECISION DEESFF, EU
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
         RETURN
   2  TT = T
      EU = DEESFF(K, X, NX, U, UOFX)
      TEMP = I1MACH(2)
      WRITE (TEMP,  3) T, EU
   3  FORMAT (15H ERROR IN U(X, , 1PE10.2, 4H ) =, 1PE10.2)
      RETURN
      END
      SUBROUTINE UOFX(X, NX, U, W)
      INTEGER NX
      DOUBLE PRECISION X(NX), U(NX), W(NX)
      COMMON /PARAM/ C
      DOUBLE PRECISION C
      COMMON /TIME/ T
      DOUBLE PRECISION T
      INTEGER I
      DOUBLE PRECISION DEXP
      DO  1 I = 1, NX
         U(I) = DEXP((-C)*X(I))
   1     CONTINUE
      RETURN
      END
