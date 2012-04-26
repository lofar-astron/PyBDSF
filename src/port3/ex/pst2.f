C$TEST PST2
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PST2
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(1100)
      EXTERNAL HANDLE, BC, AF, POSTD
      INTEGER NDX, K, IS(1000), NU, NV, NMESH
      REAL ERRPAR(2), U(200), V(1), MESH(100), DT, RS(1000)
      REAL WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST  POST ON
C      U SUB T = U SUB XX + F      ON (0,1)
C  BY SETTING U1 = U AND U2 = U1 SUB X AND SOLVING
C      U1 SUB T = U1 SUB XX + F
C                                  ON (0,1)
C      U1 SUB X = U2
C WHERE F IS CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = EXP(XT).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(1100, 4)
      NU = 2
      NV = 0
      ERRPAR(1) = 0
C ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      TSTOP = 1
      DT = 1E-2
      K = 4
C NDX UNIFORM MESH POINTS ON (0,1).
      NDX = 4
      CALL UMB(0E0, 1E0, NDX, K, MESH, NMESH)
C INITIAL CONDITIONS FOR U1.
      CALL SETR(NMESH-K, 1E0, U)
C INITIAL CONDITIONS FOR U2.
      TEMP = NMESH-K
      CALL SETR(NMESH-K, 0E0, U(TEMP+1))
      CALL POST(U, NU, K, MESH, NMESH, V, NV, 0E0, TSTOP, DT, AF, BC, 
     1   POSTD, ERRPAR, HANDLE)
C CHECK FOR ERRORS AND STACK USAGE STATISTICS.
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, X, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NX
      INTEGER NV
      REAL T, X(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(NX, NU)
      REAL V(1), VT(1), A(NX, NU), AU(NX, NU, NU), AUX(NX, NU, NU), AUT(
     1   NX, NU, NU)
      REAL AUTX(NX, NU, NU), AV(1), AVT(1), F(NX, NU), FU(NX, NU, NU), 
     1   FUX(NX, NU, NU)
      REAL FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(1), FVT(1)
      INTEGER I
      REAL EXP
      DO  1 I = 1, NX
         A(I, 1) = -U(I, 2)
         AU(I, 1, 2) = -1
         F(I, 1) = (X(I)-T**2)*EXP(X(I)*T)-UT(I, 1)
         FUT(I, 1, 1) = -1
         A(I, 2) = U(I, 1)
         AU(I, 2, 1) = 1
         F(I, 2) = U(I, 2)
         FU(I, 2, 2) = 1
   1     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU
      INTEGER NV
      REAL T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      REAL UTX(NU, 2), V(1), VT(1), B(NU, 2), BU(NU, NU, 2), BUX(NU, NU,
     1   2)
      REAL BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(1), BVT(1)
      REAL EXP
      B(1, 1) = U(1, 1)-1.
      B(1, 2) = U(1, 2)-EXP(T)
      BU(1, 1, 1) = 1
      BU(1, 1, 2) = 1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV, K, X, 
     1   NX, DT, TSTOP)
      INTEGER NXMK, NU, NX
      INTEGER NV, K
      REAL T0, U0(NXMK, NU), V0(1), T, U(NXMK, NU), V(1)
      REAL X(NX), DT, TSTOP
      COMMON /TIME/ TT
      REAL TT
      EXTERNAL U1OFX, U2OFX
      INTEGER I1MACH
      REAL EU(2), EESFF
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .EQ. T) RETURN
C U1OFX AND U2OFX NEED TIME.
      TT = T
      EU(1) = EESFF(K, X, NX, U, U1OFX)
      EU(2) = EESFF(K, X, NX, U(1, 2), U2OFX)
      TEMP = I1MACH(2)
      WRITE (TEMP,  1) T, EU
   1  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 2(1PE10.2))
      RETURN
      END
      SUBROUTINE U1OFX(X, NX, U, W)
      INTEGER NX
      REAL X(NX), U(NX), W(NX)
      COMMON /TIME/ T
      REAL T
      INTEGER I
      REAL EXP
      DO  1 I = 1, NX
         U(I) = EXP(X(I)*T)
   1     CONTINUE
      RETURN
      END
      SUBROUTINE U2OFX(X, NX, U, W)
      INTEGER NX
      REAL X(NX), U(NX), W(NX)
      COMMON /TIME/ T
      REAL T
      INTEGER I
      REAL EXP
      DO  1 I = 1, NX
         U(I) = T*EXP(X(I)*T)
   1     CONTINUE
      RETURN
      END
