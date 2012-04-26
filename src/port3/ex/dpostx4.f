C$TEST DPT4
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DPT4
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DPOST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(2000)
      EXTERNAL DEE, HANDLE, BC, AF
      INTEGER NDX, K, IS(1000), NU, NV, NMESH
      REAL ERRPAR(2), RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION U(100), V(1), MESH(100), DT, DATAN, WS(500)
      DOUBLE PRECISION TSTOP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST DPOST ON
C      U SUB T = U SUB XX - U**3 + F      ON (-PI,+PI)
C SUBJECT TO PERIODIC BOUNDARY CONDITIONS,
C WHERE F IS CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = COS(X)*SIN(T).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(2000, 4)
      NU = 1
      NV = 1
      ERRPAR(1) = 0
C ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      TSTOP = 8D0*DATAN(1D0)
      DT = 0.4
C MAKE A MESH OF NDX UNIFORM POINTS ON (-PI,+PI).
      K = 4
      NDX = 7
      CALL DUMB((-4D0)*DATAN(1D0), 4D0*DATAN(1D0), NDX, K, MESH, NMESH)
C INITIAL CONDITIONS FOR U.
      CALL SETD(NMESH-K, 0D0, U)
C INITIAL CONDITIONS FOR V.
      V(1) = 0
      CALL DPOST(U, NU, K, MESH, NMESH, V, NV, 0D0, TSTOP, DT, AF, BC, 
     1   DEE, ERRPAR, HANDLE)
C CHECK FOR ERRORS AND STACK USAGE STATISTICS.
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, X, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NV, NX
      DOUBLE PRECISION T, X(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(
     1   NX, NU)
      DOUBLE PRECISION V(NV), VT(NV), A(NX, NU), AU(NX, NU, NU), AUX(NX,
     1   NU, NU), AUT(NX, NU, NU)
      DOUBLE PRECISION AUTX(NX, NU, NU), AV(NX, NU, NV), AVT(NX, NU, NV)
     1   , F(NX, NU), FU(NX, NU, NU), FUX(NX, NU, NU)
      DOUBLE PRECISION FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(NX, NU, NV)
     1   , FVT(NX, NU, NV)
      INTEGER I
      DOUBLE PRECISION DCOS, DSIN
      DO  1 I = 1, NX
         A(I, 1) = -UX(I, 1)
         AUX(I, 1, 1) = -1
         F(I, 1) = (-UT(I, 1))-U(I, 1)**3+DCOS(X(I))*(DCOS(T)+DSIN(T)+
     1      DCOS(X(I))**2*DSIN(T)**3)
         FUT(I, 1, 1) = -1
         FU(I, 1, 1) = (-3D0)*U(I, 1)**2
   1     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU, NV
      DOUBLE PRECISION T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      DOUBLE PRECISION UTX(NU, 2), V(NV), VT(NV), B(NU, 2), BU(NU, NU, 2
     1   ), BUX(NU, NU, 2)
      DOUBLE PRECISION BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(NU, NV, 2), 
     1   BVT(NU, NV, 2)
      B(1, 1) = UX(1, 1)-V(1)
      B(1, 2) = UX(1, 2)-V(1)
      BUX(1, 1, 1) = 1
      BV(1, 1, 1) = -1
      BUX(1, 1, 2) = 1
      BV(1, 1, 2) = -1
      RETURN
      END
      SUBROUTINE DEE(T, K, X, NX, U, UT, NU, NXMK, V, VT, NV, D, 
     1   DU, DUT, DV, DVT)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      DOUBLE PRECISION T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(
     1   NV)
      DOUBLE PRECISION D(NV), DU(NV, NXMK, NU), DUT(NV, NXMK, NU), DV(
     1   NV, NV), DVT(NV, NV)
      INTEGER TEMP
C U(-PI,T) - U(+PI,T) = 0.
      TEMP = NX-K
      D(1) = U(1, 1)-U(TEMP, 1)
      DU(1, 1, 1) = 1
      TEMP = NX-K
      DU(1, TEMP, 1) = -1
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, V0, T, U, V, NU, NXMK, NV, K, X, 
     1   NX, DT, TSTOP)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      DOUBLE PRECISION T0, U0(NXMK, NU), V0(NV), T, U(NXMK, NU), V(NV)
      DOUBLE PRECISION X(NX), DT, TSTOP
      COMMON /TIME/ TT
      DOUBLE PRECISION TT
      EXTERNAL UOFX
      INTEGER I1MACH
      DOUBLE PRECISION DEESFF, EU, EV
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .EQ. T) RETURN
C UOFX NEEDS TIME.
      TT = T
      EU = DEESFF(K, X, NX, U, UOFX)
      EV = V(1)
      TEMP = I1MACH(2)
      WRITE (TEMP,  1) T, EU, EV
   1  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2, 6H   V =, 1P
     1   E10.2)
      RETURN
      END
      SUBROUTINE UOFX(X, NX, U, W)
      INTEGER NX
      DOUBLE PRECISION X(NX), U(NX), W(NX)
      COMMON /TIME/ T
      DOUBLE PRECISION T
      INTEGER I
      DOUBLE PRECISION DCOS, DSIN
      DO  1 I = 1, NX
         U(I) = DCOS(X(I))*DSIN(T)
   1     CONTINUE
      RETURN
      END
