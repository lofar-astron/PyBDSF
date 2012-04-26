C$TEST PST4
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PST4
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(2000)
      EXTERNAL DEE, HANDLE, BC, AF
      INTEGER NDX, K, IS(1000), NU, NV, NMESH
      REAL ERRPAR(2), U(100), V(1), ATAN, MESH(100), DT
      REAL RS(1000), WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST  POST ON
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
      TSTOP = 8.*ATAN(1E0)
      DT = 0.4
C MAKE A MESH OF NDX UNIFORM POINTS ON (-PI,+PI).
      K = 4
      NDX = 7
      CALL UMB((-4.)*ATAN(1E0), 4.*ATAN(1E0), NDX, K, MESH, NMESH)
C INITIAL CONDITIONS FOR U.
      CALL SETR(NMESH-K, 0E0, U)
C INITIAL CONDITIONS FOR V.
      V(1) = 0
      CALL POST(U, NU, K, MESH, NMESH, V, NV, 0E0, TSTOP, DT, AF, BC, 
     1   DEE, ERRPAR, HANDLE)
C CHECK FOR ERRORS AND STACK USAGE STATISTICS.
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, X, NX, U, UX, UT, UTX, NU, V, VT, NV, A, 
     1   AU, AUX, AUT, AUTX, AV, AVT, F, FU, FUX, FUT, FUTX, FV, FVT)
      INTEGER NU, NV, NX
      REAL T, X(NX), U(NX, NU), UX(NX, NU), UT(NX, NU), UTX(NX, NU)
      REAL V(NV), VT(NV), A(NX, NU), AU(NX, NU, NU), AUX(NX, NU, NU), 
     1   AUT(NX, NU, NU)
      REAL AUTX(NX, NU, NU), AV(NX, NU, NV), AVT(NX, NU, NV), F(NX, NU),
     1   FU(NX, NU, NU), FUX(NX, NU, NU)
      REAL FUT(NX, NU, NU), FUTX(NX, NU, NU), FV(NX, NU, NV), FVT(NX, 
     1   NU, NV)
      INTEGER I
      REAL COS, SIN
      DO  1 I = 1, NX
         A(I, 1) = -UX(I, 1)
         AUX(I, 1, 1) = -1
         F(I, 1) = (-UT(I, 1))-U(I, 1)**3+COS(X(I))*(COS(T)+SIN(T)+COS(X
     1      (I))**2*SIN(T)**3)
         FUT(I, 1, 1) = -1
         FU(I, 1, 1) = (-3.)*U(I, 1)**2
   1     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, L, R, U, UX, UT, UTX, NU, V, VT, NV, B, BU,
     1   BUX, BUT, BUTX, BV, BVT)
      INTEGER NU, NV
      REAL T, L, R, U(NU, 2), UX(NU, 2), UT(NU, 2)
      REAL UTX(NU, 2), V(NV), VT(NV), B(NU, 2), BU(NU, NU, 2), BUX(NU, 
     1   NU, 2)
      REAL BUT(NU, NU, 2), BUTX(NU, NU, 2), BV(NU, NV, 2), BVT(NU, NV, 2
     1   )
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
      REAL T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(NV)
      REAL D(NV), DU(NV, NXMK, NU), DUT(NV, NXMK, NU), DV(NV, NV), DVT(
     1   NV, NV)
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
      REAL T0, U0(NXMK, NU), V0(NV), T, U(NXMK, NU), V(NV)
      REAL X(NX), DT, TSTOP
      COMMON /TIME/ TT
      REAL TT
      EXTERNAL UOFX
      INTEGER I1MACH
      REAL EU, EESFF, EV
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .EQ. T) RETURN
C UOFX NEEDS TIME.
      TT = T
      EU = EESFF(K, X, NX, U, UOFX)
      EV = V(1)
      TEMP = I1MACH(2)
      WRITE (TEMP,  1) T, EU, EV
   1  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2, 6H   V =, 1P
     1   E10.2)
      RETURN
      END
      SUBROUTINE UOFX(X, NX, U, W)
      INTEGER NX
      REAL X(NX), U(NX), W(NX)
      COMMON /TIME/ T
      REAL T
      INTEGER I
      REAL COS, SIN
      DO  1 I = 1, NX
         U(I) = COS(X(I))*SIN(T)
   1     CONTINUE
      RETURN
      END
