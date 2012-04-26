C$TEST PST3
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PST3
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(2000)
      EXTERNAL DEE, HANDLE, BC, AF
      INTEGER NDX, K, IS(1000), NU, NV, NMESH
      REAL ERRPAR(2), U(100), V(1), MESH(100), DT, RS(1000)
      REAL WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST  POST ON
C      U SUB T = U SUB XX + V + F      ON (0,1)
C        V SUB T = U( 1/2, T )
C WHERE F IS CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = COS(XT)   AND    V(T) = 2 SIN(T/2).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(2000, 4)
      NU = 1
      NV = 1
      ERRPAR(1) = 1E-2
C ESSENTIALLY RELATIVE ERROR.
      ERRPAR(2) = 1E-6
      TSTOP = 1
      DT = 1E-6
      K = 4
      NDX = 4
C NDX UNIFORM MESH POINTS ON (0,1).
      CALL UMB(0E0, 1E0, NDX, K, MESH, NMESH)
C INITIAL CONDITIONS FOR U.
      CALL SETR(NMESH-K, 1E0, U)
C INITIAL VALUE FOR V.
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
         F(I, 1) = V(1)-UT(I, 1)-X(I)*SIN(X(I)*T)+T**2*COS(X(I)*T)-2.*
     1      SIN(T/2.)
         FUT(I, 1, 1) = -1
         FV(I, 1, 1) = 1
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
      REAL COS
      B(1, 1) = U(1, 1)-1.
      B(1, 2) = U(1, 2)-COS(T)
      BU(1, 1, 1) = 1
      BU(1, 1, 2) = 1
      RETURN
      END
      SUBROUTINE DEE(T, K, X, NX, U, UT, NU, NXMK, V, VT, NV, D, 
     1   DU, DUT, DV, DVT)
      INTEGER NXMK, NU, NV, NX
      INTEGER K
      REAL T, X(NX), U(NXMK, NU), UT(NXMK, NU), V(NV), VT(NV)
      REAL D(NV), DU(NV, NXMK, NU), DUT(NV, NXMK, NU), DV(NV, NV), DVT(
     1   NV, NV)
      INTEGER INTRVR, I, ILEFT
      REAL XI(1), BASIS(10)
      INTEGER TEMP
      XI(1) = 0.5E0
C FIND 0.5 IN MESH.
      ILEFT = INTRVR(NX, X, XI(1))
      IF (K .GT. 10) CALL SETERR(
     1   41HDEE - K .GT. 10, NEED MORE SPACE IN BASIS, 41, 1, 2)
C B-SPLINE BASIS AT XI(1).
      CALL BSPLN(K, X, NX, XI, 1, ILEFT, BASIS)
      D(1) = VT(1)
      DVT(1, 1) = 1
C VT(1) - U(0.5,T) = 0.
      DO  1 I = 1, K
         TEMP = ILEFT+I-K
         D(1) = D(1)-U(TEMP, 1)*BASIS(I)
         TEMP = ILEFT+I-K
         DU(1, TEMP, 1) = DU(1, TEMP, 1)-BASIS(I)
   1     CONTINUE
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
      REAL ABS, SIN, EU, EV, EESFF
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .EQ. T) RETURN
C UOFX NEEDS TIME.
      TT = T
      EU = EESFF(K, X, NX, U, UOFX)
      EV = ABS(V(1)-2.*SIN(T/2.))
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
      REAL COS
      DO  1 I = 1, NX
         U(I) = COS(X(I)*T)
   1     CONTINUE
      RETURN
      END
