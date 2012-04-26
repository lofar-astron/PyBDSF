C$TEST PST8
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PST8
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(5000)
      COMMON /TIME/ T
      REAL T
      COMMON /KMESH/ K, NMESH
      INTEGER K, NMESH
      COMMON /CMESH/ MESH
      REAL MESH(100)
      EXTERNAL DEE, HANDLE, UOFX, BC, AF
      INTEGER NDX, I, IS(1000), NU, NV
      REAL ERRPAR(2), U(100), V(100), DT, RS(1000), WS(1000)
      REAL TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO TEST  POST ON THE INTEGRO-PDE
C      U SUB T = 2 * U SUB XX - INT(0,1) EXP(X-Y)*U(Y) DY      ON (0,1)
C SUBJECT TO GIVEN DIRICHLET BCS, CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = EXP(T+X).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(5000, 4)
      NU = 1
      ERRPAR(1) = 0
C ABSOLUTE ERROR.
      ERRPAR(2) = 1E-2
      TSTOP = 1
      DT = 1E-2
      K = 4
C NDX UNIFORM MESH POINTS ON (0,1).
      NDX = 7
      CALL UMB(0E0, 1E0, NDX, K, MESH, NMESH)
      NV = NMESH-K
C UOFX NEEDS T.
      T = 0
C ICS FOR U.
      CALL L2SFF(UOFX, K, MESH, NMESH, U)
      TEMP = NMESH-K
      DO  1 I = 1, TEMP
         V(I) = U(I)
   1     CONTINUE
C ICS FOR V.
      CALL POST(U, NU, K, MESH, NMESH, V, NV, 0E0, TSTOP, DT, AF, BC, 
     1   DEE, ERRPAR, HANDLE)
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
      COMMON /KMESH/ K, NMESH
      INTEGER K, NMESH
      COMMON /CMESH/ MESH
      REAL MESH(100)
      INTEGER I
      DO  1 I = 1, NX
         A(I, 1) = 2.*UX(I, 1)
         AUX(I, 1, 1) = 2
         F(I, 1) = UT(I, 1)
         FUT(I, 1, 1) = 1
   1     CONTINUE
C GET THE INTEGRAL.
      CALL INTGRL(K, MESH, NMESH, V, X, NX, F, FV)
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
      REAL EXP
      B(1, 1) = U(1, 1)-EXP(T)
      B(1, 2) = U(1, 2)-EXP(T+1.)
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
      INTEGER I
      DO  1 I = 1, NXMK
         D(I) = U(I, 1)-V(I)
         DU(I, I, 1) = 1
         DV(I, I) = -1
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
      REAL EU, EESFF
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T0, DT
   1     FORMAT (16H RESTART FOR T =, 1PE10.2, 7H   DT =, 1PE10.2)
         RETURN
   2  TT = T
      EU = EESFF(K, X, NX, U, UOFX)
      TEMP = I1MACH(2)
      WRITE (TEMP,  3) T, EU
   3  FORMAT (14H ERROR IN U(X,, 1PE10.2, 4H ) =, 1PE10.2)
      RETURN
      END
      SUBROUTINE UOFX(X, NX, U, W)
      INTEGER NX
      REAL X(NX), U(NX), W(NX)
      COMMON /TIME/ T
      REAL T
      INTEGER I
      REAL EXP
      DO  1 I = 1, NX
         U(I) = EXP(T+X(I))
   1     CONTINUE
      RETURN
      END
      SUBROUTINE INTGRL(K, MESH, NMESH, V, X, NX, F, FV)
      INTEGER NX, NMESH
      INTEGER K
      REAL MESH(NMESH), V(1), X(NX), F(NX), FV(NX, 1)
      INTEGER MGQ, I, J, L, IX
      REAL EWE, KER, WGQ(3), XGQ(3), B(3, 4, 200), KERU
      REAL XX(3)
      LOGICAL FIRST
      INTEGER TEMP, TEMP1
      DATA FIRST/.TRUE./
C TO COMPUTE
C    F = INTEGRAL FROM MESH(1) TO MESH(NMESH)
C       KERNEL(X,Y,SUM(I=1,...,NMESH-K) V(I)*B(I,Y)) DY
C  AND
C    FV = D(F)/D(V).
C ASSUME THAT CALL KERNEL(X,Y,U,KER,KERU) RETURNS
C     KER = KERNEL(X,Y,U) AND
C     KERU = PARTIAL KERNEL / PARTIAL U.
C V(NMESH-K),FV(NX,NMESH-K)
C THE FOLLOWING DECLARATION IS SPECIFIC TO K = 4 SPLINES.
      IF (NMESH-K .GT. 200) CALL SETERR(27HINTGRL - NMESH-K .GT. NXMAX
     1   , 27, 1, 2)
C NEED MORE LOCAL SPACE.
      IF (K .NE. 4) CALL SETERR(17HINTGRL - K .NE. 4, 17, 2, 2)
C USE K-1 POINT GAUSSIAN-QUADRATURE RULE ON EACH INTERVAL.
      MGQ = K-1
      IF (FIRST) CALL GQM11(MGQ, XGQ, WGQ)
C ONLY GET GQ RULE ONCE, ITS EXPENSIVE.
C THE GAUSSIAN QUADRATURE RULE.
C DO INTEGRAL INTERVAL BY INTERVAL.
      TEMP = NMESH-K
      DO  6 I = K, TEMP
C G.Q. POINTS ON (MESH(I), MESH(I+1)).
         DO  1 J = 1, MGQ
            XX(J) = 0.5*(MESH(I+1)+MESH(I))+0.5*(MESH(I+1)-MESH(I))*XGQ(
     1         J)
   1        CONTINUE
         IF (FIRST) CALL BSPLN(K, MESH, NMESH, XX, MGQ, I, B(1, 1, I))
C ONLY GET B-SPLINE BASIS ONCE, ITS EXPENSIVE.
         DO  5 J = 1, MGQ
C GET SUM() V()*B()(XX).
            EWE = 0
            DO  2 L = 1, K
               TEMP1 = I+L-K
               EWE = EWE+V(TEMP1)*B(J, L, I)
   2           CONTINUE
            DO  4 IX = 1, NX
C GET KERNEL AND PARTIAL.
               CALL KERNEL(X(IX), XX(J), EWE, KER, KERU)
               F(IX) = F(IX)+0.5*KER*(MESH(I+1)-MESH(I))*WGQ(J)
               DO  3 L = 1, K
                  TEMP1 = I+L-K
                  FV(IX, TEMP1) = FV(IX, TEMP1)+0.5*B(J, L, I)*KERU*(
     1               MESH(I+1)-MESH(I))*WGQ(J)
   3              CONTINUE
   4           CONTINUE
   5        CONTINUE
   6     CONTINUE
      FIRST = .FALSE.
      RETURN
      END
      SUBROUTINE KERNEL(X, Y, U, KER, KERU)
      REAL X, Y, U, KER, KERU
      REAL EXP
C TO EVALUATE THE KERNEL EXP(X-Y)*U(Y) AND ITS PARTIAL WRT. U.
      KERU = EXP(X-Y)
      KER = KERU*U
      RETURN
      END
