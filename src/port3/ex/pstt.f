C$TEST PSTT
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PSTT
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM POST
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(2000)
      EXTERNAL HANDLE, BC, AF, POSTD
      INTEGER NDX, NXH, I, K, IS(1000), NU
      INTEGER NV, NX, I1MACH
      REAL ABS, ERR, ERRPAR(2), U(100), V(1), X(100)
      REAL AMAX1, DT, UE(100), EEBSF, UH(100), XH(100)
      REAL RS(1000), WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      INTEGER TEMP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO ESTIMATE X AND T ERROR AS SUM.
C      U SUB T = U SUB XX + F      ON (0,1)
C WHERE F IS CHOSEN SO THAT THE SOLUTION IS
C      U(X,T) = EXP(XT).
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(2000, 4)
      NU = 1
      NV = 0
      ERRPAR(1) = 0
      ERRPAR(2) = 1E-2
      K = 4
      NDX = 4
      TSTOP = 1
      DT = 1E-2
C CRUDE MESH.
      CALL UMB(0E0, 1E0, NDX, K, X, NX)
C INITIAL CONDITIONS FOR U.
      CALL SETR(NX-K, 1E0, U)
      TEMP = I1MACH(2)
      WRITE (TEMP,  1) 
   1  FORMAT (36H SOLVING ON CRUDE MESH USING ERRPAR.)
      CALL POST(U, NU, K, X, NX, V, NV, 0E0, TSTOP, DT, AF, BC, POSTD, 
     1   ERRPAR, HANDLE)
C GET RUN-TIME STATISTICS.
      CALL POSTX
C HALVE THE MESH SPACING.
      CALL UMB(0E0, 1E0, 2*NDX-1, K, XH, NXH)
C INITIAL CONDITIONS FOR UH.
      CALL SETR(NXH-K, 1E0, UH)
      DT = 1E-2
      TEMP = I1MACH(2)
      WRITE (TEMP,  2) 
   2  FORMAT (38H SOLVING ON REFINED MESH USING ERRPAR.)
      CALL POST(UH, NU, K, XH, NXH, V, NV, 0E0, TSTOP, DT, AF, BC, 
     1   POSTD, ERRPAR, HANDLE)
C GET RUN-TIME STATISTICS.
      CALL POSTX
C ESTIMATE U ERROR.
      ERR = EEBSF(K, X, NX, U, XH, NXH, UH)
      WRITE (6,  3) ERR
   3  FORMAT (24H U ERROR FROM U AND UH =, 1PE10.2)
C INITIAL CONDITIONS FOR UE.
      CALL SETR(NX-K, 1E0, UE)
      DT = 1E-2
      ERRPAR(1) = ERRPAR(1)/10.
      ERRPAR(2) = ERRPAR(2)/10.
      TEMP = I1MACH(2)
      WRITE (TEMP,  4) 
   4  FORMAT (39H SOLVING ON CRUDE MESH USING ERRPAR/10.)
      CALL POST(UE, NU, K, X, NX, V, NV, 0E0, TSTOP, DT, AF, BC, POSTD
     1   , ERRPAR, HANDLE)
C GET RUN-TIME STATISTICS.
      CALL POSTX
      ERR = 0
      TEMP = NX-K
      DO  5 I = 1, TEMP
         ERR = AMAX1(ERR, ABS(U(I)-UE(I)))
   5     CONTINUE
      WRITE (6,  6) ERR
   6  FORMAT (24H U ERROR FROM U AND UE =, 1PE10.2)
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
         A(I, 1) = -UX(I, 1)
         AUX(I, 1, 1) = -1
         F(I, 1) = (X(I)-T**2)*EXP(X(I)*T)-UT(I, 1)
         FUT(I, 1, 1) = -1
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
      EXTERNAL UOFX
      INTEGER I1MACH
      REAL EU, EESFF
      INTEGER TEMP
C OUTPUT AND CHECKING ROUTINE.
      IF (T0 .NE. T) GOTO 2
         TEMP = I1MACH(2)
         WRITE (TEMP,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
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
         U(I) = EXP(X(I)*T)
   1     CONTINUE
      RETURN
      END
