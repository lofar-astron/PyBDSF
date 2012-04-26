C$TEST TTG3
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE TTG3
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM TTGR
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(350000)
      EXTERNAL HANDLE, BC, AF
      INTEGER NDX, NDY, ISTKGT, I, IUMB, IMMM
      INTEGER IS(1000), IU, IX, IY, NU, KX
      INTEGER NX, KY, NY, ILUMB
      REAL ERRPAR(2), TSTART, DT, YB(4), LX, RS(1000)
      REAL RX, WS(1000), TSTOP
      LOGICAL LS(1000)
      COMPLEX CS(500)
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO SOLVE THE LAYERED HEAT EQUATION, WITH KAPPA = 1, 1/2, 1/3,
C   DIV . ( KAPPA(X,Y) * GRAD U ) = UT + G
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(350000, 4)
      CALL ENTER(1)
      NU = 1
      LX = 0
      RX = 1
      DO  1 I = 1, 4
         YB(I) = I-1
   1     CONTINUE
      KX = 2
      KY = 2
      NDX = 3
      NDY = 3
      TSTART = 0
      TSTOP = 1
      DT = 1
      ERRPAR(1) = 1E-2
      ERRPAR(2) = 1E-4
C UNIFORM GRID.
      IX = IUMB(LX, RX, NDX, KX, NX)
C UNIFORM GRID.
      IY = ILUMB(YB, 4, NDY, KY, NY)
C MAKE MULT = KY-1.
      IY = IMMM(IY, NY, YB(2), KY-1)
C MAKE MULT = KY-1.
      IY = IMMM(IY, NY, YB(3), KY-1)
C SPACE FOR THE SOLUTION.
      IU = ISTKGT(NU*(NX-KX)*(NY-KY), 3)
      CALL SETR(NU*(NX-KX)*(NY-KY), 0E0, WS(IU))
      CALL TTGR(WS(IU), NU, KX, WS(IX), NX, KY, WS(IY), NY, TSTART, 
     1   TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, X, NX, Y, NY, NU, U, UT, UX, UY, UXT, UYT
     1   , A, AU, AUT, AUX, AUY, AUXT, AUYT, F, FU, FUT, FUX, FUY, FUXT,
     2   FUYT)
      INTEGER NU, NX, NY
      REAL T, X(NX), Y(NY), U(NX, NY, NU), UT(NX, NY, NU), UX(NX, NY, 
     1   NU)
      REAL UY(NX, NY, NU), UXT(NX, NY, NU), UYT(NX, NY, NU), A(NX, NY, 
     1   NU, 2), AU(NX, NY, NU, NU, 2), AUT(NX, NY, NU, NU, 2)
      REAL AUX(NX, NY, NU, NU, 2), AUY(NX, NY, NU, NU, 2), AUXT(NX, NY
     1   , NU, NU, 2), AUYT(NX, NY, NU, NU, 2), F(NX, NY, NU), FU(NX, 
     2   NY, NU, NU)
      REAL FUT(NX, NY, NU, NU), FUX(NX, NY, NU, NU), FUY(NX, NY, NU, NU)
     1   , FUXT(NX, NY, NU, NU), FUYT(NX, NY, NU, NU)
      INTEGER I, P, Q
      REAL KAPPA
      LOGICAL TEMP
      DO  7 I = 1, NU
         DO  6 Q = 1, NY
            DO  5 P = 1, NX
               IF (Y(Q) .GE. 1.) GOTO 1
                  KAPPA = 1
                  GOTO  4
   1              IF (Y(Q) .GE. 2.) GOTO 2
                     KAPPA = 0.5
                     GOTO  3
   2                 KAPPA = 1./3E0
   3           CONTINUE
   4           A(P, Q, I, 1) = KAPPA*UX(P, Q, I)
               AUX(P, Q, I, I, 1) = KAPPA
               A(P, Q, I, 2) = KAPPA*UY(P, Q, I)
               AUY(P, Q, I, I, 2) = KAPPA
               F(P, Q, I) = UT(P, Q, I)
               FUT(P, Q, I, I) = 1
               F(P, Q, I) = F(P, Q, I)-Y(Q)/KAPPA
               TEMP = 1. .LT. Y(Q)
               IF (TEMP) TEMP = Y(Q) .LT. 2.
               IF (TEMP) F(P, Q, I) = F(P, Q, I)+1.
               TEMP = 2. .LT. Y(Q)
               IF (TEMP) TEMP = Y(Q) .LT. 3.
               IF (TEMP) F(P, Q, I) = F(P, Q, I)+3.
   5           CONTINUE
   6        CONTINUE
   7     CONTINUE
      RETURN
      END
      SUBROUTINE BC(T, X, NX, Y, NY, LX, RX, LY, RY, U, UT, UX, 
     1   UY, UXT, UYT, NU, B, BU, BUT, BUX, BUY, BUXT, BUYT)
      INTEGER NU, NX, NY
      REAL T, X(NX), Y(NY), LX, RX, LY
      REAL RY, U(NX, NY, NU), UT(NX, NY, NU), UX(NX, NY, NU), UY(NX, NY,
     1   NU), UXT(NX, NY, NU)
      REAL UYT(NX, NY, NU), B(NX, NY, NU), BU(NX, NY, NU, NU), BUT(NX, 
     1   NY, NU, NU), BUX(NX, NY, NU, NU), BUY(NX, NY, NU, NU)
      REAL BUXT(NX, NY, NU, NU), BUYT(NX, NY, NU, NU)
      INTEGER I, J
      LOGICAL TEMP
      DO  6 J = 1, NY
         DO  5 I = 1, NX
            TEMP = X(I) .EQ. LX
            IF (.NOT. TEMP) TEMP = X(I) .EQ. RX
            IF (.NOT. TEMP) GOTO 1
               BUX(I, J, 1, 1) = 1
C LEFT OR RIGHT.
C NEUMANN BCS.
               B(I, J, 1) = UX(I, J, 1)
               GOTO  4
   1           IF (Y(J) .NE. LY) GOTO 2
                  B(I, J, 1) = U(I, J, 1)
C BOTTOM.
                  BU(I, J, 1, 1) = 1
                  GOTO  3
   2              B(I, J, 1) = U(I, J, 1)-6.*T
C TOP.
                  BU(I, J, 1, 1) = 1
   3        CONTINUE
   4        CONTINUE
   5        CONTINUE
   6     CONTINUE
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, T, U, NV, DT, TSTOP)
      INTEGER NV
      REAL T0, U0(NV), T, U(NV), DT, TSTOP
      COMMON /A7TGRP/ ERRPAR, NU, MXQ, MYQ
      INTEGER NU, MXQ, MYQ
      REAL ERRPAR(2)
      COMMON /A7TGRM/ KX, IX, NX, KY, IY, NY
      INTEGER KX, IX, NX, KY, IY, NY
      IF (T0 .NE. T) GOTO 2
         WRITE (6,  1) T
   1     FORMAT (16H RESTART FOR T =, 1PE10.2)
         RETURN
   2  CALL GERR(KX, IX, NX, KY, IY, NY, U, NU, T)
      RETURN
      END
      SUBROUTINE GERR(KX, IX, NX, KY, IY, NY, U, NU, T)
      INTEGER KX, IX, NX, KY, IY, NY
      INTEGER NU
      REAL U(1), T
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IFA, ITA(2), IXA(2), NTA(2), NXA(2), IXS
      INTEGER IYS, NXS, NYS, ISTKGT, I, IEWE
      INTEGER KA(2), MA(2), IS(1000), ILUMD, I1MACH
      REAL ABS, ERRU, AMAX1, RS(1000), WS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO GET AND PRINT THE ERROR AT EACH TIME-STEP.
C U(NX-KX,NY,KY,NU).
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
C FIND THE ERROR IN THE SOLUTION AT 2*KX * 2*KY POINTS / MESH RECTANGLE.
C X SEARCH GRID.
      IXS = ILUMD(WS(IX), NX, 2*KX, NXS)
C Y SEARCH GRID.
      IYS = ILUMD(WS(IY), NY, 2*KY, NYS)
C U SEARCH GRID VALUES.
      IEWE = ISTKGT(NXS*NYS, 3)
C THE EXACT SOLUTION.
      CALL EWE(T, WS(IXS), NXS, WS(IYS), NYS, WS(IEWE), NU)
      KA(1) = KX
      KA(2) = KY
      ITA(1) = IX
      ITA(2) = IY
      NTA(1) = NX
      NTA(2) = NY
      IXA(1) = IXS
      IXA(2) = IYS
      NXA(1) = NXS
      NXA(2) = NYS
      MA(1) = 0
C GET SOLUTION.
      MA(2) = 0
C APPROXIMATE SOLUTION VALUES.
      IFA = ISTKGT(NXS*NYS, 3)
C EVALUATE THEM.
      CALL TSD1(2, KA, WS, ITA, NTA, U, WS, IXA, NXA, MA, WS(IFA))
C ERROR IN SOLUTION VALUES.
      ERRU = 0
      TEMP = NXS*NYS
      DO  1 I = 1, TEMP
         TEMP2 = IEWE+I
         TEMP1 = IFA+I
         ERRU = AMAX1(ERRU, ABS(WS(TEMP2-1)-WS(TEMP1-1)))
   1     CONTINUE
      TEMP = I1MACH(2)
      WRITE (TEMP,  2) T, ERRU
   2  FORMAT (14H ERROR IN U(.,, 1PE10.2, 3H) =, 1PE10.2)
      CALL LEAVE
      RETURN
      END
      SUBROUTINE EWE(T, X, NX, Y, NY, U, NU)
      INTEGER NU, NX, NY
      REAL T, X(NX), Y(NY), U(NX, NY, NU)
      INTEGER I, J, P
C THE EXACT SOLUTION.
      DO  7 P = 1, NU
         DO  6 I = 1, NX
            DO  5 J = 1, NY
               IF (Y(J) .GE. 1.) GOTO 1
                  U(I, J, P) = T*Y(J)
                  GOTO  4
   1              IF (Y(J) .GE. 2.) GOTO 2
                     U(I, J, P) = 2.*T*Y(J)-T
                     GOTO  3
   2                 U(I, J, P) = 3.*T*Y(J)-3.*T
   3           CONTINUE
   4           CONTINUE
   5           CONTINUE
   6        CONTINUE
   7     CONTINUE
      RETURN
      END
