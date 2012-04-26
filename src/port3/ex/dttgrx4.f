C$TEST DTG4
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DTG4
C***********************************************************************
C
C  EXAMPLE OF USE OF PORT PROGRAM DTTGR
C
C***********************************************************************
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(350000)
      EXTERNAL HANDLE, BC, AF
      INTEGER NDX, NDY, ISTKGT, IS(1000), IU, IX
      INTEGER IY, NU, KX, NX, KY, NY
      INTEGER IDUMB
      REAL ERRPAR(2), RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION TSTART, DT, LX, LY, RX, RY
      DOUBLE PRECISION WS(500), TSTOP
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO SOLVE THE LINEAR HEAT EQUATION
C   GRAD . ( UX - 0.1 * UY , 0.1*UX +  UY ) = UT - X*Y
C WITH SOLUTION U == T*X*Y ON [0,+1]**2, EXACT FOR K = 4,
C WITH TILTED TOP AND BOTTOM, NORMAL BCS THERE.
C THE PORT LIBRARY STACK AND ITS ALIASES.
C INITIALIZE THE PORT LIBRARY STACK LENGTH.
      CALL ISTKIN(350000, 4)
      CALL ENTER(1)
      NU = 1
      LX = 0
      RX = 1
      LY = 0
      RY = 1
      KX = 4
      KY = 4
      NDX = 3
      NDY = 3
      TSTART = 0
      TSTOP = 1
      DT = 1
      ERRPAR(1) = 1E-2
      ERRPAR(2) = 1E-4
C UNIFORM GRID.
      IX = IDUMB(LX, RX, NDX, KX, NX)
C UNIFORM GRID.
      IY = IDUMB(LY, RY, NDY, KY, NY)
C SPACE FOR THE SOLUTION.
      IU = ISTKGT(NU*(NX-KX)*(NY-KY), 4)
      CALL SETD(NU*(NX-KX)*(NY-KY), 0D0, WS(IU))
      CALL DTTGR(WS(IU), NU, KX, WS(IX), NX, KY, WS(IY), NY, TSTART, 
     1   TSTOP, DT, AF, BC, ERRPAR, HANDLE)
      CALL LEAVE
      CALL WRAPUP
      STOP 
      END
      SUBROUTINE AF(T, XI, NX, YI, NY, NU, U, UT, UX, UY, UXT, 
     1   UYT, A, AU, AUT, AUX, AUY, AUXT, AUYT, F, FU, FUT, FUX, FUY, 
     2   FUXT, FUYT)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, XI(NX), YI(NY), U(NX, NY, NU), UT(NX, NY, NU),
     1   UX(NX, NY, NU)
      DOUBLE PRECISION UY(NX, NY, NU), UXT(NX, NY, NU), UYT(NX, NY, NU),
     1   A(NX, NY, NU, 2), AU(NX, NY, NU, NU, 2), AUT(NX, NY, NU, NU, 2)
      DOUBLE PRECISION AUX(NX, NY, NU, NU, 2), AUY(NX, NY, NU, NU, 2), 
     1   AUXT(NX, NY, NU, NU, 2), AUYT(NX, NY, NU, NU, 2), F(NX, NY, NU)
     2   , FU(NX, NY, NU, NU)
      DOUBLE PRECISION FUT(NX, NY, NU, NU), FUX(NX, NY, NU, NU), FUY(NX,
     1   NY, NU, NU), FUXT(NX, NY, NU, NU), FUYT(NX, NY, NU, NU)
      EXTERNAL BT, LR
      INTEGER I, P, Q
      DOUBLE PRECISION D(600), X, Y, XX(100), YY(100)
      INTEGER TEMP
      IF (NX*NY .GT. 100) CALL SETERR(19HAF - NX*NY .GT. 100, 19, 1, 2)
      CALL DBTMAP(T, XI, YI, NX, NY, LR, BT, XX, YY, D)
C MAP INTO (X,Y).
      CALL DTTGRU(NX, NY, D, UX, UY, UT, NU)
      DO  3 I = 1, NU
         DO  2 Q = 1, NY
            DO  1 P = 1, NX
               TEMP = P+(Q-1)*NX
               X = XX(TEMP)
               TEMP = P+(Q-1)*NX
               Y = YY(TEMP)
               A(P, Q, I, 1) = UX(P, Q, I)-.1*UY(P, Q, I)
               A(P, Q, I, 2) = UY(P, Q, I)+.1*UX(P, Q, I)
               AUX(P, Q, I, I, 1) = 1
               AUY(P, Q, I, I, 2) = 1
               AUY(P, Q, I, I, 1) = -.1
               AUX(P, Q, I, I, 2) = .1
               F(P, Q, 1) = UT(P, Q, 1)-X*Y
               FUT(P, Q, 1, 1) = 1
   1           CONTINUE
   2        CONTINUE
   3     CONTINUE
C MAP INTO (XI,ETA).
      CALL DTTGRG(NX, NY, D, NU, A, AU, AUX, AUY, F, FU, FUX, FUY)
      RETURN
      END
      SUBROUTINE BC(T, XI, NX, YI, NY, LX, RX, LY, RY, U, UT, UX
     1   , UY, UXT, UYT, NU, B, BU, BUT, BUX, BUY, BUXT, BUYT)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, XI(NX), YI(NY), LX, RX, LY
      DOUBLE PRECISION RY, U(NX, NY, NU), UT(NX, NY, NU), UX(NX, NY, NU)
     1   , UY(NX, NY, NU), UXT(NX, NY, NU)
      DOUBLE PRECISION UYT(NX, NY, NU), B(NX, NY, NU), BU(NX, NY, NU, 
     1   NU), BUT(NX, NY, NU, NU), BUX(NX, NY, NU, NU), BUY(NX, NY, NU
     2   , NU)
      DOUBLE PRECISION BUXT(NX, NY, NU, NU), BUYT(NX, NY, NU, NU)
      EXTERNAL BT, LR
      INTEGER I, J
      DOUBLE PRECISION D(600), X, Y, XX(100), YY(100)
      INTEGER TEMP1
      LOGICAL TEMP
      IF (NX*NY .GT. 100) CALL SETERR(19HBC - NX*NY .GT. 100, 19, 1, 2)
      CALL DBTMAP(T, XI, YI, NX, NY, LR, BT, XX, YY, D)
C MAP INTO (X,Y).
      CALL DTTGRU(NX, NY, D, UX, UY, UT, NU)
      DO  6 J = 1, NY
         DO  5 I = 1, NX
            TEMP1 = I+(J-1)*NX
            X = XX(TEMP1)
            TEMP1 = I+(J-1)*NX
            Y = YY(TEMP1)
            TEMP = XI(I) .EQ. LX
            IF (.NOT. TEMP) TEMP = XI(I) .EQ. RX
            IF (.NOT. TEMP) GOTO 1
               BU(I, J, 1, 1) = 1
C LEFT OR RIGHT.
               B(I, J, 1) = U(I, J, 1)-T*X*Y
               GOTO  4
   1           IF (YI(J) .NE. LY) GOTO 2
                  B(I, J, 1) = (UX(I, J, 1)-T*Y)-(UY(I, J, 1)-T*X)
C BOTTOM.
                  BUX(I, J, 1, 1) = 1
C NORMAL IS (1,-1).
                  BUY(I, J, 1, 1) = -1
                  GOTO  3
   2              B(I, J, 1) = (UY(I, J, 1)-T*X)-(UX(I, J, 1)-T*Y)
C TOP.
                  BUX(I, J, 1, 1) = -1
C NORMAL IS (-1,1).
                  BUY(I, J, 1, 1) = 1
   3        CONTINUE
   4        CONTINUE
   5        CONTINUE
   6     CONTINUE
C MAP INTO (XI,ETA).
      CALL DTTGRB(NX, NY, D, NU, BUX, BUY, BUT)
      RETURN
      END
      SUBROUTINE HANDLE(T0, U0, T, U, NV, DT, TSTOP)
      INTEGER NV
      DOUBLE PRECISION T0, U0(NV), T, U(NV), DT, TSTOP
      COMMON /D7TGRP/ ERRPAR, NU, MXQ, MYQ
      INTEGER NU, MXQ, MYQ
      REAL ERRPAR(2)
      COMMON /D7TGRM/ KX, IX, NX, KY, IY, NY
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
      DOUBLE PRECISION U(1), T
      COMMON /CSTAK/ DS
      DOUBLE PRECISION DS(500)
      INTEGER IFA, ITA(2), IXA(2), NTA(2), NXA(2), IDLUMD
      INTEGER IXS, IYS, NXS, NYS, ISTKGT, I
      INTEGER IEWE, KA(2), MA(2), IS(1000), I1MACH
      REAL RS(1000)
      LOGICAL LS(1000)
      COMPLEX CS(500)
      DOUBLE PRECISION DABS, ERRU, DMAX1, WS(500)
      INTEGER TEMP, TEMP1, TEMP2
      EQUIVALENCE (DS(1), CS(1), WS(1), RS(1), IS(1), LS(1))
C TO GET AND PRINT THE ERROR AT EACH TIME-STEP.
C U(NX-KX,NY-KY,NU).
C THE PORT LIBRARY STACK AND ITS ALIASES.
      CALL ENTER(1)
C FIND THE ERROR IN THE SOLUTION AT 2*KX * 2*KY POINTS / MESH RECTANGLE.
C X SEARCH GRID.
      IXS = IDLUMD(WS(IX), NX, 2*KX, NXS)
C Y SEARCH GRID.
      IYS = IDLUMD(WS(IY), NY, 2*KY, NYS)
C U SEARCH GRID VALUES.
      IEWE = ISTKGT(NXS*NYS, 4)
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
      IFA = ISTKGT(NXS*NYS, 4)
C EVALUATE THEM.
      CALL DTSD1(2, KA, WS, ITA, NTA, U, WS, IXA, NXA, MA, WS(IFA))
C ERROR IN SOLUTION VALUES.
      ERRU = 0
      TEMP = NXS*NYS
      DO  1 I = 1, TEMP
         TEMP2 = IEWE+I
         TEMP1 = IFA+I
         ERRU = DMAX1(ERRU, DABS(WS(TEMP2-1)-WS(TEMP1-1)))
   1     CONTINUE
      TEMP = I1MACH(2)
      WRITE (TEMP,  2) T, ERRU
   2  FORMAT (14H ERROR IN U(.,, 1PE10.2, 3H) =, 1PE10.2)
      CALL LEAVE
      RETURN
      END
      SUBROUTINE EWE(T, XI, NX, YI, NY, U, NU)
      INTEGER NU, NX, NY
      DOUBLE PRECISION T, XI(NX), YI(NY), U(NX, NY, NU)
      EXTERNAL BT, LR
      INTEGER I, J, P
      DOUBLE PRECISION D(6000), X, Y, XX(1000), YY(1000)
C THE EXACT SOLUTION.
      IF (NY .GT. 1000) CALL SETERR(18HEWE - NY .GT. 1000, 18, 1, 2)
      DO  3 P = 1, NU
         DO  2 I = 1, NX
            CALL DBTMAP(T, XI(I), YI, 1, NY, LR, BT, XX, YY, D)
            DO  1 J = 1, NY
               X = XX(J)
               Y = YY(J)
               U(I, J, P) = T*X*Y
   1           CONTINUE
   2        CONTINUE
   3     CONTINUE
      RETURN
      END
      SUBROUTINE LR(T, LX, RX, LXT, RXT)
      DOUBLE PRECISION T, LX, RX, LXT, RXT
C TO GET THE L AND R END-POINTS OF THE MAPPING IN X.
      LX = 0
      RX = 1
      LXT = 0
      RXT = 0
      RETURN
      END
      SUBROUTINE BT(T, X, F, G, FX, GX, FT, GT)
      DOUBLE PRECISION T, X, F, G, FX, GX
      DOUBLE PRECISION FT, GT
C TO GET THE BOTTOM AND TOP OF MAPPING IN Y.
      F = X-1D0
      G = X
      FT = 0
      GT = 0
      FX = 1
      GX = 1
      RETURN
      END
