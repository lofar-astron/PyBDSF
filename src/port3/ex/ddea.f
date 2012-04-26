C$TEST DDEA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE DDEA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM IODE
C
C***********************************************************************
      REAL TSTOP,V(2),DT
      REAL ERRPAR(2)
      INTEGER NV
      EXTERNAL DEE,HANDLE
C
      NV = 2
C
C  SET FOR 1E-2 ABSOLUTE ERROR.
C
      ERRPAR(1) = 0
      ERRPAR(2) = 1E-2
C
      TSTOP = 1E+20
      DT = 1E-7
C
C  INITIAL CONDITIONS FOR V.
C
      V(1) = 1
      V(2) = 1
C
      CALL IODE (V,NV,
     *           0E0,TSTOP,DT,
     *           DEE,
     *           ERRPAR,
     *           HANDLE)
C
      STOP
C
      END
      SUBROUTINE DEE(T,
     *               V,VT,NV,
     *               D,DV,DVT)
C
      REAL T,V(NV),VT(NV),D(NV),DV(NV,NV),DVT(NV,NV)
      INTEGER NV
C
      D(1) = VT(1)+2E0*VT(2) + V(1) + 2E+6*V(2)
      D(2) = 3E0*VT(1)+VT(2) + 3E0*V(1) + 1E+6*V(2)
C
      DVT(1,1) = 1
      DVT(1,2) = 2
      DV(1,1) = 1
      DV(1,2) = 2E+6
C
      DVT(2,1) = 3
      DVT(2,2) = 1
      DV(2,1) = 3
      DV(2,2) = 1E+6
C
      RETURN
C
      END
      SUBROUTINE HANDLE(T0,V0,T,V,NV,DT,TSTOP)
C
C OUTPUT AND CHECKING ROUTINE.
C
      REAL T0,V0(NV),T,V(NV),DT,TSTOP
      INTEGER NV
C
      REAL EV(2)
      INTEGER I1MACH
C
      IF ( T0 .EQ. T ) RETURN
C
      EV(1) = V(1) - EXP(-T)
      EV(2) = V(2) - EXP(-1E+6*T)
C
      IWUNIT = I1MACH(2)
      WRITE(IWUNIT,9000) T,EV(1),EV(2)
 9000 FORMAT(13H ERROR IN V( ,1P1E10.2,4H ) =,1P2E10.2)
C
      RETURN
C
      END
