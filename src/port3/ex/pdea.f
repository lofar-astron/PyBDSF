C$TEST PDEA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE PDEA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM POST
C
C***********************************************************************
      REAL TSTOP,V( 1 ),DT,MESH( 100 ),U( 100 )
      REAL ERRPAR( 2 )
      INTEGER K,NMESH,NDX,NU,NV
      EXTERNAL AF,BC,DEE,HANDLE,UOFX
C
      COMMON/TIME/TT
      REAL TT
      COMMON/CSTAK/DS( 2000 )
      DOUBLEPRECISION DS
      REAL WS( 1000 )
      REAL RS( 1000 )
      INTEGER IS( 1000 )
      LOGICAL LS( 1000 )
      EQUIVALENCE( DS( 1 ),WS( 1 ),RS( 1 ),IS( 1 ),LS( 1 ) )
C
C  INITIALIZE THE PORT STACK LENGTH
C
      CALL ISTKIN( 2000,4 )
C
      NU = 1
      NV = 1
C
C  SET THE ERROR CRITERION FOR ABSOLUTE ERROR
C
      ERRPAR( 1 ) = 0
      ERRPAR( 2 ) = 1.E-2
C
      TSTOP = 8.*ATAN( 1.E0 )
      DT = 0.4
C
C MAKE A MESH OF NDX UNIFORM POINTS ON (-PI, +PI)
C
      K = 4
      NDX = 7
      CALL UMB(  - 4.*ATAN( 1.E0 ), + 4.*ATAN( 1.E0 ),NDX,K,MESH,NMESH )
      TT = 0
C
C  SET THE INITIAL CONDITIONS FOR U
C
      CALL L2SFF( UOFX,K,MESH,NMESH,U )
C
C  SET THE INITIAL CONDITIONS FOR V
C
      V( 1 ) =  - 1.
C
      CALL POST( U,NU,K,MESH,NMESH,V,NV,0E0,TSTOP,DT,AF,BC,DEE,ERRPAR,HA
     *NDLE )
C
      STOP
      END
      SUBROUTINE AF( T,X,NX,U,UX,UT,UTX,NU,V,VT,NV,A,AU,AUX,AUT,AUTX,AV,
     *AVT,F,FU,FUX,FUT,FUTX,FV,FVT )
      REAL T,X( NX ),U( NX,NU ),UX( NX,NU ),UT( NX,NU ),UTX( NX,NU ),V(
     *NV ),VT( NV ),A( NX,NU ),AU( NX,NU,NU ),AUX( NX,NU,NU ),AUT( NX,NU
     *,NU ),AUTX( NX,NU,NU ),AV( NX,NU,NV ),AVT( NX,NU,NV ),F( NX,NU ),F
     *U( NX,NU,NU ),FUX( NX,NU,NU ),FUT( NX,NU,NU ),FUTX( NX,NU,NU ),FV(
     * NX,NU,NV ),FVT( NX,NU,NV )
      INTEGER NU,NV,NX
      INTEGER I
      DO 23000 I = 1,NX
      A( I,1 ) =  - UX( I,1 )
      AUX( I,1,1 ) =  - 1
      F( I,1 ) =  - UT( I,1 ) - U( I,1 )**3 + SIN( X( I ) )*( COS( T ) -
     * SIN( T ) + SIN( X( I ) )**2*COS( T )**3 )
      FUT( I,1,1 ) =  - 1
      FU( I,1,1 ) =  - 3*U( I,1 )**2
23000 CONTINUE
      RETURN
      END
      SUBROUTINE BC( T,L,R,U,UX,UT,UTX,NU,V,VT,NV,B,BU,BUX,BUT,BUTX,BV,B
     *VT )
      REAL T,L,R,U( NU,2 ),UX( NU,2 ),UT( NU,2 ),UTX( NU,2 ),V( NV ),VT(
     * NV ),B( NU,2 ),BU( NU,NU,2 ),BUX( NU,NU,2 ),BUT( NU,NU,2 ),BUTX(
     *NU,NU,2 ),BV( NU,NV,2 ),BVT( NU,NV,2 )
      INTEGER NU,NV
      B( 1,1 ) = UX( 1,1 ) - V( 1 )
      B( 1,2 ) = UX( 1,2 ) - V( 1 )
      BUX( 1,1,1 ) = 1
      BV( 1,1,1 ) =  - 1
      BUX( 1,1,2 ) = 1
      BV( 1,1,2 ) =  - 1
      RETURN
      END
      SUBROUTINE DEE( T,K,X,NX,U,UT,NU,NXMK,V,VT,NV,D,DU,DUT,DV,DVT )
      REAL T,X( NX ),U( NXMK,NU ),UT( NXMK,NU ),V( NV ),VT( NV ),D( NV )
     *,DU( NV,NXMK,NU ),DUT( NV,NXMK,NU ),DV( NV,NV ),DVT( NV,NV )
      INTEGER K,NX,NU,NXMK,NV
      D( 1 ) = U( 1,1 ) - U( NX - K,1 )
      DU( 1,1,1 ) = 1
      DU( 1,NX - K,1 ) =  - 1
      RETURN
      END
      SUBROUTINE HANDLE( T0,U0,V0,T,U,V,NU,NXMK,NV,K,X,NX,DT,TSTOP )
      REAL T0,U0( NXMK,NU ),V0( NV ),T,U( NXMK,NU ),V( NV ),X( NX ),DT,T
     *STOP
      INTEGER NU,NXMK,NV,K,NX
      COMMON/TIME/TT
      REAL TT
      REAL EU,EESFF,EV
      INTEGER I1MACH
      EXTERNAL UOFX
      IF( T0 .EQ. T )GO TO 23002
      GO TO 23003
23002 CONTINUE
      RETURN
23003 CONTINUE
      TT = T
      EU = EESFF( K,X,NX,U,UOFX )
      EV = V( 1 ) + COS( T )
      IWUNIT = I1MACH( 2 )
      WRITE( IWUNIT,9001 )T,EU,EV
9001  FORMAT( 14H ERROR IN U(X,,1P1E10.2,4H ) =,1P1E10.2,6H   V =,1P4E10
     *.2 )
      RETURN
      END
      SUBROUTINE UOFX( X,NX,U,W )
      REAL X( NX ),U( NX ),W( NX )
      INTEGER NX
      COMMON/TIME/T
      REAL T
      INTEGER I
      DO 23005 I = 1,NX
      U( I ) = SIN( X( I ) )*COS( T )
23005 CONTINUE
      RETURN
      END
