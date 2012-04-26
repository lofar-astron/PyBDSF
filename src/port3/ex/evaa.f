C$TEST EVAA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE EVAA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM EIGEN
C
C***********************************************************************
      REAL A(4,4),ORT(4),Z(4,4)
      REAL H(4,4),WR(4),WI(4)
C
      DATA A(1,1),A(1,2),A(1,3),A(1,4) / 3., 1., 2., 5. /
      DATA A(2,1),A(2,2),A(2,3),A(2,4) / 2., 1., 3., 7. /
      DATA A(3,1),A(3,2),A(3,3),A(3,4) / 3., 1., 2., 4. /
      DATA A(4,1),A(4,2),A(4,3),A(4,4) / 4., 1., 3., 2. /
C
      NM=4
      N=4
C
C     SET OUTPUT WRITE UNIT
C
       IWUNIT=I1MACH(2)
C
      CALL EIGEN(NM,N,A,WR,WI,Z)
C
      WRITE (IWUNIT,96)
  96  FORMAT (22H0THE EIGENVALUES ARE -/)
C
      WRITE (IWUNIT,97) (WR(J),WI(J),J=1,N)
  97  FORMAT (/1X,2E20.8)
C
      DO 20 K=1,N
      SCALE=AMAX1(ABS(Z(1,K)),ABS(Z(2,K)),ABS(Z(3,K)),ABS(Z(4,K)))
      DO 20 J=1,N
 20   Z(J,K)=Z(J,K)/SCALE
C
      WRITE (IWUNIT,98)
  98  FORMAT (30H0THE SCALED EIGENVECTORS ARE -//)
C
      WRITE (IWUNIT,99) ((Z(J,K),K=1,N),J=1,N)
  99  FORMAT (1X,1P4E18.8/)
C
      STOP
      END
