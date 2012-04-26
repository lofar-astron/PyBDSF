C$TEST LLZA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LLZA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM LZ
C
C***********************************************************************
       COMPLEX B(5,5),A(5,5),EIGA(5),EIGB(5),X,EIG
       IIN=I1MACH(1)
       IOUT=I1MACH(2)
C
C READ IN MATRICES
C
       READ(IIN,10)((A(I,J),J=1,5),I=1,5)
       READ(IIN,10)((B(I,J),J=1,5),I=1,5)
  10   FORMAT(10F6.0)
C
C PRINT MATRICES
C
       WRITE(IOUT,20)
  20   FORMAT(13H THE A MATRIX)
       WRITE(IOUT,30)((A(I,J),J=1,5),I=1,5)
  30   FORMAT(5(F6.0,2H+ ,F6.0,1HI))
       WRITE(IOUT,40)
  40   FORMAT(13H THE B MATRIX)
       WRITE(IOUT,30)((B(I,J),J=1,5), I=1,5)
C
C SOLVE THE EIGENVALUE PROBLEM
C
       CALL LZ(5,A,5,B,5,X,1,.FALSE.,EIGA,EIGB)
       WRITE(IOUT,50)
  50   FORMAT(10X,4HEIGA,16X,4HEIGB,22X,10HEIGENVALUE)
       DO 60 I=1,5
          EIG=CMPLX(R1MACH(2),R1MACH(2))
          IF(REAL(EIGB(I)).NE.0.0.OR.AIMAG(EIGB(I)).NE.0.0)
     1    EIG=EIGA(I)/EIGB(I)
          WRITE(IOUT,70)EIGA(I),EIGB(I),EIG
  60   CONTINUE
  70   FORMAT(1H ,2E10.3,2X,2E10.3,2X,2E16.8)
       STOP
       END
C
C DATA FOR THE EXAMPLE IN THE PORT SHEET...  (REMOVE THE C
C IN COLUMN 1 BEFORE FEEDING THIS DATA TO THE PROGRAM ABOVE.)
C$DATA
C   41. -369. -143. -747.  -20.-1368.   20.  486.  104. -432.
C  148.  261.  144.  666.   -6.-1152.  -78.   45.    8. -540.
C  -19.  819.   87.  243.    4. 1548.  -56. -954. -164.  180.
C  -60. -945.  -81. -279.   99.  171.   34.  441.   84. -144.
C    1. -468.  133.  747.  132.  774.  -46.  -45.  -12. -216.
C   90.  161.  180.  335.   36.  182.  -90. -162.  -72.  -36.
C -105. -169. -210. -322.  -42.   24.  105.  167.   84.  204.
C  -90. -211. -180. -307.  -36. -160.   90.  186.   72.   36.
C   75.  205.  150.  215.   30.   45.  -75. -165.  -60.  -80.
C  -75.  -48. -150. -299.  -30. -102.   75.   89.   60.   88.
