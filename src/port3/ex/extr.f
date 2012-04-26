C$TEST EXTR
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE EXTR
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM EXTRMX
C
C***********************************************************************
       INTEGER  IWRITE,IEXT(100),NEX,IMAX,IMIN,IMAG
       INTEGER  I1MACH,I,J
       REAL     PI,STEP,X,F(100)
C
       IWRITE = I1MACH(2)
       PI = 3.1415926532
       STEP = 2.0*PI/99.0
       DO 10 I=1,100
          X = STEP*FLOAT(I-1)
  10      F(I) = EXP(-X)*COS(X)
C
       CALL EXTRMR(100,F,NEX,IEXT,IMAX,IMIN,IMAG)
C
       WRITE(IWRITE,20)
  20   FORMAT(6X,9HEXTREMALS/5X,1HX,10X,4HF(X))
       DO 30 J=1,NEX
          I =IEXT(J)
          X = STEP*FLOAT(I-1)
  30      WRITE(IWRITE,40) X,F(I)
  40   FORMAT(2F10.5)
       STOP
       END
