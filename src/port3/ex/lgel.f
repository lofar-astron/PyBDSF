C$TEST LGEL
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE LGEL
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM GEBS
C
C***********************************************************************
        INTEGER N, I, J, IWRITE, I1MACH
        REAL A(15,15), B(15)
        N=15
C
C FORM THE MATRIX AND SET THE RIGHT-HAND SIDE
C TO THE LAST COLUMN OF THE IDENTITY MATRIX
        DO 20 I=1,N
           DO 10 J=I,N
              A(I,J) = -1.0
  10       CONTINUE
           A(I,I) = 1.0
           B(I)   = 0.0
  20    CONTINUE
        B(N)=1.0
C FIND THE LAST COLUMN OF THE INVERSE MATRIX
        CALL GEBS(N,A,15,B,N,1)
        IWRITE=I1MACH(2)
        WRITE(IWRITE,21)(I,B(I),I=1,N)
  21    FORMAT(3H B(,I3,3H )=,F15.4)
        STOP
        END
