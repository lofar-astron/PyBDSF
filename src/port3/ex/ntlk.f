C$TEST NTLK
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE NTLK
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM SMNH
C
C***********************************************************************
        INTEGER N
        EXTERNAL ROSN,ROSGH
        REAL X(2)
        N=2
C INITIALIZE X
        X(1)=-1.2
        X(2)=1.0
C
C SOLVE THE PROBLEM
C
        CALL SMNH(N, X, ROSN, ROSGH, 100, 1.E-4)
C       PRINT RESULTS ON STANDARD OUTPUT UNIT
        IWRITE=I1MACH(2)
        WRITE(IWRITE,10)(X(I),I=1,N)
 10     FORMAT(10H SOLUTION-,5E15.5)
        STOP
        END
        SUBROUTINE ROSN(N,X,NF,F)
C THIS SUBROUTINE COMPUTES THE  FUNCTION
        INTEGER N, NF
        REAL X(N), F
        F=100.0*(X(2)-X(1)*X(1))**2 + (1.0 - X(1))**2
        RETURN
        END
        SUBROUTINE ROSGH(N,X,NF,G,H)
C THIS SUBROUTINE COMPUTES THE GRADIENT AND THE HESSIAN
        INTEGER N,NF
        REAL X(N), G(N), H(1)
        G(1)=200.0*(X(2)-X(1)*X(1))*(-2.0)*X(1) - 2.0*(1-X(1))
        G(2)=200.0*(X(2)-X(1)*X(1))
C H(1) HAS THE (1,1) ELEMENT, H(2) HAS THE (2,1) ELEMENT,
C H(3) HAS THE (2,2) ELEMENT OF THE MATRIX OF SECOND PARTIALS
        H(1)=200.0*(X(2)-X(1)*X(1))*(-2.0)+800.0*X(1)*X(1)+2.0
        H(2)=-400.0*X(1)
        H(3)=200.0
        RETURN
        END
