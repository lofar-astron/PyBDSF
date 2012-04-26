C$TEST MNNA
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE MNNA
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM FMIN
C
C***********************************************************************
       EXTERNAL F
       INTEGER IWRITE,I1MACH
       REAL A,B,T,ANS,X
       IWRITE = I1MACH(2)
       A    =  .8
       B    = 1.2
       T    =  .0000001
       ANS  = FMIN(F,X,A,B,T)
       WRITE (IWRITE,9999) A,B,T
9999      FORMAT (5H A = ,1PE14.8,5H B = ,1PE14.8,5H T = ,1PE9.3)
       WRITE (IWRITE,9998) ANS
9998      FORMAT(16H THE MINIMUM IS ,1PE16.8)
       WRITE (IWRITE,9997) X
9997      FORMAT(14H IT OCCURS AT ,1PE18.8)
       STOP
       END
       FUNCTION F(X)
       F = -X * EXP(-X)
       RETURN
       END
