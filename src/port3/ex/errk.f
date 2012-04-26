C$TEST ERRK
C TO RUN AS A MAIN PROGRAM REMOVE NEXT LINE
      SUBROUTINE ERRK
C***********************************************************************
C
C  EXAMPLE OF USE OF THE PORT PROGRAM STKDMP
C
C***********************************************************************
C     SAMPLE USE OF THE STACK DUMP
      INTEGER IPTR, ISTKGT
C
      COMMON  /CSTAK/   DSTAK
      DOUBLE PRECISION  DSTAK(500)
      INTEGER           ISTAK(1000)
      LOGICAL           LSTAK(1000)
      REAL              RSTAK(1000)
      COMPLEX           CMSTAK(500)
      EQUIVALENCE (DSTAK(1), ISTAK(1))
      EQUIVALENCE (DSTAK(1), LSTAK(1))
      EQUIVALENCE (DSTAK(1), RSTAK(1))
      EQUIVALENCE (DSTAK(1), CMSTAK(1))
C
      IPTR = ISTKGT(25, 1)
      CALL SETL(25, .FALSE., LSTAK(IPTR))
      IPTR = ISTKGT(25, 2)
      CALL SETI(25, -1, ISTAK(IPTR))
      IPTR = ISTKGT(25, 3)
      CALL SETR(25, 1.0, RSTAK(IPTR))
      IPTR = ISTKGT(25, 4)
      CALL SETD(25, 1.0D0, DSTAK(IPTR))
      IPTR = ISTKGT(25, 5)
      CALL SETC(25, CMPLX(1.0, -1.0), CMSTAK(IPTR))
      IPTR = ISTKGT(25, 5)
      CALL SETC(25, CMPLX(1.0, -1.0), CMSTAK(IPTR))
      CALL ISTKRL(1)
C
      CALL STKDMP
      STOP
      END
