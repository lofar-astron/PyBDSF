      SUBROUTINE I0TK01
C
      LOGICAL DONE
C
      DATA DONE /.FALSE./
C
      IF(DONE) RETURN
      DONE = .TRUE.
      IUNIT = I1MACH(4)
C
      WRITE( IUNIT, 100)
      WRITE( IUNIT, 200)
      WRITE( IUNIT, 300)
C
      RETURN
C
 100  FORMAT (1H1,
     *62H YOU HAVE USED, DIRECTLY OR INDIRECTLY, ONE OF THE STORAGE AL-/
     *62H LOCATION  PROGRAMS  IALLOC, DALLOC, STINIT, NIRALL, MTSTAK OR/
     *62H SRECAP.  THESE ARE BASED ON THE ASSUMPTION THAT ONE -UNIT- OF/
     *62H STORAGE  IS  ALLOCATED  TO  DATA OF TYPE LOGICAL, INTEGER AND/
     *62H REAL AND THAT TWO -UNITS- OF STORAGE ARE ALLOCATED TO DATA OF/
     *62H TYPE  DOUBLE PRECISION AND COMPLEX.  THIS ASSUMPTION PREVENTS/
     *62H MOVING PORT TO MANY MINI-COMPUTERS.                          /
     *62H                                                              /
     *62H TO OVERCOME THIS DIFFICULTY, THE PACKAGE HAS  BEEN  REWRITTEN/
     *62H WITH  NEW  NAMES AND SIMILAR CALLING SEQUENCES.  CALLS TO THE/
     *62H OLD SUBPROGRAMS SHOULD  BE  REPLACED  BY  CALLS  TO  THE  NEW/
     *62H PACKAGE  WHEN  CONVENIENT.   TO AVOID OBSOLETING OLD PROGRAMS/
     *62H THE OLD CALLING SEQUENCES WILL CONTINUE TO BE SUPPORTED.     /
     *62H                                                              /
     *)
C
 200  FORMAT(
     *62H THE OLD AND NEW CALLING SEQUENCES ARE AS FOLLOWS-            /
     *62H                                                              /
     *62H FUNCTION   OLD                       NEW                     /
     *62H                                                              /
     *62H GET        IX = IALLOC(NDATA,ISIZE)  IX = ISTKGT(NDATA,ITYPE)/
     *62H RELEASE    CALL DALLOC(NFRAMES)      CALL ISTKRL(NFRAMES)    /
     *62H INITIALIZE CALL STINIT(NDATA,ISIZE)  CALL ISTKIN(NDATA,ITYPE)/
     *62H MODIFY     IX = MTSTAK(NDATA)        IX = ISTKMD(NDATA)      /
     *62H STATISTICS CALL SRECAP(IUNIT)        - NO EQUIVALENT -       /
     *62H QUERY      N  = NIRALL(ISIZE)        N  = ISTKQU(ITYPE)      /
     *62H                                                              /
     *)
C
 300  FORMAT(
     *62H IN THE ABOVE ITYPE IS AS FOLLOWS-                            /
     *62H                                                              /
     *62H 1          LOGICAL                                           /
     *62H 2          INTEGER                                           /
     *62H 3          REAL                                              /
     *62H 4          DOUBLE PRECISION                                  /
     *62H 5          COMPLEX                                           /
     *62H                                                              /
     *62H NOTE ALSO  THAT ALLOCATIONS SHOULD  NOT  BE SPLIT INTO SUBAL-/
     *62H LOCATIONS  OF  DIFFERENT  TYPE  AS THIS ALSO COMPROMISES POR-/
     *62H TABILITY.                                                    /
     *)
C
      END
