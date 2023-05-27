!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! PROGRAM: SH23
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Author: Tianhong Huang
! Date:   MAY 2023
! Version: hb_1km_v5.1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PROGRAM SH23

    IMPLICIT NONE

    INCLUDE 'integer_types.h'
    INCLUDE 'real_types.h'

    REAL(dp) :: startTime
    REAL(dp) :: endTime

    print *, NEW_LINE(' ')
    print *, NEW_LINE(' ')

    CALL CPU_TIME(startTime)
    CALL MAIN
    CALL CPU_TIME(endTime)

    print '(a)', repeat('-', 60)
    print *, NEW_LINE(' '), 'Execution time: ', endTime - startTime
    print '(a)', repeat('-', 60)
    print *, NEW_LINE(' ')
    print *, NEW_LINE(' ')

CONTAINS

SUBROUTINE MAIN

    USE MPI
    USE INITIALIZE
    USE ITERATIONS
    USE WRITE_OUTPUTS

    IMPLICIT NONE

    INTEGER(qb) :: ierror

    CALL MPI_INIT(ierror)

    CALL INIT_PARAM

    CALL INIT_GRIDS

    CALL INIT_OPDIR

    CALL INIT_EIGEN

    CALL INIT_RANDN

    CALL INIT_DEBUG

    CALL iterate_steps

    CALL MPI_FINALIZE(ierror)

END SUBROUTINE MAIN

END PROGRAM SH23
