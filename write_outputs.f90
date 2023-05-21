!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MODULE: write_outputs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Version: hb_1km_v5.1
! Update:  new qf eqtn solving by divergence/advection form
! Update:  setting qf as qfhat constant in fluid core
! Update:  adjusting threshold of zero Fmodes to 1E-14
! Update:  adding histogram statistics
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Program: SH23
! Author:  Tianhong Huang
! Date:    MAY 2023
! SBRs:    output_model_varbs
!          output_statistics
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE write_outputs

IMPLICIT NONE

PRIVATE

INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: output_model_varbs
PUBLIC :: output_statistics

CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: output_model_varbs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[in] step_num: iteration steps number.
! @param[in] op_idx:   index noting model variables op status.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE output_model_varbs(step_num)

    USE MPI
    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: step_num
    INTEGER(qb) :: i, j
    CHARACTER(LEN=66) :: filename_1
    CHARACTER(LEN=66) :: filename_2
    CHARACTER(LEN=66) :: filename_3
    CHARACTER(LEN=66) :: filename_4
    CHARACTER(LEN=66) :: filename_5
    CHARACTER(LEN=66) :: filename_6
    CHARACTER(LEN=66) :: filename_7
    CHARACTER(LEN=66) :: filename_8
    CHARACTER(LEN=66) :: filename_9
    CHARACTER(LEN=66) :: filename_10
    CHARACTER(LEN=66) :: filename_11
    CHARACTER(LEN=66) :: filename_12
    CHARACTER(LEN=66) :: filename_13
    CHARACTER(LEN=66) :: filename_14
    CHARACTER(LEN=66) :: filename_15
    CHARACTER(LEN=66) :: filename_101

    !< save ub_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_1,'(A,I0.4,A)') './output_serial/ub/',step_num,'.csv'
    ELSE
        WRITE(filename_1,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/ub/', step_num, '.csv'
    END IF
    OPEN(101,file=filename_1,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(101,'(E32.16,A,1x)',ADVANCE='NO') ub(i,j,op_idx), ','
    END DO
    WRITE(101,'(1x)')
    END DO
    CLOSE(101)

    !< save vb_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_2,'(A,I0.4,A)') './output_serial/vb/',step_num,'.csv'
    ELSE
        WRITE(filename_2,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/vb/', step_num, '.csv'
    END IF
    OPEN(102,file=filename_2,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(102,'(E32.16,A,1x)',ADVANCE='NO') vb(i,j,op_idx), ','
    END DO
    WRITE(102,'(1x)')
    END DO
    CLOSE(102)

    !< save u0_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_3,'(A,I0.4,A)') './output_serial/u0/',step_num,'.csv'
    ELSE
        WRITE(filename_3,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/u0/', step_num, '.csv'
    END IF
    OPEN(103,file=filename_3,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(103,'(E32.16,A,1x)',ADVANCE='NO') u0(i,j,op_idx), ','
    END DO
    WRITE(103,'(1x)')
    END DO
    CLOSE(103)

    !< save v0_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_4,'(A,I0.4,A)') './output_serial/v0/',step_num,'.csv'
    ELSE
        WRITE(filename_4,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/v0/', step_num, '.csv'
    END IF
    OPEN(104,file=filename_4,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(104,'(E32.16,A,1x)',ADVANCE='NO') v0(i,j,op_idx), ','
    END DO
    WRITE(104,'(1x)')
    END DO
    CLOSE(104)

    !< save u1_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_5,'(A,I0.4,A)') './output_serial/u1/',step_num,'.csv'
    ELSE
        WRITE(filename_5,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/u1/', step_num, '.csv'
    END IF
    OPEN(105,file=filename_5,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
       WRITE(105,'(E32.16,A,1x)',ADVANCE='NO') u1(i,j,op_idx), ','
    END DO
    WRITE(105,'(1x)')
    END DO
    CLOSE(105)

    !< save v1_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_6,'(A,I0.4,A)') './output_serial/v1/',step_num,'.csv'
    ELSE
        WRITE(filename_6,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/v1/', step_num, '.csv'
    END IF
    OPEN(106,file=filename_6,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
       WRITE(106,'(E32.16,A,1x)',ADVANCE='NO') v1(i,j,op_idx), ','
    END DO
    WRITE(106,'(1x)')
    END DO
    CLOSE(106)

    !< save qtb_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_7,'(A,I0.4,A)') './output_serial/qtb/',step_num,'.csv'
    ELSE
        WRITE(filename_7,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/qtb/', step_num, '.csv'
    END IF
    OPEN(107,file=filename_7,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
       WRITE(107,'(E32.16,A,1x)',ADVANCE='NO') qtb(i,j,op_idx), ','
    END DO
    WRITE(107,'(1x)')
    END DO
    CLOSE(107)

    !< save qf_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_8,'(A,I0.4,A)') './output_serial/qf/',step_num,'.csv'
    ELSE
        WRITE(filename_8,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/qf/', step_num, '.csv'
    END IF
    OPEN(108,file=filename_8,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(108,'(E32.16,A,1x)',ADVANCE='NO') qf(i,j,op_idx), ','
    END DO
    WRITE(108,'(1x)')
    END DO
    CLOSE(108)

    !< save Scloud_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_9,'(A,I0.4,A)') './output_serial/shallow-cloud/',step_num,'.csv'
    ELSE
        WRITE(filename_9,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/shallow-cloud/', step_num, '.csv'
    END IF
    OPEN(109,file=filename_9,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(109,'(E32.16,A,1x)',ADVANCE='NO') sig_c(i,j), ','
    END DO
    WRITE(109,'(1x)')
    END DO
    CLOSE(109)

    !< save To_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_10,'(A,I0.4,A)') './output_serial/To/',step_num,'.csv'
    ELSE
        WRITE(filename_10,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/To/', step_num, '.csv'
    END IF
    OPEN(110,file=filename_10,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(110,'(E32.16,A,1x)',ADVANCE='NO') To(i,j,op_idx), ','
    END DO
    WRITE(110,'(1x)')
    END DO
    CLOSE(110)

    !< save Tb_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_11,'(A,I0.4,A)') './output_serial/Tb/',step_num,'.csv'
    ELSE
        WRITE(filename_11,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/Tb/', step_num, '.csv'
    END IF
    OPEN(111,file=filename_11,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(111,'(E32.16,A,1x)',ADVANCE='NO') Tb(i,j), ','
    END DO
    WRITE(111,'(1x)')
    END DO
    CLOSE(111)

    !< save thetaeb_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_12,'(A,I0.4,A)') './output_serial/thetaeb/',step_num,'.csv'
    ELSE
        WRITE(filename_12,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/thetaeb/', step_num, '.csv'
    END IF
    OPEN(112,file=filename_12,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(112,'(E32.16,A,1x)',ADVANCE='NO') thetaeb(i,j,op_idx), ','
    END DO
    WRITE(112,'(1x)')
    END DO
    CLOSE(112)


    !< save Tf_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_13,'(A,I0.4,A)') './output_serial/Tf/',step_num,'.csv'
    ELSE
        WRITE(filename_13,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/Tf/', step_num, '.csv'
    END IF
    OPEN(113,file=filename_13,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(113,'(E32.16,A,1x)',ADVANCE='NO') Tf(i,j), ','
    END DO
    WRITE(113,'(1x)')
    END DO
    CLOSE(113)

    !< save theta1_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_14,'(A,I0.4,A)') './output_serial/theta1/',step_num,'.csv'
    ELSE
        WRITE(filename_14,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/theta1/', step_num, '.csv'
    END IF
    OPEN(114,file=filename_14,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(114,'(E32.16,A,1x)',ADVANCE='NO') theta1(i,j,op_idx), ','
    END DO
    WRITE(114,'(1x)')
    END DO
    CLOSE(114)

    !< save Dcloud_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_15,'(A,I0.4,A)') './output_serial/deep-cloud/',step_num,'.csv'
    ELSE
        WRITE(filename_15,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/deep-cloud/', step_num, '.csv'
    END IF
    OPEN(115,file=filename_15,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(115,'(E32.16,A,1x)',ADVANCE='NO') sig_f(i,j), ','
    END DO
    WRITE(115,'(1x)')
    END DO
    CLOSE(115)

    !< save Bcloud_step_num.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_101,'(A,I0.4,A)') './output_serial/shallow-cloud-ast/',step_num,'.csv'
    ELSE
        WRITE(filename_101,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/shallow-cloud-ast/', step_num, '.csv'
    END IF
    OPEN(201,file=filename_101,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
    WRITE(201,'(E32.16,A,1x)',ADVANCE='NO') sig_b(i,j), ','
    END DO
    WRITE(201,'(1x)')
    END DO
    CLOSE(201)


END SUBROUTINE output_model_varbs


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: output_statistics
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE output_statistics

    USE MPI
    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb) :: i, j, t
    CHARACTER(LEN=66) :: filename_16
    CHARACTER(LEN=66) :: filename_17
    CHARACTER(LEN=66) :: filename_18
    CHARACTER(LEN=66) :: filename_19
    CHARACTER(LEN=66) :: filename_20
    CHARACTER(LEN=66) :: filename_21
    CHARACTER(LEN=66) :: filename_22
    CHARACTER(LEN=66) :: filename_23
    CHARACTER(LEN=66) :: filename_24
    CHARACTER(LEN=66) :: filename_26
    CHARACTER(LEN=66) :: filename_27
    CHARACTER(LEN=66) :: filename_28
    CHARACTER(LEN=66) :: filename_29
    CHARACTER(LEN=66) :: filename_30
    CHARACTER(LEN=66) :: filename_31
    CHARACTER(LEN=66) :: filename_32
    CHARACTER(LEN=66) :: filename_35
    CHARACTER(LEN=66) :: filename_36
    CHARACTER(LEN=66) :: filename_37
    CHARACTER(LEN=66) :: filename_38
    CHARACTER(LEN=66) :: filename_39

    !< save shallow_cloud_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_16,'(A,I0.4,A)') './output_serial/stats/shallow_cloud_tavg.csv'
    ELSE
        WRITE(filename_16,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/shallow_cloud_tavg.csv'
    END IF
    OPEN(116,file=filename_16,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(116,'(E32.16,A,1x)',ADVANCE='NO') Scloud_tavg(i,j), ','
    END DO
    WRITE(116,'(1x)')
    END DO
    CLOSE(116)

    !< save deep_cloud_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_17,'(A,I0.4,A)') './output_serial/stats/deep_cloud_tavg.csv'
    ELSE
        WRITE(filename_17,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/deep_cloud_tavg.csv'
    END IF
    OPEN(117,file=filename_17,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(117,'(E32.16,A,1x)',ADVANCE='NO') Dcloud_tavg(i,j), ','
    END DO
    WRITE(117,'(1x)')
    END DO
    CLOSE(117)

    !< save shallow_cloud_ast_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_38,'(A,I0.4,A)') './output_serial/stats/shallow_cloud_ast_tavg.csv'
    ELSE
        WRITE(filename_38,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/shallow_cloud_ast_tavg.csv'
    END IF
    OPEN(138,file=filename_38,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(138,'(E32.16,A,1x)',ADVANCE='NO') Bcloud_tavg(i,j), ','
    END DO
    WRITE(138,'(1x)')
    END DO
    CLOSE(138)

    !< save shallow_cloud_rate.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_18,'(A,I0.4,A)') './output_serial/stats/shallow_cloud_rate.csv'
    ELSE
        WRITE(filename_18,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/shallow_cloud_rate.csv'
    END IF
    OPEN(118,file=filename_18,form='formatted')
    DO t = 1, numSteps
       WRITE(118,'(E32.16,A,1x)',ADVANCE='NO') Scloud_rate(t), ','
       WRITE(118,'(1x)')
    END DO
    CLOSE(118)

    !< save deep_cloud_rate.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_19,'(A,I0.4,A)') './output_serial/stats/deep_cloud_rate.csv'
    ELSE
        WRITE(filename_19,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/deep_cloud_rate.csv'
    END IF
    OPEN(119,file=filename_19,form='formatted')
    DO t = 1, numSteps
       WRITE(119,'(E32.16,A,1x)',ADVANCE='NO') Dcloud_rate(t), ','
       WRITE(119,'(1x)')
    END DO
    CLOSE(119)

    !< save shallow_cloud_ast_rate.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_39,'(A,I0.4,A)') './output_serial/stats/shallow_cloud_ast_rate.csv'
    ELSE
        WRITE(filename_39,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/shallow_cloud_ast_rate.csv'
    END IF
    OPEN(139,file=filename_39,form='formatted')
    DO t = 1, numSteps
       WRITE(139,'(E32.16,A,1x)',ADVANCE='NO') Bcloud_rate(t), ','
       WRITE(139,'(1x)')
    END DO
    CLOSE(139)

    !< save u1_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_20,'(A,I0.4,A)') './output_serial/stats/u1_tavg.csv'
    ELSE
        WRITE(filename_20,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/u1_tavg.csv'
    END IF
    OPEN(120,file=filename_20,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(120,'(E32.16,A,1x)',ADVANCE='NO') u1_tavg(i,j), ','
    END DO
    WRITE(120,'(1x)')
    END DO
    CLOSE(120)

    !< save v1_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_21,'(A,I0.4,A)') './output_serial/stats/v1_tavg.csv'
    ELSE
        WRITE(filename_21,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/v1_tavg.csv'
    END IF
    OPEN(121,file=filename_21,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(121,'(E32.16,A,1x)',ADVANCE='NO') v1_tavg(i,j), ','
    END DO
    WRITE(121,'(1x)')
    END DO
    CLOSE(121)

    !< save qtb_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_22,'(A,I0.4,A)') './output_serial/stats/qtb_tavg.csv'
    ELSE
        WRITE(filename_22,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/qtb_tavg.csv'
    END IF
    OPEN(122,file=filename_22,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(122,'(E32.16,A,1x)',ADVANCE='NO') qtb_tavg(i,j), ','
    END DO
    WRITE(122,'(1x)')
    END DO
    CLOSE(122)

    !< save qf_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_23,'(A,I0.4,A)') './output_serial/stats/qf_tavg.csv'
    ELSE
        WRITE(filename_23,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/qf_tavg.csv'
    END IF
    OPEN(123,file=filename_23,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(123,'(E32.16,A,1x)',ADVANCE='NO') qf_tavg(i,j), ','
    END DO
    WRITE(123,'(1x)')
    END DO
    CLOSE(123)

    !< save Tfrad_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_24,'(A,I0.4,A)') './output_serial/stats/Tfrad_tavg.csv'
    ELSE
        WRITE(filename_24,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Tfrad_tavg.csv'
    END IF
    OPEN(124,file=filename_24,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(124,'(E32.16,A,1x)',ADVANCE='NO') Tfrad_tavg(i,j), ','
    END DO
    WRITE(124,'(1x)')
    END DO
    CLOSE(124)

    !< save Tfconv_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_26,'(A,I0.4,A)') './output_serial/stats/Tfconv_tavg.csv'
    ELSE
        WRITE(filename_26,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Tfconv_tavg.csv'
    END IF
    OPEN(126,file=filename_26,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(126,'(E32.16,A,1x)',ADVANCE='NO') Tfconv_tavg(i,j), ','
    END DO
    WRITE(126,'(1x)')
    END DO
    CLOSE(126)

    !< save Torad_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_27,'(A,I0.4,A)') './output_serial/stats/Torad_tavg.csv'
    ELSE
        WRITE(filename_27,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Torad_tavg.csv'
    END IF
    OPEN(127,file=filename_27,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(127,'(E32.16,A,1x)',ADVANCE='NO') Torad_tavg(i,j), ','
    END DO
    WRITE(127,'(1x)')
    END DO
    CLOSE(127)

    !< save Tbrad_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_28,'(A,I0.4,A)') './output_serial/stats/Tbrad_tavg.csv'
    ELSE
        WRITE(filename_28,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Tbrad_tavg.csv'
    END IF
    OPEN(128,file=filename_28,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(128,'(E32.16,A,1x)',ADVANCE='NO') Tbrad_tavg(i,j), ','
    END DO
    WRITE(128,'(1x)')
    END DO
    CLOSE(128)

    !< save To_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_35,'(A,I0.4,A)') './output_serial/stats/To_tavg.csv'
    ELSE
        WRITE(filename_35,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/To_tavg.csv'
    END IF
    OPEN(135,file=filename_35,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(135,'(E32.16,A,1x)',ADVANCE='NO') To_tavg(i,j), ','
    END DO
    WRITE(135,'(1x)')
    END DO
    CLOSE(135)

    !< save Tb_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_36,'(A,I0.4,A)') './output_serial/stats/Tb_tavg.csv'
    ELSE
        WRITE(filename_36,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Tb_tavg.csv'
    END IF
    OPEN(136,file=filename_36,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(136,'(E32.16,A,1x)',ADVANCE='NO') Tb_tavg(i,j), ','
    END DO
    WRITE(136,'(1x)')
    END DO
    CLOSE(136)

    !< save Tf_tavg.csv
    IF (sch%commSize.EQ.1) THEN
        WRITE(filename_37,'(A,I0.4,A)') './output_serial/stats/Tf_tavg.csv'
    ELSE
        WRITE(filename_37,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/Tf_tavg.csv'
    END IF
    OPEN(137,file=filename_37,form='formatted')
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        WRITE(137,'(E32.16,A,1x)',ADVANCE='NO') Tf_tavg(i,j), ','
    END DO
    WRITE(137,'(1x)')
    END DO
    CLOSE(137)

    !< save RainSize_hist.csv
    IF (sch%commSize .EQ. 1) THEN
        WRITE(filename_29,'(A,I0.4,A)') './output_serial/stats/RainSize_hist.csv'
    ELSE
        WRITE(filename_29,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/RainSize_hist.csv'
    END IF
    OPEN(129,file=filename_29,form='formatted')
    DO i = 0, RSBinNum-1
       WRITE(129,'(E32.16,A,1x)',ADVANCE='NO') RainSizeHist(i), ','
       WRITE(129,'(1x)')
    END DO
    CLOSE(129)

    !< save RainTime_hist.csv
    IF (sch%commSize .EQ. 1) THEN
        WRITE(filename_30,'(A,I0.4,A)') './output_serial/stats/RainTime_hist.csv'
    ELSE
        WRITE(filename_30,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/RainTime_hist.csv'
    END IF
    OPEN(130,file=filename_30,form='formatted')
    DO i = 0, RTBinNum-1
       WRITE(130,'(E32.16,A,1x)',ADVANCE='NO') RainTimeHist(i), ','
       WRITE(130,'(1x)')
    END DO
    CLOSE(130)

    !< save DryTime_hist.csv
    IF (sch%commSize .EQ. 1) THEN
        WRITE(filename_31,'(A,I0.4,A)') './output_serial/stats/DryTime_hist.csv'
    ELSE
        WRITE(filename_31,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/DryTime_hist.csv'
    END IF
    OPEN(131,file=filename_31,form='formatted')
    DO i = 0, DTBinNum-1
       WRITE(131,'(E32.16,A,1x)',ADVANCE='NO') DryTimeHist(i), ','
       WRITE(131,'(1x)')
    END DO
    CLOSE(131)

    !< save SigcTime_hist.csv
    IF (sch%commSize .EQ. 1) THEN
        WRITE(filename_32,'(A,I0.4,A)') './output_serial/stats/SigcTime_hist.csv'
    ELSE
        WRITE(filename_32,'(A,I0.2,A,I0.2,A,I0.4,A)') './output_parallel/MAX_', &
        sch%commSize, '_ID_', sch%procID, '/stats/SigcTime_hist.csv'
    END IF
    OPEN(132,file=filename_32,form='formatted')
    DO i = 0, STBinNum-1
       WRITE(132,'(E32.16,A,1x)',ADVANCE='NO') SigcTimeHist(i), ','
       WRITE(132,'(1x)')
    END DO
    CLOSE(132)

END SUBROUTINE output_statistics

END MODULE write_outputs


