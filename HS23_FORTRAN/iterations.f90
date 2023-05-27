!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! MOUDLE: ITERATIONS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Version: hb_1km_v5.1
! Update:  new qf eqtn solving by divergence/advection form
! Update:  setting qf as qfhat constant in fluid core
! Update:  adjusting threshold of zero Fmodes to 1E-14
! Update:  adding histogrma data statistics
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Program: SH23
! Author:  Tianhong Huang
! Date:    MAY 2023
! Intro:   module handle SH23 model simulations, involving 4 schemes.
!          physical ODE, fluid dynamics, stochastic equation and cloud indicator.
! SBRs:    iterate_steps
!          scheme_physc
!          scheme_fluid
!          scheme_stoch
!          scheme_randn
!          scheme_cloud
!          update_varbs
!          lnsys_solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE ITERATIONS

IMPLICIT NONE

PRIVATE

INCLUDE 'integer_types.h'
INCLUDE 'real_types.h'

PUBLIC :: iterate_steps

CONTAINS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: iterate_steps
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE iterate_steps

    USE MPI
    USE INITIALIZE
    USE WRITE_OUTPUTS
    USE PARALLEL_TOOLS
    USE PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER(qb) :: i
    INTEGER(qb) :: step
    INTEGER(qb) :: hrs_idx
    INTEGER(qb) :: step_num
    INTEGER(qb) :: time_idx
    REAL(dp) :: simulation_T
    REAL(dp) :: thresholds_T
    CHARACTER(LEN=66) :: time_sequence_name

    hrs_idx  = 0_qb
    step_num = 0_qb
    time_idx = 1_qb

    !< save initial condition
    CALL output_model_varbs(step_num)

    DO step = 1, numSteps

        !< physical ODE scheme
        IF (flag_physc .EQ. 1) THEN
            CALL scheme_physc(step)
        END IF

        !< fluid dynamics scheme
        IF (flag_fluid .EQ. 1) THEN
            CALL scheme_fluid
        END IF

        !< stochastic equation scheme
        IF (flag_stoch .EQ. 1) THEN
            CALL scheme_randn
        END IF

        !< update cloud status
        IF (flag_physc .EQ. 1) THEN
            CALL scheme_cloud
        END IF

        !< update model variables
        CALL update_varbs(step)

        !< output model variables at specific time steps
        IF (MOD(step, op_Freq) .EQ. 0) THEN
            step_num = step_num + 1_qb
            time_idx = time_idx + 1_qb
            op_time(time_idx) = step*dt/3600.0_dp/24.0_dp
            CALL output_model_varbs(step_num)
        END IF

        !< display simulation time info
        simulation_T = dt * (step+1_qb)
        thresholds_T = hrs_idx * 3600.0_dp * 24.0_dp
        IF ((simulation_T.GT.thresholds_T) .AND. (sch%procID.EQ.0)) THEN
            print '(a)', repeat('-', 60)
            print *, NEW_LINE(' '), 'model simulation time on No.0 thread (days): ', hrs_idx
            hrs_idx = hrs_idx + (hrs_freq/24_qb)
        END IF

    END DO

    !< save time-avg debugging output
        !< Tf physical terms:
        Tfconv_tavg = Tfconv_tavg / numSteps / (1-statistic_tic)
        !< radiation components:
        Torad_tavg  = Torad_tavg / numSteps / (1-statistic_tic)
        Tbrad_tavg  = Tbrad_tavg / numSteps / (1-statistic_tic)
        Tfrad_tavg  = Tfrad_tavg / numSteps / (1-statistic_tic)

    CALL output_statistics

    !< write time sequence and display info when simulation finished
    IF (sch%procID .EQ. 0) THEN
        IF (sch%commSize .EQ. 1) THEN
            WRITE(time_sequence_name,'(A)') './output_serial/time_sequence.csv'
        ELSE
            WRITE(time_sequence_name,'(A)') './output_parallel/time_sequence.csv'
        END IF
        OPEN(115,file=time_sequence_name,form='formatted')
        DO i = 1, op_Max
           WRITE(115,'(E32.16,A,1x)',ADVANCE='NO') op_time(i), ','
           WRITE(115,'(1x)')
        END DO
        CLOSE(115)

        print '(a)', repeat('-', 60)
        print *, NEW_LINE(' '), 'num of sub-gird rows on thread: ', xLen
        print *, NEW_LINE(' '), 'num of sub-gird cols on thread: ', yLen
        print *, NEW_LINE(' '), 'model simulation time (days): ', INT(op_time(time_idx))
        print '(a)', repeat('-', 60)
        print *, NEW_LINE(' '), 'stored time sequential info, closing program.'
        print '(a)', repeat('-', 60)

    END IF

END SUBROUTINE iterate_steps

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: scheme_physc
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @varbs[in]: Tb, Tf, sig_c, qvb
! @varbs[io]: To, thetaeb, theta1, qtb, qf
! @varbs[io]: ub, vb, u0, v0, u1, v1
! @simulation_process: 0 >>> 1
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE scheme_physc(step_num)

    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb) :: i, j, k, l
    REAL(dp) :: pos_x, pos_y
    REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp)

    INTEGER(qb), INTENT(IN) :: step_num

    !< define physical terms which need to be re-calculated during each iterations
    REAL(dp) :: qc(0:xLen-1,0:yLen-1)        !< ft convection threshold
    REAL(dp) :: al(0:xLen-1,0:yLen-1)        !< ABL long wave absp
    REAL(dp) :: alf(0:xLen-1,0:yLen-1)       !< ft long wave absp
    REAL(dp) :: F_Tf(0:xLen-1,0:yLen-1)      !< extra force in FT's radiation
    REAL(dp) :: qsat(0:xLen-1,0:yLen-1)      !< saturation value of sea surface
    REAL(dp) :: Ecloud(0:xLen-1,0:yLen-1)    !< cloud top exchange coeff
    REAL(dp) :: Frad_o(0:xLen-1,0:yLen-1)    !< ocean radiation term
    REAL(dp) :: Frad_a(0:xLen-1,0:yLen-1)    !< ABL radiation term
    REAL(dp) :: Frad_f(0:xLen-1,0:yLen-1)    !< FT radiation term
    REAL(dp) :: ss_evap(0:xLen-1,0:yLen-1)   !< sea surface evaporation
    REAL(dp) :: ex_moist(0:xLen-1,0:yLen-1)  !< exchange of moisture quantity
    REAL(dp) :: sens_heat(0:xLen-1,0:yLen-1) !< sensitive heating between ABL & FT
    REAL(dp) :: ft_convec(0:xLen-1,0:yLen-1) !< ft strong convection quantity

    REAL(dp) :: SigcTemp(0:RainSizexLen-1,0:RainSizeyLen-1)
    REAL(dp) :: RainSizeTemp(0:RainSizexLen-1,0:RainSizeyLen-1)

    DO i = 0, xLen-1
    DO j = 0, yLen-1

        !< spatial position
        pos_x = xRef + i*dx
        pos_y = yRef + j*dy

        !< update extra physical variables
        Tf(i,j)        = T0f + T1f*theta1(i,j,0)
        qc(i,j)        = q0f + q1f*(Tf(i,j)-T0c)
        qsat(i,j)      = q0b + q1b*To(i,j,0)
        qsat_b(i,j)    = q0b + q1b*Tb(i,j)
        sig_f(i,j)     = heaviside1(qf(i,j,0)-qc(i,j))
        sig_b(i,j) = sig_c(i,j)*(1-sig_f(i,j))

        !< update solar terms
            F_Tf(i,j) = 0.0_dp
            al(i,j)   = al0 + al1*((qvb(i,j)/qsat_b(i,j))*(1.0_dp-sig_c(i,j))+sig_c(i,j))
            alf(i,j)  = alf0 + alf1*((qf(i,j,0)/qc(i,j))*(1.0_dp-sig_f(i,j))+sig_f(i,j))
            Frad_o(i,j) = S*(1.0_dp-Ac*sig_c(i,j))*(1.0_dp-Af*sig_f(i,j))*(1.0_dp-asf)*(1.0_dp-as) &
                            + boltz*alf(i,j)*(1.0_dp-al(i,j))*Tf(i,j)**4_qb &
                            + boltz*al(i,j)*Tb(i,j)**4_qb &
                            - boltz*To(i,j,0)**4_qb
            Frad_a(i,j) = S*(1.0_dp-Ac*sig_c(i,j))*(1.0_dp-Af*sig_f(i,j))*(1.0_dp-asf)*as &
                            + boltz*alf(i,j)*al(i,j)*Tf(i,j)**4_qb &
                            + boltz*al(i,j)*To(i,j,0)**4_qb &
                            - 2.0_dp*boltz*al(i,j)*Tb(i,j)**4_qb
            Frad_f(i,j) = S*asf*(1.0_dp-Af*sig_f(i,j)) &
                            + S*Ac*sig_c(i,j)*(1.0_dp-Af*sig_f(i,j))*asf*(1.0_dp-asf) &
                            + boltz*alf(i,j)*(1.0_dp-al(i,j))*To(i,j,0)**4_qb &
                            + boltz*alf(i,j)*al(i,j)*Tb(i,j)**4_qb &
                            - 2.0_dp*boltz*alf(i,j)*Tf(i,j)**4_qb + F_Tf(i,j)

        !< solve physical processings
        ss_evap(i,j)   = (qsat(i,j)-qtb(i,j,0)) / tau_e
        sens_heat(i,j) = (To(i,j,0)-Tb(i,j)) / tau_s
        ft_convec(i,j) = sig_f(i,j) * (qf(i,j,0)-qc(i,j)) / tau_q

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !< [2022.Mar]: optimal cloud top mixing methods
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !Ecloud(i,j)    = sig_c(i,j)/tau_t1
            !Ecloud(i,j)    = (sig_c(i,j)/tau_t1) + (sig_f(i,j)/tau_t2)
             Ecloud(i,j)    = sig_c(i,j) * ((1/tau_t1)+(sig_f(i,j)/tau_t2))
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            !ex_moist(i,j)  = Ecloud(i,j) * (qtb(i,j,0)-CTex_qf*qf(i,j,0))
            !ex_moist(i,j)  = Ecloud(i,j) * (qvb(i,j)-CTex_qf*qf(i,j,0))
             ex_moist(i,j)  = Ecloud(i,j) * MAX(qtb(i,j,0)-CTex_qf*qf(i,j,0), 0.0_dp)
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        !< update model variables forward in time
        qtb(i,j,1)     = dt * ( - ex_moist(i,j) + ss_evap(i,j) ) + qtb(i,j,0)
        qf(i,j,1)      = dt * ( + ex_moist(i,j) - ft_convec(i,j) ) + qf(i,j,0)
        To(i,j,1)      = dt * ( - (Lvo/cs)*ss_evap(i,j) - (Cb/Co)*sens_heat(i,j) &
                                + Frad_o(i,j)/Co ) + To(i,j,0)
        thetaeb(i,j,1) = dt * ( - (Lvb/cp)*ex_moist(i,j) + sens_heat(i,j) + (Lvb/cp)*ss_evap(i,j) &
                                + Frad_a(i,j)/Cb ) + thetaeb(i,j,0)
        theta1(i,j,1)  = (dt/T1f) * ( + LHM*(Lvf/cp)*ft_convec(i,j) + Frad_f(i,j)/Cf ) &
                                + theta1(i,j,0)
        ub(i,j,1) = dt * (-(Eu/hb)*(ub(i,j,0)-u0(i,j,0)-SQRT(2.0_dp)*u1(i,j,0))) + ub(i,j,0)
        vb(i,j,1) = dt * (-(Eu/hb)*(vb(i,j,0)-v0(i,j,0)-SQRT(2.0_dp)*v1(i,j,0))) + vb(i,j,0)
        u0(i,j,1) = dt * (+(Eu/HT)*(ub(i,j,0)-u0(i,j,0)-SQRT(2.0_dp)*u1(i,j,0))) + u0(i,j,0)
        v0(i,j,1) = dt * (+(Eu/HT)*(vb(i,j,0)-v0(i,j,0)-SQRT(2.0_dp)*v1(i,j,0))) + v0(i,j,0)
        u1(i,j,1) = dt * (-u1(i,j,0)/tau_R) + u1(i,j,0)
        v1(i,j,1) = dt * (-v1(i,j,0)/tau_R) + v1(i,j,0)

        !< add spatial varying forcing
        To(i,j,1) = dt * (-Fo*SIN(scale_y*pos_y)) + To(i,j,1)

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !< debug stats: radiation segments and Tf EQ details
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            IF (step_num.GT.numSteps*statistic_tic) THEN
                Torad_tavg(i,j) = Torad_tavg(i,j) + Frad_o(i,j)/Co
                Tbrad_tavg(i,j) = Tbrad_tavg(i,j) + Frad_a(i,j)/Cb
                Tfrad_tavg(i,j) = Tfrad_tavg(i,j) + Frad_f(i,j)/Cf
                Tfconv_tavg(i,j) = Tfconv_tavg(i,j) + LHM*(Lvf/cp)*ft_convec(i,j)
            END IF
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    END DO
    END DO

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !< [2022.Aug]: histogram data
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !< [2022.Aug]: rainfall size statistics
    !< [2022.Sep]: rain/dry events duration statistics
    !< [2022.Sep]: shallow cluster duration statistics
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RainSizeTemp(:,:) = 0.0_dp
        DO i = 0, RainSizexLen-1
        DO j = 0, RainSizeyLen-1
        !< solving rainfall size average over unit boxes at CURRENT time step
            DO k = 0, RainBoxxLen-1
            DO l = 0, RainBoxyLen-1
                SigcTemp(i,j) = SigcTemp(i,j) + sig_c(i*RainBoxxLen+k,j*RainBoxyLen+l)
                RainSizeTemp(i,j) = RainSizeTemp(i,j) + &
                    dt*ft_convec(i*RainBoxxLen+k,j*RainBoxyLen+l)
            END DO
            END DO
            SigcTemp(i,j) = SigcTemp(i,j)/RainBoxxLen/RainBoxyLen
            RainSizeTemp(i,j) = RainSizeTemp(i,j)/RainBoxxLen/RainBoxyLen
        !< Histogram data: most and dry events
            !< case 1: continuing rainfall
            IF ((RainSizeTemp(i,j).GT.0.0_dp).AND.(RainTrigger(i,j).EQ.1_qb)) THEN
                DryTime(i,j) = 0.0_dp
                RainTime(i,j) = RainTime(i,j) + 1.0_dp
                RainSize(i,j) = RainSize(i,j) + RainSizeTemp(i,j)
                RainTrigger(i,j) = 1_qb
            !< case 2: rainfall stopped
            ELSEIF ((RainSizeTemp(i,j).LE.0.0_dp).AND.(RainTrigger(i,j).EQ.1_qb)) THEN
                RainTime(i,j) = RainTime(i,j) + 1.0_dp
                RainSize(i,j) = RainSize(i,j) + RainSizeTemp(i,j)
                !< moist spell duration
                IF (RainTime(i,j).GT.RainTimeMIN) THEN
                    RTBinIdx = INT((LOG(RainTime(i,j))-LOG(RainTimeMIN))/RTLogBinSize, qb)
                    RTBinIdx = MIN(RTBinIdx, RTBinNum-1_qb)
                    RainTimeHist(RTBinIdx) = RainTimeHist(RTBinIdx) + 1.0_dp
                END IF
                !< rainfall events size
                IF (RainSize(i,j).GT.RainSizeMIN) THEN
                    !< rainfall size lower threshold
                    IF (RainSize(i,j).LT.RainSizeMAX) THEN
                        RSBinIdx = INT((LOG(RainSize(i,j))-LOG(RainSizeMIN))/RSLogBinSize, qb)
                        !RainBinIdx = INT((RainSize(i,j)-RainSizeMIN)/BinSize, qb)
                    !< rainfall size upper threshold
                    ELSEIF (RainSize(i,j).GE.RainSizeMAX) THEN
                        RSBinIdx = RSBinNum-1_qb
                    END IF
                    !< update rainfall size histogram data
                    RainSizeHist(RSBinIdx) = RainSizeHist(RSBinIdx) + 1.0_dp
                END IF
                !< reset moist spell stats
                RTBinIdx = 0_qb
                RSBinIdx = 0_qb
                DryTime(i,j) = 1.0_dp
                RainTime(i,j) = 0.0_dp
                RainSize(i,j) = 0.0_dp
                RainTrigger(i,j) = 0_qb
            !< case 3: rainfall starts
            ELSEIF ((RainSizeTemp(i,j).GT.0.0_dp).AND.(RainTrigger(i,j).EQ.0_qb)) THEN
                DryTime(i,j) = DryTime(i,j) + 1.0_dp
                !< dry spell duration
                IF (DryTime(i,j).GT.DryTimeMIN) THEN
                    !< ONLY collect data over the warm pool quator
                    IF ((sch%procID.GE.2_qb).AND.(sch%procID.LE.7_dp)) THEN
                        DTBinIdx = INT((LOG(DryTime(i,j))-LOG(DryTimeMIN))/DTLogBinSize, qb)
                        DTBinIdx = MIN(DTBinIdx, DTBinNum-1_qb)
                        DryTimeHist(DTBinIdx) = DryTimeHist(DTBinIdx) + 1.0_dp
                    END IF
                END IF
                !< reset of dry spell stats
                DTBinIdx = 0_qb
                DryTime(i,j) = 0.0_dp
                !< moist trigger on
                RainTime(i,j) = 1.0_dp
                RainSize(i,j) = RainSizeTemp(i,j)
                RainTrigger(i,j) = 1_qb
            !< case 4: continuing drying
            ELSEIF ((RainSizeTemp(i,j).LE.0.0_dp).AND.(RainTrigger(i,j).EQ.0_qb)) THEN
                DryTime(i,j) = DryTime(i,j) + 1.0_dp
                RainTime(i,j) = 0.0_dp
                RainSize(i,j) = 0.0_dp
                RainTrigger(i,j) = 0_qb
            END IF

            !< Histogrma data: shallow cloud duration
            !< ONLY collect data over the cool pool quator
            IF ((sch%procID.GE.12_qb).AND.(sch%procID.LE.17_dp)) THEN
                IF ((SigcTemp(i,j).GT.0_qb).AND.(SigcTrigger(i,j).EQ.0_qb)) THEN
                    !SigcTime(i,j) = 1.0_dp
                    SigcTrigger(i,j) = 1_qb
                ELSEIF ((SigcTemp(i,j).GT.0_qb).AND.(SigcTrigger(i,j).EQ.1_qb)) THEN
                    SigcTime(i,j) = SigcTime(i,j) + 1.0_dp
                    SigcTrigger(i,j) = 1_qb
                ELSEIF ((SigcTemp(i,j).LE.0_qb).AND.(SigcTrigger(i,j).EQ.0_qb)) THEN
                    SigcTime(i,j) = 0.0_dp
                    SigcTrigger(i,j) = 0_qb
                ELSEIF ((SigcTemp(i,j).LE.0_qb).AND.(SigcTrigger(i,j).EQ.1_qb)) THEN
                    SigcTime(i,j) = SigcTime(i,j) + 1.0_dp
                    IF (SigcTime(i,j).GT.SigcTimeMIN) THEN
                        STBinIdx = INT((LOG(SigcTime(i,j))-LOG(SigcTimeMIN))/STLogBinSize, qb)
                        STBinIdx = MIN(STBinIdx, STBinNum-1_qb)
                        SigcTimeHist(STBinIdx) = SigcTimeHist(STBinIdx) + 1.0_dp
                    END IF
                    STBinIdx = 0_qb
                    SigcTime(i,j) = 0.0_dp
                    SigcTrigger(i,j) = 0_qb
                END IF
            END IF

        END DO
        END DO
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

END SUBROUTINE scheme_physc


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: scheme_fluid
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @varbs[in]: Tb
! @varbs[io]: ub, vb, u0, v0, u1, v1, theta1, qf
! @varbs[io]: To, thetaeb (diff-only)
! @simulation_process: 1 >>> 2
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE scheme_fluid

    USE MPI
    USE INITIALIZE
    USE PARALLEL_TOOLS
    USE PARALLEL_STRUCTS

    IMPLICIT NONE

    INTEGER(qb) :: i, j

    !< define spectral variables' array
    COMPLEX(dp), DIMENSION(3,1) :: Lnsys_varb
    COMPLEX(dp), DIMENSION(3,1) :: Lnsys_force
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1) :: Freq_Tb
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1) :: Freq_To
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1) :: Freq_thetaeb
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_ub
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_vb
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_u0
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_v0
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_u1
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_v1
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_qf
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_qb
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_theta1
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: psi_0
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: psi_b
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: phi_0
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: phi_b
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: phi_max
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: phi_min
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1,2) :: Freq_w1

    !< advection part in nonlinear(linearized as qfhat) surface flux
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1) :: Freq_advqf1
    COMPLEX(dp), DIMENSION(0:xLen-1,0:yLen-1) :: Freq_advqf2

    Freq_Tb(:,:) = Tb(:,:)
    Freq_To(:,:) = To(:,:,1)
    Freq_ub(:,:,1) = ub(:,:,1)
    Freq_vb(:,:,1) = vb(:,:,1)
    Freq_u0(:,:,1) = u0(:,:,1)
    Freq_v0(:,:,1) = v0(:,:,1)
    Freq_u1(:,:,1) = u1(:,:,1)
    Freq_v1(:,:,1) = v1(:,:,1)
    Freq_qf(:,:,1) = qf(:,:,1)
    Freq_qb(:,:,1) = qtb(:,:,1)
    Freq_thetaeb(:,:) = thetaeb(:,:,1)
    Freq_theta1(:,:,1) = theta1(:,:,1)

    ! 2d zfft for fluid dynamics varbs
    CALL PARALLEL_ZFFT(Freq_To(:,:), sch)
    CALL PARALLEL_ZFFT(Freq_Tb(:,:), sch)
    CALL PARALLEL_ZFFT(Freq_ub(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_vb(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_u0(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_v0(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_u1(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_v1(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_qf(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_qb(:,:,1), sch)
    CALL PARALLEL_ZFFT(Freq_thetaeb(:,:), sch)
    CALL PARALLEL_ZFFT(Freq_theta1(:,:,1), sch)

    DO i = 0, xLen-1
    DO j = 0, yLen-1

        !< special modification for zero Fourier modes
        IF ((ZABS(spec_dx(i,j)).LT.1E-14).AND.(ZABS(spec_dy(i,j)).LT.1E-14)) THEN

            Freq_w1(i,j,1) = DCMPLX(0.0, 0.0)
            Freq_ub(i,j,2) = Freq_ub(i,j,1) * fluid_coeffs_1
            Freq_vb(i,j,2) = Freq_vb(i,j,1) * fluid_coeffs_1
            Freq_u0(i,j,2) = Freq_u0(i,j,1)
            Freq_v0(i,j,2) = Freq_v0(i,j,1)
            Freq_u1(i,j,2) = Freq_u1(i,j,1)
            Freq_v1(i,j,2) = Freq_v1(i,j,1)
            Freq_theta1(i,j,2) = Freq_theta1(i,j,1)

        !< general Fourier modes
        ELSE

            !< psi = dy*u - dx*v / dxx + dyy
            psi_0(i,j,1) = (spec_dx(i,j)*Freq_v0(i,j,1)-spec_dy(i,j)*Freq_u0(i,j,1)) &
                            / -(spec_dxx(i,j)+spec_dyy(i,j))
            psi_b(i,j,1) = (spec_dx(i,j)*Freq_vb(i,j,1)-spec_dy(i,j)*Freq_ub(i,j,1)) &
                            / -(spec_dxx(i,j)+spec_dyy(i,j))

            !< phi = dx*u + dy*v / dxx + dyy
            phi_0(i,j,1) = (spec_dx(i,j)*Freq_u0(i,j,1)+spec_dy(i,j)*Freq_v0(i,j,1)) &
                            / +(spec_dxx(i,j)+spec_dyy(i,j))
            phi_b(i,j,1) = (spec_dx(i,j)*Freq_ub(i,j,1)+spec_dy(i,j)*Freq_vb(i,j,1)) &
                            / +(spec_dxx(i,j)+spec_dyy(i,j))

            !< c.o.v: phi_+ / phi_- / Freq_w1
            phi_max(i,j,1) = hb*phi_b(i,j,1) + HT*phi_0(i,j,1)
            phi_min(i,j,1) = hb*phi_b(i,j,1) - HT*phi_0(i,j,1)
            Freq_w1(i,j,1) = spec_dx(i,j)*Freq_u1(i,j,1)+spec_dy(i,j)*Freq_v1(i,j,1)

            !< variables and forcing terms in spectral linear system
            Lnsys_varb(1,1)  = Freq_w1(i,j,1)
            Lnsys_varb(2,1)  = Freq_theta1(i,j,1)
            Lnsys_varb(3,1)  = phi_min(i,j,1)
            Lnsys_force(1,1) = (0.0_dp, 0.0_dp)
            Lnsys_force(2,1) = -fluid_coeffs_2*phi_max(i,j,1)*(spec_dxx(i,j)+spec_dyy(i,j))
            Lnsys_force(3,1) = fluid_coeffs_3*phi_max(i,j,1)

            !< solve spectral linear system semi-analytically
            CALL lnsys_solver(Lnsys_varb, Lnsys_force, i, j)

            !< update stream and potential arraies
            phi_min(i,j,2) = Lnsys_varb(3,1)
            phi_max(i,j,2) = phi_max(i,j,1)
            psi_0(i,j,2)   = psi_0(i,j,1)
            psi_b(i,j,2)   = psi_b(i,j,1) * fluid_coeffs_1
            phi_0(i,j,2)   = (phi_max(i,j,2)-phi_min(i,j,2))/2.0_dp/HT
            phi_b(i,j,2)   = (phi_max(i,j,2)+phi_min(i,j,2))/2.0_dp/hb

            !< 2d inverse Helmholtz decomposition
            Freq_ub(i,j,2) = spec_dx(i,j)*phi_b(i,j,2) + spec_dy(i,j)*psi_b(i,j,2)
            Freq_vb(i,j,2) = spec_dy(i,j)*phi_b(i,j,2) - spec_dx(i,j)*psi_b(i,j,2)
            Freq_u0(i,j,2) = spec_dx(i,j)*phi_0(i,j,2) + spec_dy(i,j)*psi_0(i,j,2)
            Freq_v0(i,j,2) = spec_dy(i,j)*phi_0(i,j,2) - spec_dx(i,j)*psi_0(i,j,2)
            Freq_w1(i,j,2) = Lnsys_varb(1,1)
            Freq_u1(i,j,2) = (spec_dx(i,j)*Freq_w1(i,j,2) + Freq_u1(i,j,1)*(spec_dy(i,j)**2_qb) &
                              - Freq_v1(i,j,1)*spec_dx(i,j)*spec_dy(i,j)) &
                             / (spec_dx(i,j)**2_qb+spec_dy(i,j)**2_qb)
            Freq_v1(i,j,2) = (spec_dy(i,j)*Freq_w1(i,j,2) + Freq_v1(i,j,1)*(spec_dx(i,j)**2_qb) &
                              - Freq_u1(i,j,1)*spec_dx(i,j)*spec_dy(i,j)) &
                             / (spec_dx(i,j)**2_qb+spec_dy(i,j)**2_qb)
            Freq_theta1(i,j,2) = Lnsys_varb(2,1)

        END IF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !< evolve in frequency space, setting qf as qfhat
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Freq_advqf1(i,j) = spec_dx(i,j)*Freq_u0(i,j,1)+spec_dy(i,j)*Freq_v0(i,j,1)
        Freq_advqf2(i,j) = spec_dx(i,j)*Freq_u1(i,j,1)+spec_dy(i,j)*Freq_v1(i,j,1)
        Freq_qf(i,j,2) = Freq_qf(i,j,1) - dt*qfhat*alpha_adv*Freq_advqf2(i,j) &
            - dt*qfhat*(1.0_dp-(1/Q0_hat))*Freq_advqf1(i,j)
        Freq_qb(i,j,2) = Freq_qb(i,j,1)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    END DO
    END DO

    !< add eddy-diffusion to model variables (except moisture variables)
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        Freq_To(i,j) = Freq_To(i,j)*ZEXP(-diff_coeff_To(i,j)*dt)
        Freq_u0(i,j,2) = Freq_u0(i,j,2)*ZEXP(-diff_coeff_u0(i,j)*dt)
        Freq_v0(i,j,2) = Freq_v0(i,j,2)*ZEXP(-diff_coeff_v0(i,j)*dt)
        Freq_u1(i,j,2) = Freq_u1(i,j,2)*ZEXP(-diff_coeff_u1(i,j)*dt)
        Freq_v1(i,j,2) = Freq_v1(i,j,2)*ZEXP(-diff_coeff_v1(i,j)*dt)
        Freq_ub(i,j,2) = Freq_ub(i,j,2)*ZEXP(-diff_coeff_ub(i,j)*dt)
        Freq_vb(i,j,2) = Freq_vb(i,j,2)*ZEXP(-diff_coeff_vb(i,j)*dt)
        Freq_thetaeb(i,j) = Freq_thetaeb(i,j)*ZEXP(-diff_coeff_thetaeb(i,j)*dt)
        Freq_theta1(i,j,2) = Freq_theta1(i,j,2)*ZEXP(-diff_coeff_theta1(i,j)*dt)
    END DO
    END DO

    !< inverse 2dfft
    CALL PARALLEL_IFFT(Freq_To(:,:), sch)
    CALL PARALLEL_IFFT(Freq_ub(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_vb(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_u0(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_v0(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_u1(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_v1(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_qf(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_qb(:,:,2), sch)
    CALL PARALLEL_IFFT(Freq_thetaeb(:,:), sch)
    CALL PARALLEL_IFFT(Freq_theta1(:,:,2), sch)

    !< update fluid dynamics variables
    To(:,:,2) = REAL(Freq_To(:,:), dp)
    qf(:,:,2) = REAL(Freq_qf(:,:,2), dp)
    ub(:,:,2) = REAL(Freq_ub(:,:,2), dp)
    vb(:,:,2) = REAL(Freq_vb(:,:,2), dp)
    u0(:,:,2) = REAL(Freq_u0(:,:,2), dp)
    v0(:,:,2) = REAL(Freq_v0(:,:,2), dp)
    u1(:,:,2) = REAL(Freq_u1(:,:,2), dp)
    v1(:,:,2) = REAL(Freq_v1(:,:,2), dp)
    qtb(:,:,2) = REAL(Freq_qb(:,:,2), dp)
    theta1(:,:,2) = REAL(Freq_theta1(:,:,2), dp)
    thetaeb(:,:,2) = REAL(Freq_thetaeb(:,:), dp)


END SUBROUTINE scheme_fluid


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: scheme_randn
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @varbs[io]: qtb, qf
! @simulation_process: 2 >>> 3
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE scheme_randn

    USE MPI
    USE INITIALIZE
    USE PARALLEL_TOOLS
    USE PARALLEL_STRUCTS

    IMPLICIT NONE

    REAL(dp), PARAMETER :: pi_dp = 4.0_dp * ATAN(1.0_dp)

    INTEGER(qb) :: i, j
    REAL(dp) :: rand1, rand2
    REAL(dp) :: temp0_dp, temp1_dp, temp2_dp

    REAL(dp) :: phys_noise_FT(0:xLen-1,0:yLen-1)
    REAL(dp) :: phys_noise_ABL(0:xLen-1,0:yLen-1)
    COMPLEX(dp) :: freq_noise_FT(0:xLen-1,0:yLen-1)
    COMPLEX(dp) :: freq_noise_ABL(0:xLen-1,0:yLen-1)
    COMPLEX(dp) :: white_noise_1(0:xLen-1,0:yLen-1)
    COMPLEX(dp) :: white_noise_2(0:xLen-1,0:yLen-1)
    COMPLEX(dp) :: Freq_qtb(0:xLen-1,0:yLen-1,2:3)
    COMPLEX(dp) :: Freq_qf(0:xLen-1,0:yLen-1,2:3)

    !< 2d zfft moisture variables
    Freq_qtb(:,:,2) = qtb(:,:,2)
    Freq_qf(:,:,2)  = qf(:,:,2)
    CALL PARALLEL_ZFFT(Freq_qtb(:,:,2), sch)
    CALL PARALLEL_ZFFT(Freq_qf(:,:,2), sch)

    ! white noise in ABL equation: using Box-Mueller transform generate Gaussian randomness
    ! flag_randn = 1 >>> physical white noise / flag_randn = 2 >>> spectral white noise
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        CALL RANDOM_NUMBER(rand1)
        CALL RANDOM_NUMBER(rand2)
        IF (flag_randn.EQ.1) THEN
            temp0_dp = SQRT(-2.0_dp*LOG(rand1))*COS(2.0_dp*pi_dp*rand2)
            phys_noise_ABL(i,j) = 0.0_dp + temp0_dp*SQRT(phys_var_ABL)
        ELSE IF (flag_randn.EQ.2) THEN
            temp1_dp = SQRT(-2.0_dp*LOG(rand1))*COS(2.0_dp*pi_dp*rand2)
            temp2_dp = SQRT(-2.0_dp*LOG(rand1))*SIN(2.0_dp*pi_dp*rand2)
            temp1_dp = 0.0_dp + temp1_dp*SQRT(REALPART(freq_var_ABL(i,j)))
            temp2_dp = 0.0_dp + temp2_dp*SQRT(REALPART(freq_var_ABL(i,j)))
            freq_noise_ABL(i,j) = 1.0_dp/SQRT(2.0_dp)*CMPLX(temp1_dp, temp2_dp, dp)
        END IF
    END DO
    END DO

    ! white noise in FT equation: using Box-Mueller transform generate Gaussian randomness
    ! flag_randn = 1 >>> physical white noise / flag_randn = 2 >>> spectral white noise
    DO i = 0, xLen-1
    DO j = 0, yLen-1
        CALL RANDOM_NUMBER(rand1)
        CALL RANDOM_NUMBER(rand2)
        IF (flag_randn.EQ.1) THEN
            temp0_dp = SQRT(-2.0_dp*LOG(rand1))*COS(2.0_dp*pi_dp*rand2)
            phys_noise_FT(i,j) = 0.0_dp + temp0_dp*SQRT(phys_var_FT)
        ELSE IF (flag_randn.EQ.2) THEN
            temp1_dp = SQRT(-2.0_dp*LOG(rand1))*COS(2.0_dp*pi_dp*rand2)
            temp2_dp = SQRT(-2.0_dp*LOG(rand1))*SIN(2.0_dp*pi_dp*rand2)
            temp1_dp = 0.0_dp + temp1_dp*SQRT(REALPART(freq_var_FT(i,j)))
            temp2_dp = 0.0_dp + temp2_dp*SQRT(REALPART(freq_var_FT(i,j)))
            freq_noise_FT(i,j) = 1.0_dp/SQRT(2.0_dp)*CMPLX(temp1_dp, temp2_dp, dp)
        END IF
    END DO
    END DO

    !< physical white noise case
    IF (flag_randn.EQ.1) THEN
        white_noise_1(:,:) = phys_noise_ABL(:,:)
        white_noise_2(:,:) = phys_noise_FT(:,:)
        CALL PARALLEL_ZFFT(white_noise_1, sch)
        CALL PARALLEL_ZFFT(white_noise_2, sch)
        DO i = 0, xLen-1
        DO j = 0, yLen-1
            IF (ZABS(diff_coeff_qtb(i,j)).LT.1E-15) THEN
                white_noise_1(i,j) = white_noise_1(i,j) * dt
            ELSE
                white_noise_1(i,j) = (1.0_dp - ZEXP(-diff_coeff_qtb(i,j)*dt)) &
                                     * white_noise_1(i,j) / diff_coeff_qtb(i,j)
            END IF
            IF (ZABS(diff_coeff_qf(i,j)).LT.1E-15) THEN
                white_noise_2(i,j) = white_noise_2(i,j) * dt
            ELSE
                white_noise_2(i,j) = (1.0_dp - ZEXP(-diff_coeff_qf(i,j)*dt)) &
                                     * white_noise_2(i,j) / diff_coeff_qf(i,j)
            END IF
            Freq_qtb(i,j,3) = Freq_qtb(i,j,2)*ZEXP(-diff_coeff_qtb(i,j)*dt) + white_noise_1(i,j)
            Freq_qf(i,j,3)  = Freq_qf(i,j,2)*ZEXP(-diff_coeff_qf(i,j)*dt) + white_noise_2(i,j)
        END DO
        END DO
    END IF

    !< spectral white noise case
    IF (flag_randn.EQ.2) THEN
        white_noise_1(:,:) = freq_noise_ABL(:,:)
        white_noise_2(:,:) = freq_noise_FT(:,:)
        DO i = 0, xLen-1
        DO j = 0, yLen-1
            Freq_qtb(i,j,3) = Freq_qtb(i,j,2)*ZEXP(-diff_coeff_qtb(i,j)*dt) + white_noise_1(i,j)
            Freq_qf(i,j,3)  = Freq_qf(i,j,2)*ZEXP(-diff_coeff_qf(i,j)*dt) + white_noise_2(i,j)
        END DO
        END DO
    END IF

    !< update moisture variables back in physical space
    CALL PARALLEL_IFFT(Freq_qf(:,:,3), sch)
    CALL PARALLEL_IFFT(Freq_qtb(:,:,3), sch)
    qf(:,:,3)  = REAL(Freq_qf(:,:,3), dp)
    qtb(:,:,3) = REAL(Freq_qtb(:,:,3), dp)

    !< do NOT update other model variables
    ub(:,:,3) = ub(:,:,2)
    vb(:,:,3) = vb(:,:,2)
    u0(:,:,3) = u0(:,:,2)
    v0(:,:,3) = v0(:,:,2)
    u1(:,:,3) = u1(:,:,2)
    v1(:,:,3) = v1(:,:,2)
    To(:,:,3) = To(:,:,2)
    theta1(:,:,3) = theta1(:,:,2)
    thetaeb(:,:,3) = thetaeb(:,:,2)

END SUBROUTINE scheme_randn

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: scheme_cloud
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @varbs[in]: thetaeb, qtb
! @varbs[ot]: sig_c, qvb, Tb, Tf
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE scheme_cloud

    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb) :: i, j

    REAL(dp) :: qcloud_implict(0:xLen-1,0:yLen-1)

    DO i = 0, xLen-1
    DO j = 0, yLen-1
        qcloud_implict(i,j) = (q0b+q1b*thetaeb(i,j,op_idx)) / (1.0_dp+(Lvb/cp)*q1b)
        sig_c(i,j) = cloud_type(qtb(i,j,op_idx), qcloud_implict(i,j))
        qvb(i,j) = sig_c(i,j)*qcloud_implict(i,j)+(1.0_dp-sig_c(i,j))*qtb(i,j,op_idx)
        Tb(i,j) = thetaeb(i,j,op_idx) - (Lvb/cp)*qvb(i,j)
    END DO
    END DO

END SUBROUTINE scheme_cloud

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: update_varbs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @model_varbs(ip_idx) = model_varbs(op_idx)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE update_varbs(step_num)

    USE INITIALIZE

    IMPLICIT NONE

    INTEGER(qb) :: i, j
    INTEGER(qb), INTENT(IN) :: step_num

    ub(:,:,ip_idx) = ub(:,:,op_idx)
    vb(:,:,ip_idx) = vb(:,:,op_idx)
    u0(:,:,ip_idx) = u0(:,:,op_idx)
    v0(:,:,ip_idx) = v0(:,:,op_idx)
    u1(:,:,ip_idx) = u1(:,:,op_idx)
    v1(:,:,ip_idx) = v1(:,:,op_idx)
    To(:,:,ip_idx) = To(:,:,op_idx)
    qf(:,:,ip_idx) = qf(:,:,op_idx)
    qtb(:,:,ip_idx) = qtb(:,:,op_idx)
    theta1(:,:,ip_idx) = theta1(:,:,op_idx)
    thetaeb(:,:,ip_idx) = thetaeb(:,:,op_idx)

    IF (step_num.GT.numSteps*statistic_tic) THEN

        To_tavg = (To_tavg*tavg_count+To(:,:,op_idx))/(tavg_count+1_qb)
        Tb_tavg = (Tb_tavg*tavg_count+Tb)/(tavg_count+1_qb)
        Tf_tavg = (Tf_tavg*tavg_count+Tf)/(tavg_count+1_qb)
        u1_tavg = (u1_tavg*tavg_count+u1(:,:,op_idx))/(tavg_count+1_qb)
        v1_tavg = (v1_tavg*tavg_count+v1(:,:,op_idx))/(tavg_count+1_qb)
        qf_tavg = (qf_tavg*tavg_count+qf(:,:,op_idx))/(tavg_count+1_qb)
        qtb_tavg = (qtb_tavg*tavg_count+qtb(:,:,op_idx))/(tavg_count+1_qb)
        Scloud_tavg = (Scloud_tavg*tavg_count+sig_c)/(tavg_count+1_qb)
        Dcloud_tavg = (Dcloud_tavg*tavg_count+sig_f)/(tavg_count+1_qb)
        Bcloud_tavg = (Bcloud_tavg*tavg_count+sig_b)/(tavg_count+1_qb)

        tavg_count = tavg_count + 1_qb

    END IF

    DO i = 0, xLen-1
    DO j = 0, yLen-1
        Scloud_rate(step_num) = Scloud_rate(step_num) + sig_c(i,j)
        Dcloud_rate(step_num) = Dcloud_rate(step_num) + sig_f(i,j)
        Bcloud_rate(step_num) = Bcloud_rate(step_num) + sig_b(i,j)
    END DO
    END DO
    Scloud_rate(step_num) = Scloud_rate(step_num) / xLen / yLen
    Dcloud_rate(step_num) = Dcloud_rate(step_num) / xLen / yLen
    Bcloud_rate(step_num) = Bcloud_rate(step_num) / xLen / yLen

END SUBROUTINE update_varbs

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: Lnsys_solver
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @param[in] mode_x: index of Fourier modes in 1st dimension.
! @param[in] mode_y: index of Fourier modes in 2nd dimension.
! @varbs[in] lnsys_force: RHS forcing term in linear system.
! @varbs[io] lnsys_varb:  model variables updated by linear system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE Lnsys_solver(lnsys_varb, lnsys_force, mode_x, mode_y)

    USE INITIALIZE

    IMPLICIT NONE

    !< mode_x/y: index number for specific Fourier modes
    INTEGER(qb), INTENT(IN) :: mode_x, mode_y
    !< Lnsys_varb: [Freq_w1, Freq_theta1, phi_min]
    COMPLEX(dp), DIMENSION(3,1), INTENT(INOUT) :: Lnsys_varb
    COMPLEX(dp), DIMENSION(3,1), INTENT(IN) :: Lnsys_force
    !< Lnsys_semi: semi-analytical solution in frequence space
    COMPLEX(dp), DIMENSION(3,1) :: Lnsys_semi
    COMPLEX(dp), DIMENSION(3,1) :: Lnsys_semi1
    COMPLEX(dp), DIMENSION(3,1) :: Lnsys_semi2

    !< semi_solution I = V * SIG1 * inv(V) * Force
    Lnsys_semi1 = MATMUL(Evec_inv(mode_x,mode_y,:,:), Lnsys_force)
    Lnsys_semi1 = MATMUL(Eval_ex1(mode_x,mode_y,:,:), Lnsys_semi1)
    Lnsys_semi1 = MATMUL(Evec(mode_x,mode_y,:,:), Lnsys_semi1)

    !< semi_solution II = V * SIG2 * inv(V) * varb
    Lnsys_semi2 = MATMUL(Evec_inv(mode_x,mode_y,:,:), Lnsys_varb)
    Lnsys_semi2 = MATMUL(Eval_ex2(mode_x,mode_y,:,:), Lnsys_semi2)
    Lnsys_semi2 = MATMUL(Evec(mode_x,mode_y,:,:), Lnsys_semi2)

    !< update lnr system variables based on semi-analytical solution
    Lnsys_semi = Lnsys_semi1 + Lnsys_semi2
    Lnsys_varb = Lnsys_semi

END SUBROUTINE Lnsys_solver

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! [SBR]: cloud_type
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! @varbs[in]: qtb, qcloud_implict
! @varbs[ot]: sig_c
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION cloud_type(qtb_1, qcloud_1) RESULT(output0)

    USE INITIALIZE

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: qtb_1, qcloud_1

    output0 = 0.0_dp

    IF (flag_cloud.EQ.1_qb) THEN
        output0 = heaviside1(qtb_1-qcloud_1)
    ELSE IF (flag_cloud.EQ.2_qb) THEN
        output0 = EXP((qtb_1-qcloud_1)*scale_cloud_sig)/(1.0_dp+EXP((qtb_1-qcloud_1)*scale_cloud_sig))
    END IF

END FUNCTION cloud_type

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  [FUNCTION]: heaviside
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
REAL(dp) FUNCTION heaviside1(x) RESULT(output1)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: x

    IF (x .ge. 0.0_dp) THEN
        output1 = 1.0_dp
    ELSE
        output1 = 0.0_dp
    END IF

END FUNCTION heaviside1


END MODULE ITERATIONS
