! chiralpl
! version 0
! 2024 Louis Minion
! This program uses the Frenkel-Holstein Hamiltonian to model a 3D lattice of organic chromophores
! It uses the multiparticle basis set of Philpott, truncating to one- and two-particle states only (two-particle approximation)
! The Hamiltonian is first constructed, then projected into its eigenbasis.
! Properties are then calculated from the eigenstates and vectors.
! Hartree units are used throughout.


! CHECK IF OSC STRENGTHS STORE THEIR EVAL INSTEAD OF JUST INDEX, COULD BE PLOTTING AT INCORRECT ENERGY

program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    use index
    use franckcondon
    use hamiltonian
    use disorder
    use spectra
    implicit none
    integer(wp) :: j, o
    real(wp) :: estimated_RAM,begin_disorder_avgtime,end_disorder_avgtime
    character*8  :: date_now
    character*10 :: time_now
    external :: dsyevr, dlamch
    !declare variables


    print*, '*****************************************'
    print*, '                 chiralpl'
    print*, '*****************************************'
    print*, 'Louis Minion 2024'
    call date_and_time(DATE=date_now,TIME=time_now)
    write(*,'((a),(a))') 'Program Started ', date_now
    ! do parameter setting

    print*, 'Reading input'
    call readInput()
    call LatticeIndex()
    if ( bool_one_particle_states ) call oneParticleIndex()
    if ( bool_two_particle_states ) call twoParticleIndex()
    print*, 'Precalculating vibrational overlap integrals'
    call calcFranckCondonTables()
    call dipole_moment() ! precalc dipole moment vectors for each molecule
    ! Export dipole_moment array to draw pics of system
    write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_dipoles.dat')
    eval_out_f = trim(eval_out_f)
    open(unit=33, file=eval_out_f)
    write(33,'(A,I0,A,I0,A,I0)') 'LATTICE DIMENSIONS (X,Y,Z) ',lattice_dimx,',',lattice_dimy,',',lattice_dimz
    do j = 1, size(mu_xyz,dim=1)
        write(33, '(*(F12.8 : ", "))') mu_xyz(j, 1),mu_xyz(j, 2),mu_xyz(j, 3)
    end do
    close(33)
    print*,'Hamiltonian size()','(',general_counter,general_counter,')'
    estimated_RAM = ((((1.0_wp*general_counter)**2)*64)/(8.0_wp*1.0E9_wp))
    if ( estimated_RAM> 8.0_wp) then
        print*,estimated_RAM,'GB virtual mem requested, higher than the available 8GB'
        write(*,*) '<STOP> Estimated need for RAM greater than that available, reduce size of basis.'
        STOP
    else
        print*,estimated_RAM, 'GB RAM requested'
    end if
    if ( .not. allocated(H)) then
         allocate(H(general_counter,general_counter))
    end if
    if ( .not. allocated( EVAL ) ) then
        allocate( EVAL( general_counter ) )
    end if
    allocate(diagonal_disorder_offsets(lattice_count))
    call construct_covariance_matrix()
    call cholesky_decomp(A_covar,lattice_count)
    allocate(xpl_osc(max_vibs+1,size(EVAL)))
    allocate(ypl_osc(max_vibs+1,size(EVAL)))
    allocate(abs_osc_strengths_x(general_counter))
    allocate(abs_osc_strengths_y(general_counter))
    allocate(xpl_osc_configavg(max_vibs+1,size(EVAL)))
    allocate(ypl_osc_configavg(max_vibs+1,size(EVAL)))
    allocate(abs_osc_strengths_x_configavg(general_counter))
    allocate(abs_osc_strengths_y_configavg(general_counter))
    call cpu_time(begin_disorder_avgtime)
    do o = 1, configs
        EVAL = 0.0_wp
        H = 0.0_wp
        diagonal_disorder_offsets = 0.0_wp
        call draw_multivar_distr(lattice_count,diagonal_disorder_offsets,sigma,A_covar)

        if ( bool_one_particle_states ) call build1particleHamiltonian()
        if ( bool_two_particle_states ) call build2particleHamiltonian()
        if ( bool_one_particle_states .and. bool_two_particle_states ) call build1particle2particleHamiltonian()
        
        if (sigma .eq. 0.0_wp) then
            if (H_out .eq. .true.) then
                write(*,*) 'Writing out Hamiltonian to File'
                open(unit=10, file='h_.csv')
                do j = 1, size(H,dim=1)
                    write(10, '(*(F12.8 : ", "))') H(:, j)
                end do
                close(10)
            else
                write(*,*) 'Skipping saving Hamiltonian' 
            end if
        end if
                
        call Diagonalize(H,'A',general_counter, EVAL, EVAL_COUNT, IU)

    ! Save eigenvalue logic (only if sigma 0)
        if (sigma .eq. 0.0_wp) then
            if (save_evals .eq. .true.) then
                write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_EVALS.csv')
                eval_out_f = trim(eval_out_f)
                write(*,'(a,a)') 'Writing eigenvalues to ', eval_out_f
                open(unit=11, file=eval_out_f)
                do j = 1, size(EVAL,dim=1)
                    write(11, '(*(F12.8 : ", "))') EVAL(j)
                end do
                close(11)
            end if
            if (save_evecs .eq. .true.) then
                write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_EVECS.csv')
                eval_out_f = trim(eval_out_f)
                write(*,*) 'Writing out eigenvectors to File'
                open(unit=15, file=eval_out_f)
                do j = 1, size(H,dim=1)
                    write(15, '(*(F12.8 : ", "))') H(:, j) ! H now matrix of eigenvectors
                end do
                close(15)
            else
                write(*,*) 'Skipping saving eigenvectors' 
            end if
        end if

        ! print*,'Calculating PL oscillator strengths for iteration', o
        ! write(*,'(A,F14.5,A)') '0K: EMISSION FROM LOWEST EXCITON WITH ENERGY:',(EVAL(1)*hw),'cm^-1'
        call pl()
        xpl_osc_configavg = xpl_osc_configavg + xpl_osc
        ypl_osc_configavg = ypl_osc_configavg + ypl_osc
        call absorption()
        abs_osc_strengths_x_configavg = abs_osc_strengths_x_configavg + abs_osc_strengths_x
        abs_osc_strengths_y_configavg = abs_osc_strengths_y_configavg + abs_osc_strengths_y
    end do
    call cpu_time(end_disorder_avgtime)
    print*, 'Configurational average finished in',end_disorder_avgtime-begin_disorder_avgtime,'seconds'
    xpl_osc_configavg = xpl_osc_configavg/configs
    ypl_osc_configavg = ypl_osc_configavg/configs
    abs_osc_strengths_x_configavg = abs_osc_strengths_x_configavg/configs
    abs_osc_strengths_y_configavg = abs_osc_strengths_y_configavg/configs
    call pl_output()
    write(*,*)'Calculating PL spectrum'
    call calc_pl_spec()
    call write_pl()
    ! print*,'Calculating Abs oscillator strengths'
    write(*,*)'Calculating Abs spectrum'
    call calc_abs_spec()
    call write_abs()
    ! call absorption

    ! ! Need to setup dipoles in chiral manner
    ! call photoluminescence
    call cpl()
    call calc_cpl_spec()
    call write_cpl()
    ! call circularlypolarizedluminescence


    ! call outputData

end program



subroutine readInput()
    use variables
    character*255 :: fname
    character*100 :: buffer, label
    integer(wp) :: errstat
    logical :: exists
    integer(wp) :: io_stat, file_no, line_no, pos_space
    parameter (file_no=87)
    call get_command_argument(1, fname, status=errstat )
    if (errstat .ne. 0) then
        write(*,*) 'NO INPUT FILE SPECIFIED'
        write(*,*) '<STOP> Input file required.'
        STOP
    end if
    inquire( file=trim(fname), exist=exists)
    if ( .not. exists ) then
        write(*,*) 'INPUT FILE DOES NOT EXIST'
        write(*,*) '<STOP> Non-existent input file.'
        STOP
    end if
    write(*,*) trim(fname)
    io_stat = 0
    line_no = 0
    open( unit = file_no, file = fname, status = 'old', action = 'read' )
    do while (io_stat .eq. 0)
        read( file_no, '(a)', iostat=io_stat ) buffer
        if (io_stat .eq. 0) then
            line_no = line_no + 1
            pos_space = scan( buffer, ' ' )
            label = buffer( 1:pos_space )
            buffer = buffer( pos_space + 1:)
            if ( label(1:1) == '#' ) cycle ! Enable comments with line start #
            write(*,"(*(a))") '->',trim(label), ' ', trim(buffer)
            select case ( label )
            case ( 'INPUT_NAME' )
                read(buffer, *, iostat=io_stat) INPUT_NAME
            case ( 'LATTICE_DIMX' )
                read(buffer, *, iostat=io_stat) lattice_dimx
            case ( 'LATTICE_DIMY' )
                read(buffer, *, iostat=io_stat) lattice_dimy
            case ( 'LATTICE_DIMZ' )
                read(buffer, *, iostat=io_stat) lattice_dimz
            case ( 'X_SPACING')
                read(buffer, *, iostat=io_stat) x_spacing
            case ( 'Y_SPACING')
                read(buffer, *, iostat=io_stat) y_spacing
            case ( 'Z_SPACING')
                read(buffer, *, iostat=io_stat) z_spacing
            case ( 'MANUAL_COUPLING' )
                read(buffer, *, iostat=io_stat) manual_coupling
            case ( 'NEAREST_NEIGHBOURS' )
                read(buffer, *, iostat=io_stat) n_nearest_neighbour
            case ( 'TWIST_ANGLE' )
                read(buffer, *, iostat=io_stat) phi
            case ( 'ONE_PARTICLE_STATES' )
                read(buffer, *, iostat=io_stat) bool_one_particle_states
            case ( 'TWO_PARTICLE_STATES' )
                read(buffer, *, iostat=io_stat) bool_two_particle_states
            case ( 'MAX_VIB')
                read(buffer, *, iostat=io_stat) max_vibs
            case ( 'HUANG_RHYS_NEUTRAL')
                read(buffer, *, iostat=io_stat) lambda_neutral
            case ( 'MONOMER_TRANSITION')
                read(buffer, *, iostat=io_stat) w00
            case ( 'TRANSITION_DIPOLE' )
                read(buffer, *, iostat=io_stat) mu_0
            case ( 'TEMPERATURE' )
                read(buffer, *, iostat=io_stat) temp
            case ( 'JCOULX' )
                read(buffer, *, iostat=io_stat) JCoulx
            case ( 'JCOULY' )
                read(buffer, *, iostat=io_stat) JCouly
            case ( 'JCOULZ' )
                read(buffer, *, iostat=io_stat) JCoulz
            case ( 'LINEWIDTH' )
                read(buffer, *,iostat=io_stat) lw
            case ('WRITE_HAMILTONIAN_OUT')
                read(buffer, *, iostat=io_stat) H_out
            case ('SAVE_EIGENVALS')
                read(buffer, *, iostat=io_stat) save_evals
            case ('SAVE_EIGENVECS')
                read(buffer, *, iostat=io_stat) save_evecs
            case ( 'CORRELATION_LENGTH' )
                read(buffer, *, iostat=io_stat) l0
            case ( 'NUMBER_CONFIGURATIONS' )
                read(buffer, *, iostat=io_stat) configs
            case ( 'DISORDER_WIDTH' )
                read(buffer, *, iostat=io_stat) sigma
            case default
                write(*,"(*(a))") 'Unable to assign input variable: ', label
            end select

        end if
    end do

    close(unit =  file_no)

    ! Normalise units to hw
    w00 = w00/hw
    lw = lw/hw
    ! hw = hw/hw
    JCoulx = JCoulx * eV / hw
    JCouly = JCouly * eV / hw
    JCoulz = JCoulz * eV / hw


    phi = pi*(phi/180.0_wp)
    mu_0 = mu_0/Debye
    k = w00/c
end subroutine







