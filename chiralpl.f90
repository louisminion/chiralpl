! chiralpl
! version 0
! 2024 Louis Minion
! This program uses the Frenkel-Holstein Hamiltonian to model a 3D lattice of organic chromophores
! It uses the multiparticle basis set of Philpott, truncating to one- and two-particle states only (two-particle approximation)
! The Hamiltonian is first constructed, then projected into its eigenbasis.
! Properties are then calculated from the eigenstates and vectors.
! Hartree units are used throughout.
program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    use index
    use disorder
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

subroutine calcFranckCondonTables()
    ! precalculates vibrational overlap integrals into table
    ! starting with just ground-state to frenkel exciton type <m|n> factors
    use variables
    implicit none
    integer(wp) :: ground_vib, exc_vib
    real(wp) :: vibolap
    if (.not. allocated(fc_ground_to_neutral)) then
        allocate(fc_ground_to_neutral(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
    end if

    do ground_vib=0,max_vibs
        do exc_vib=0,max_vibs
            call calcVibrationalOverlap(0.0_wp,ground_vib,lambda_neutral,exc_vib,vibolap)
            fc_ground_to_neutral(ground_vib, exc_vib) = vibolap
        end do
    end do
    !FCWRITE


end subroutine


subroutine calcVibrationalOverlap(lambda_1,n,lambda_2,m,vibolap)
    ! calc <m|n> via the recursion formula
    use variables
    implicit none
    real(kind=wp), external :: factorial
    integer(wp) , intent (in) :: n, m
    integer(wp) :: j
    real(kind=wp), intent(in):: lambda_1, lambda_2
    real(kind=wp) ::  vibolap
    real(kind=wp) :: lambda_hr

    lambda_hr = lambda_2 - lambda_1
    vibolap = 0.0_wp
    do j=0,min(n,m)

        vibolap = vibolap + (((-1.0_wp)**(m-j))/(factorial(n-j)*factorial(j)*factorial(m-j)))*   &
        (lambda_hr**(n+m-2*j))

    end do
    vibolap = vibolap*dsqrt(1.0_wp*factorial(n)*    &
    factorial(m))*dexp(-1.0_wp*lambda_hr**2/2.0_wp)
end subroutine


function factorial(j) result(fac)
    ! calc factorial(j) via gamma(j+1) = j!
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    integer(wp), intent(in) :: j
    real(kind=wp) :: fac
    if (j < 0) then
        write(*,*)  '<ERROR> Factorial undefined for arg', j, '<=0'
        write(*,*) ' <STOP> Bad argument to factorial'
        stop
    end if
    fac = gamma(real(j+1, kind=wp))

end function factorial

! function dipole_moment(x,y,z) result(mu)
!     use variables
!     implicit none
!     real(wp) , dimension(3) :: mu
!     integer, intent(in) :: x,y,z
!     mu = 0.0_wp
!     mu(1) = cos(z*phi)
!     mu(2) = sin(z*phi)
! end function


subroutine absorption()
    use variables
    implicit none
    integer(wp) i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer(wp) j
    integer(wp) :: h_i
    real(wp) :: c_nv

    abs_osc_strengths_x = complex_zero
    abs_osc_strengths_y = complex_zero

    do j=1,general_counter
        do i_x=1,lattice_dimx  ! sum over all n,v-tilde
            do i_y=1,lattice_dimy
                do i_z=1,lattice_dimz
                    i_xyz = lattice_index_arr(i_x,i_y,i_z) ! get n
                    do vib=0,max_vibs
                        h_i = one_particle_index_arr( i_xyz, vib )
                        if ( h_i == empty ) cycle
                        c_nv = H(h_i,j)
                        abs_osc_strengths_y(j) = abs_osc_strengths_y(j) + c_nv*fc_ground_to_neutral(0,vib)*mu_xyz(i_xyz,2)
                        abs_osc_strengths_x(j) = abs_osc_strengths_x(j) + c_nv*fc_ground_to_neutral(0,vib)*mu_xyz(i_xyz,1)
                    end do
                end do
            end do
        end do
        abs_osc_strengths_x(j) = abs_osc_strengths_x(j)*conjg(abs_osc_strengths_x(j))
        abs_osc_strengths_y(j) = abs_osc_strengths_y(j)*conjg(abs_osc_strengths_y(j))
    end do
end subroutine

subroutine calc_abs_spec()
    use variables
    implicit none
    integer(wp) :: spec_point, j
    real(wp) :: energy, eigenstate_energy
    real(wp) :: lineshape,sum_wx, sum_wy,step
    allocate(abs_specx(spec_steps))
    allocate(abs_specy(spec_steps))
    abs_specx = 0.0_wp
    abs_specy = 0.0_wp
    step = (min(maxval(EVAL),minval(EVAL) + 10000/hw) - minval(EVAL) + 8.0_wp*lw)/(1.0_wp*spec_steps) ! restrict spectrum to full or 10000cm-1 window
    energy = minval(eval)-4.0_wp*lw
    do spec_point=1,spec_steps
        energy = energy + step ! each spec_point add same amount to energy
        sum_wx = 0.0_wp
        sum_wy = 0.0_wp
        do j=1,general_counter
            eigenstate_energy = EVAL(j)
            lineshape = dexp(-(energy - eigenstate_energy)**2/(2.0_wp*(lw**2)))/dsqrt(2.0_wp*lw**2*pi)
            sum_wx = sum_wx + lineshape*abs_osc_strengths_x_configavg(j)
            sum_wy = sum_wy + lineshape*abs_osc_strengths_y_configavg(j)
        end do
        sum_wx = sum_wx*(1.0_wp/(1.0_wp*lattice_count))
        sum_wy = sum_wy*(1.0_wp/(1.0_wp*lattice_count)) ! 1/N normalisation
        abs_specx(spec_point) = sum_wx
        abs_specy(spec_point) = sum_wy
    end do

end subroutine


subroutine write_abs()
    use variables
    implicit none
    integer(wp) :: spec_point
    real(wp) :: step, energy

    write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_abs.csv')
    eval_out_f = trim(eval_out_f)
    write(*,('(a,a)')) 'Writing Absorption to ',eval_out_f
    108 format(*(F14.7, :, ","))
    open(unit=669, file=eval_out_f, action='write')
    109 format(*(A14, :, ","))
    write(669,*) trim(INPUT_NAME)
    write(669,*) 'hw',hw
    write(669,*) 'MAX_VIBS',max_vibs
    write(669,109) 'Energy','ABSX','ABSY'
    
    step = (min(maxval(EVAL),minval(EVAL) + 10000/hw) - minval(EVAL) + 8.0_wp*lw)/(1.0_wp*spec_steps) ! restrict spectrum to full or 10000cm-1 window
    energy = minval(eval)-4.0_wp*lw
    do spec_point=1,spec_steps
        energy = energy + step
        write( 669, 108 ) energy*hw, abs_specx(spec_point),abs_specy(spec_point)
    end do
    close(669)
end subroutine


subroutine pl()
    use variables
    implicit none
    integer(wp) i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer(wp) i_x2,i_y2,i_z2,i_xyz2,vib2 ! indices for ground state vibrational excitation (for two particle states)
    integer(wp) :: h_i, j
    real(wp) :: c_nv
    real(wp) :: c_nvmv, c_mvnv
    complex(kind=wp) :: I_from_n, I_twoparticle
    complex(kind=wp) :: I_from_nx, I_twoparticlex
    complex(kind=wp) :: I_from_ny, I_twoparticley
    complex(kind=wp) :: I_00
    complex(kind=wp) :: I_01
    complex(kind=wp) :: I_02
    complex(kind=wp) :: I_03
    complex(kind=wp) :: I_04

    complex(kind=wp) :: I_00x
    complex(kind=wp) :: I_01x
    complex(kind=wp) :: I_02x
    complex(kind=wp) :: I_03x
    complex(kind=wp) :: I_04x

    complex(kind=wp) :: I_00y
    complex(kind=wp) :: I_01y
    complex(kind=wp) :: I_02y
    complex(kind=wp) :: I_03y
    complex(kind=wp) :: I_04y

    xpl_osc = 0.0_wp
    ypl_osc = 0.0_wp

    do j=1,general_counter
    !!!!!! 0-0 intensity calculation !!!!!!
    I_00 = complex_zero
    I_00x = complex_zero
    I_00y = complex_zero
    do i_x=1,lattice_dimx  ! sum over all n,v-tilde
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                do vib=0,max_vibs
                    i_xyz = lattice_index_arr(i_x,i_y,i_z) ! get n
                    h_i = one_particle_index_arr( i_xyz, vib )
                    c_nv = H(h_i,j) ! get coefficient of one-particle state in eigenbasis
                    I_00 = I_00 + c_nv*fc_ground_to_neutral(vib,0) !*dipole_moment(x,y,z) !implement variable dipole later
                    I_00x = I_00x + c_nv*fc_ground_to_neutral(vib,0)*mu_xyz(i_xyz,1)
                    I_00y = I_00y + c_nv*fc_ground_to_neutral(vib,0)*mu_xyz(i_xyz,2)
                end do
            end do
        end do
    end do
    I_00 = I_00*conjg(I_00) ! dconjg is obselete, conjg will generally know the kind of its variable
    I_00x = I_00x*conjg(I_00x)
    I_00y = I_00y*conjg(I_00y)
    xpl_osc(1,j) = I_00x%re
    ypl_osc(1,j) = I_00y%re
    if (max_vibs < 1) then
        return
    end if
    !!!!!! 0-1 intensity calculation !!!!!!
    I_01 = complex_zero
    I_01x = complex_zero
    I_01y = complex_zero
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                I_from_n = complex_zero
                I_from_nx = complex_zero
                I_from_ny = complex_zero
                do vib=0,max_vibs
                    h_i = one_particle_index_arr(n,vib)
                    c_nv = H(h_i,j)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(1,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(1,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(1,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                if (bool_two_particle_states .eq. .false.) cycle
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,1) !
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(0,vib2)
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,2)
                            end do
                        end do
                    end do
                end do
                I_from_n = I_from_n*conjg(I_from_n)
                I_from_nx = I_from_nx*conjg(I_from_nx)
                I_from_ny = I_from_ny*conjg(I_from_ny)
                I_01 = I_01 + I_from_n
                I_01x = I_01x + I_from_nx
                I_01y = I_01y + I_from_ny
            end do
        end do
    end do
    xpl_osc(2,j) = I_01x%re
    ypl_osc(2,j) = I_01y%re
    if (max_vibs < 2) then
        return
    end if
    !!!!!! 0-2 intensity calculation !!!!!!
    I_02 = complex_zero
    I_02x = complex_zero
    I_02y = complex_zero
    ! one-particle terms
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                I_from_n = complex_zero
                I_from_nx = complex_zero
                I_from_ny = complex_zero
                do vib=0,max_vibs
                    h_i = one_particle_index_arr(n,vib)
                    c_nv = H(h_i,j)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(2,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(2,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(2,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                if (bool_two_particle_states .eq. .false.) cycle
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,2) !
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(0,vib2)
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,2)
                            end do
                        end do
                    end do
                end do
                I_from_n = I_from_n*conjg(I_from_n)
                I_from_nx = I_from_nx*conjg(I_from_nx)
                I_from_ny = I_from_ny*conjg(I_from_ny)
                I_02 = I_02 + I_from_n
                I_02x = I_02x + I_from_nx
                I_02y = I_02y + I_from_ny
            end do
        end do
    end do
    ! two particle terms
    I_twoparticle = complex_zero
    I_twoparticlex = complex_zero
    I_twoparticley = complex_zero
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            if (bool_two_particle_states .eq. .false.) cycle
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,1)
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,1)
                                c_mvnv = H(H_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(1,vib) + c_mvnv*fc_ground_to_neutral(1,vib) !*dipole moments of n
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(1,vib)*mu_xyz(n,1) + c_mvnv*fc_ground_to_neutral(1,vib)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(1,vib)*mu_xyz(n,2) + c_mvnv*fc_ground_to_neutral(1,vib)*mu_xyz(i_xyz2,2)
                            end do
                            I_from_n  =  I_from_n*conjg(I_from_n)
                            I_from_nx  =  I_from_nx*conjg(I_from_nx)
                            I_from_ny  =  I_from_ny*conjg(I_from_ny)
                            I_twoparticle = I_twoparticle + I_from_n
                            I_twoparticlex = I_twoparticlex + I_from_nx
                            I_twoparticley = I_twoparticley + I_from_ny
                        end do
                    end do
                end do
            end do
        end do
    end do
    I_twoparticle = 0.5_wp*I_twoparticle
    I_twoparticlex = 0.5_wp*I_twoparticlex
    I_twoparticley = 0.5_wp*I_twoparticley
    I_02 = I_02 + I_twoparticle
    I_02x = I_02x + I_twoparticlex
    I_02y = I_02y + I_twoparticley
    xpl_osc(3,j) = I_02x%re
    ypl_osc(3,j) = I_02y%re
    ! write(*,*) I_02
    if (max_vibs < 3) then
        return
    end if
    !!!!!! 0-3 intensity calculation !!!!!!

    I_03 = complex_zero
    I_03x = complex_zero
    I_03y = complex_zero
    ! one-particle terms
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                I_from_n = complex_zero
                I_from_nx = complex_zero
                I_from_ny = complex_zero
                do vib=0,max_vibs
                    h_i = one_particle_index_arr(n,vib)
                    c_nv = H(h_i,j)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(3,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(3,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(3,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                if (bool_two_particle_states .eq. .false.) cycle
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,3) !
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(0,vib2)
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,2)
                            end do
                        end do
                    end do
                end do
                I_from_n = I_from_n*conjg(I_from_n)
                I_from_nx = I_from_nx*conjg(I_from_nx)
                I_from_ny = I_from_ny*conjg(I_from_ny)
                I_03 = I_03 + I_from_n
                I_03x = I_03x + I_from_nx
                I_03y = I_03y + I_from_ny
            end do
        end do
    end do
    ! two particle terms
    I_twoparticle = complex_zero
    I_twoparticlex = complex_zero
    I_twoparticley = complex_zero
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            if (bool_two_particle_states .eq. .false.) cycle
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,1)
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,2)
                                if (h_i .eq. empty) cycle
                                c_mvnv = H(H_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(1,vib) + c_mvnv*fc_ground_to_neutral(2,vib) !*dipole moments of n
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(1,vib)*mu_xyz(n,1) + c_mvnv*fc_ground_to_neutral(2,vib)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(1,vib)*mu_xyz(n,2) + c_mvnv*fc_ground_to_neutral(2,vib)*mu_xyz(i_xyz2,2)
                            end do
                            I_from_n  =  I_from_n*conjg(I_from_n)
                            I_from_nx  =  I_from_nx*conjg(I_from_nx)
                            I_from_ny  =  I_from_ny*conjg(I_from_ny)
                            I_twoparticle = I_twoparticle + I_from_n
                            I_twoparticlex = I_twoparticlex + I_from_nx
                            I_twoparticley = I_twoparticley + I_from_ny
                        end do
                    end do
                end do
            end do
        end do
    end do
    I_twoparticle = 0.5_wp*I_twoparticle
    I_twoparticlex = 0.5_wp*I_twoparticlex
    I_twoparticley = 0.5_wp*I_twoparticley
    I_03 = I_03 + I_twoparticle
    I_03x = I_03x + I_twoparticlex
    I_03y = I_03y + I_twoparticley
    xpl_osc(4,j) = I_03x%re
    ypl_osc(4,j) = I_03y%re
    if (max_vibs < 4) then
        return
    end if

    !!!!!! 0-4 intensity calculation !!!!!!

    I_04 = complex_zero
    I_04x = complex_zero
    I_04y = complex_zero
    ! one-particle terms
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                I_from_n = complex_zero
                I_from_nx = complex_zero
                I_from_ny = complex_zero
                do vib=0,max_vibs
                    h_i = one_particle_index_arr(n,vib)
                    c_nv = H(h_i,j)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(4,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(4,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(4,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                if (bool_two_particle_states .eq. .false.) cycle
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,4) !
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(0,vib2)
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(0,vib2)*mu_xyz(i_xyz2,2)
                            end do
                        end do
                    end do
                end do
                I_from_n = I_from_n*conjg(I_from_n)
                I_from_nx = I_from_nx*conjg(I_from_nx)
                I_from_ny = I_from_ny*conjg(I_from_ny)
                I_04 = I_04 + I_from_n
                I_04x = I_04x + I_from_nx
                I_04y = I_04y + I_from_ny
            end do
        end do
    end do
    ! two particle terms
    I_twoparticle = complex_zero
    I_twoparticlex = complex_zero
    I_twoparticley = complex_zero
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            if (bool_two_particle_states .eq. .false.) cycle
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,2)
                                if (h_i .eq. empty) cycle
                                c_nvmv = H(h_i,j)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,2)
                                if (h_i .eq. empty) cycle
                                c_mvnv = H(H_i,j)
                                I_from_n = I_from_n + c_nvmv*fc_ground_to_neutral(2,vib) + c_mvnv*fc_ground_to_neutral(2,vib) !*dipole moments of n
                                I_from_nx = I_from_nx + c_nvmv*fc_ground_to_neutral(2,vib)*mu_xyz(n,1) + c_mvnv*fc_ground_to_neutral(2,vib)*mu_xyz(i_xyz2,1)
                                I_from_ny = I_from_ny + c_nvmv*fc_ground_to_neutral(2,vib)*mu_xyz(n,2) + c_mvnv*fc_ground_to_neutral(2,vib)*mu_xyz(i_xyz2,2)
                            end do
                            I_from_n  =  I_from_n*conjg(I_from_n)
                            I_from_nx  =  I_from_nx*conjg(I_from_nx)
                            I_from_ny  =  I_from_ny*conjg(I_from_ny)
                            I_twoparticle = I_twoparticle + I_from_n
                            I_twoparticlex = I_twoparticlex + I_from_nx
                            I_twoparticley = I_twoparticley + I_from_ny
                        end do
                    end do
                end do
            end do
        end do
    end do
    I_twoparticle = 0.5_wp*I_twoparticle
    I_twoparticlex = 0.5_wp*I_twoparticlex
    I_twoparticley = 0.5_wp*I_twoparticley
    I_04 = I_04 + I_twoparticle
    I_04x = I_04x + I_twoparticlex
    I_04y = I_04y + I_twoparticley
    xpl_osc(5,j) = I_04x%re
    ypl_osc(5,j) = I_04y%re
    end do

end subroutine


subroutine pl_output()
    use variables
    implicit none
    integer(wp) :: vt
    character*256 :: peak
    print*, '*****************************************'
    print*, '         PL Oscillator Strengths'
    print*, '     (For the lowest energy exciton)'
    print*, '          Configurational average   '
    print*, '*****************************************'
    write(*,'(A5,A14,A14)') 'PEAKS','X','Y'
    do vt=1,max_vibs+1
        write(peak,'(A4,I1)'), 'I_0-',(vt-1)
        write(*,'(A5,F14.7,F14.7)'), peak, xpl_osc_configavg(vt,1),ypl_osc_configavg(vt,1)
    end do

end subroutine


subroutine calc_pl_spec()
    use variables
    implicit none
    integer(wp) :: spec_point, vt, j
    real(wp) :: spectrum_start, spectrum_end, energy
    real(wp) :: lineshape, sum_w,sum_wx, sum_wy
    real(wp) :: exciton_energy, boltzfactor, boltzsum
    real(wp), allocatable :: occupation_factors(:)
    allocate(pl_specx(spec_steps))
    allocate(pl_specy(spec_steps))
    allocate(occupation_factors(general_counter))
    pl_specx = 0.0_wp
    pl_specy = 0.0_wp
    spectrum_start = 0.0_wp
    spectrum_end = w00+20.0_wp*lw
    boltzsum = 0.0_wp
    if (temp .eq. 0.0_wp) then
        go to 967
    end if
    do j=1,general_counter
        exciton_energy = EVAL(j)
        boltzfactor = dexp((-1.0_wp*(exciton_energy-EVAL(1)))/(kB*temp))
        boltzsum = boltzsum + boltzfactor
        occupation_factors(j) = boltzfactor
    end do
    occupation_factors = occupation_factors/boltzsum
    967 do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        sum_w = 0.0_wp
        sum_wx = 0.0_wp
        sum_wy = 0.0_wp
        do vt=0,max_vibs
            boltzsum = 0.0_wp
            do j=1, general_counter
                exciton_energy = EVAL(j)
                if (temp .ne. 0.0_wp) then
                    if (occupation_factors(j) < 0.0001) cycle
                    lineshape = dexp(-(energy - exciton_energy + (vt*1.0_wp))**2/(2.0_wp*(lw**2)))/dsqrt(2.0_wp*lw**2*pi)
                    sum_wx = sum_wx + lineshape*xpl_osc_configavg((vt+1),j)*((exciton_energy-(vt*1.0_wp))**3)*occupation_factors(j)
                    sum_wy = sum_wy + lineshape*ypl_osc_configavg((vt+1),j)*((exciton_energy-(vt*1.0_wp))**3)*occupation_factors(j)
                else 
                    if (j .ne. 1) cycle
                    lineshape = dexp(-(energy - exciton_energy + (vt*1.0_wp))**2/(2.0_wp*(lw**2)))/dsqrt(2.0_wp*lw**2*pi)
                    sum_wx = sum_wx + lineshape*xpl_osc_configavg((vt+1),j)*((exciton_energy-(vt*1.0_wp))**3)
                    sum_wy = sum_wy + lineshape*ypl_osc_configavg((vt+1),j)*((exciton_energy-(vt*1.0_wp))**3)
                end if

            end do
            ! write(*,*) 'BOLTZMANN SUM', boltzsum, vt
            ! sum_wx = sum_wx*(1/boltzsum)
            ! sum_wy = sum_wy*(1/boltzsum)
        end do
        pl_specx(spec_point) = sum_wx
        pl_specy(spec_point) = sum_wy
    end do

end subroutine


subroutine write_pl()
    use variables
    implicit none
    integer(wp) :: spec_point
    real(wp) :: spectrum_start, spectrum_end, energy
    write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_pl.csv')
    eval_out_f = trim(eval_out_f)
    write(*,('(a,a)')) 'Writing PL to ',eval_out_f
    101 format(*(F20.13, :, ","))
    spectrum_start = 0.0_wp
    spectrum_end = w00+20.0_wp*lw
    open(unit=668, file=eval_out_f, action='write')
    102 format(*(A20, :, ","))
    write(668,*) trim(INPUT_NAME)
    write(668,*) 'ENERGY OF EMITTING EXCITON',EVAL(1)*hw
    write(668,*) 'hw',hw
    write(668,*) 'MAX_VIBS',max_vibs
    write(668,102) 'Energy','PLX','PLY'
    do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        write( 668, 101 ) energy*hw, pl_specx(spec_point),pl_specy(spec_point)
    end do
    close(668)
end subroutine



subroutine cpl()
    use variables
    implicit none
    integer(wp) i_x, i_y, i_z, n, vib ! indices for vibronic excitations
    integer(wp) i_x2,i_y2,i_z2,m,vib2, n1 ! indices for ground state vibrational excitation (for two particle states)
    integer(wp) i_x3,i_y3,i_z3,n2,vib3 ! indices for third iteration
    integer(wp) :: h_i, j
    real(wp),dimension(3) :: mu_n, mu_m ! dipole moment vectors for m and n 
    real(wp),dimension(3) :: r_m, r_n, rdiff ! position vectors for m and n, and between them
    real(wp) :: c_nv, c_mv
    real(wp) :: c_nvmv, c_mvnv
    real(wp) :: R_00, R_01
    real(wp), dimension(3) :: crossproduct
    write(*,*) 'CPL calculation'
    r_m = 0.0_wp
    r_n = 0.0_wp
    rdiff = 0.0_wp
    R_00 = 0.0_wp
    do i_x=1,lattice_dimx  ! sum over all n,v-tilde
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                r_n(1) = i_x*x_spacing
                r_n(2) = i_y*y_spacing
                r_n(3) = i_z*z_spacing
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            m = lattice_index_arr(i_x2,i_y2,i_z2)
                            if (n .eq. m) cycle
                            !<psi(em)|mu(n)|g,0> X <psi(em)|mu(m)|g,0>.(rm-rn)
                            mu_m = mu_xyz(m,:)
                            mu_n = mu_xyz(n,:)
                            crossproduct = cross(mu_n,mu_m)
                            r_m(1) = i_x2*x_spacing
                            r_m(2) = i_y2*y_spacing
                            r_m(3) = i_z2*z_spacing
                            rdiff = r_n-r_m
                            do vib=0,max_vibs
                                do vib2=0,max_vibs
                                    h_i = one_particle_index_arr( n, vib )
                                    c_nv = H(h_i,1)
                                    h_i = one_particle_index_arr( m, vib2 )
                                    c_mv = H(h_i,1)
                                    R_00 = R_00 + c_nv*c_mv*fc_ground_to_neutral(0,vib)*fc_ground_to_neutral(vib2,0)*dot_product(crossproduct,rdiff)  
                                end do

                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    R_00 = R_00*(k)/(mu_0**2)
    write(*,*) 'R_00', R_00
    rot_strengths(1) = R_00
    R_01 = 0.0_wp
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                r_n(1) = i_x*x_spacing
                r_n(2) = i_y*y_spacing
                r_n(3) = i_z*z_spacing
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            m = lattice_index_arr(i_x2,i_y2,i_z2)
                            if (n .eq. m) cycle
                            mu_m = mu_xyz(m,:)
                            mu_n = mu_xyz(n,:)
                            crossproduct = cross(mu_n,mu_m)
                            r_m(1) = i_x2*x_spacing
                            r_m(2) = i_y2*y_spacing
                            r_m(3) = i_z2*z_spacing
                            rdiff = r_n-r_m
                            do vib=0,max_vibs
                                do vib2=0,max_vibs
                                    h_i = one_particle_index_arr( n, vib )
                                    c_nv = H(h_i,1)
                                    h_i = two_particle_index_arr( m, vib2 ,n,1)
                                    if (h_i .eq. empty) cycle
                                    c_nvmv = H(h_i,1)
                                    R_01 = R_01 + fc_ground_to_neutral(vib,1)*fc_ground_to_neutral(0,vib2)*c_nv*c_nvmv*dot_product(crossproduct,rdiff)
                                end do
                            end do

                        end do
                    end do
                end do
            end do
        end do
    end do
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                n = lattice_index_arr(i_x,i_y,i_z)
                do i_x2=1,lattice_dimx
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            n1 = lattice_index_arr(i_x2,i_y2,i_z2)
                            r_n(1) = i_x2*x_spacing
                            r_n(2) = i_y2*y_spacing
                            r_n(3) = i_z2*z_spacing
                            do i_x3=1,lattice_dimx
                                do i_y3=1,lattice_dimy
                                    do i_z3=1,lattice_dimz
                                        n2 = lattice_index_arr(i_x3,i_y3,i_z3)
                                        if (n1 .eq. n2) cycle
                                        mu_m = mu_xyz(n2,:)
                                        mu_n = mu_xyz(n1,:)
                                        crossproduct = cross(mu_n,mu_m)
                                        r_m(1) = i_x3*x_spacing
                                        r_m(2) = i_y3*y_spacing
                                        r_m(3) = i_z3*z_spacing
                                        rdiff = r_n-r_m
                                        do vib2=0,max_vibs
                                            do vib3=0,max_vibs
                                                h_i = two_particle_index_arr(n1,vib2,n,1)
                                                if (h_i .eq. empty) cycle
                                                c_nv = H(h_i,1)
                                                h_i = two_particle_index_arr(n2,vib3,n,1)
                                                if (h_i .eq. empty) cycle
                                                c_nvmv = H(h_i,1)
                                                R_01 = R_01 + fc_ground_to_neutral(vib2,0)*fc_ground_to_neutral(vib3,0)*c_nv*c_nvmv*dot_product(crossproduct,rdiff)
                                            end do
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    R_01 = R_01*(k)/(mu_0**2)

    write(*,*) 'R_01', R_01
    rot_strengths(2) = R_01


    ! Add 0-2,0-3,0-4 peaks

    ! do j=1,2
    !     (1/(mu_0**2))*(xpl_osc(j,1) + ypl_osc(j,1))*rot_strengths(j) ! pl osc strength
    ! end do

end subroutine

subroutine calc_cpl_spec()
    use variables
    implicit none
    integer(wp) :: spec_point, vt, j
    real(wp) :: spectrum_start, spectrum_end, energy
    real(wp) :: lineshape, sum_I, sum_R
    real(wp) :: exciton_energy, boltzfactor, boltzsum
    real,dimension(spec_steps) :: smI, smR
    allocate(cpl_spec(spec_steps))
    cpl_spec = 0.0_wp
    spectrum_start = 0.0_wp
    spectrum_end = w00+20.0_wp*lw
    exciton_energy = EVAL(1)

    do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        sum_I = 0.0_wp
        sum_R = 0.0_wp
        do vt=0,1
            ! (1/(mu_0**2))*(xpl_osc(j,1) + ypl_osc(j,1))
            lineshape = dexp(-(energy - exciton_energy + (vt*1.0_wp))**2/(2.0_wp*(lw**2)))/dsqrt(2.0_wp*lw**2*pi)
            sum_I = sum_I + (lineshape*(xpl_osc((vt+1),1)+ypl_osc((vt+1),1))*((exciton_energy-(vt*1.0_wp))**3))*(1/(mu_0**2))
            sum_R = sum_R + (lineshape*(rot_strengths((vt+1)))*((exciton_energy-(vt*1.0_wp))**3))
        end do
        smR(spec_point) = sum_R
        smI(spec_point) = sum_I
        if (sum_I .eq. 0.0_wp) then
            cpl_spec(spec_point) = 0.0_wp
        else
            cpl_spec(spec_point) = sum_R/sum_I
        end if
    end do
    write(34,*) smI
    write(56,*) smR

end subroutine

subroutine write_cpl()
    use variables
    implicit none
    integer(wp) :: spec_point
    real(wp) :: spectrum_start, spectrum_end, energy
    write(eval_out_f,'(a,a)') trim(INPUT_NAME), trim('_cpl.csv')
    eval_out_f = trim(eval_out_f)
    write(*,('(a,a)')) 'Writing CPL to ',eval_out_f
    105 format(*(F20.13, :, ","))
    spectrum_start = 0.0_wp
    spectrum_end = w00+20.0_wp*lw
    open(unit=670, file=eval_out_f, action='write')
    106 format(*(A20, :, ","))
    write(670,*) trim(INPUT_NAME)
    write(670,*) 'ENERGY OF EMITTING EXCITON',EVAL(1)*hw
    write(670,*) 'hw',hw
    write(670,*) 'MAX_VIBS',max_vibs
    write(670,106) 'Energy','GLUM'
    do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        write( 670, 105 ) energy*hw, cpl_spec(spec_point)
    end do
    close(670)
end subroutine

subroutine cd()
    use variables
    implicit none
    ! integer i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer(wp) nx,ny,nz,n,mx,my,mz,m
    integer(wp) j
    integer(wp) :: h_i
    real(wp) :: c_nv
    real(wp), dimension(3) :: mu_m, mu_n
    real(wp), dimension(3) :: crossproduct    

    do j=1,general_counter
        do nx=1,lattice_dimx
            do ny=1,lattice_dimy
                do nz=1,lattice_dimz
                    n = lattice_index_arr(nx,ny,nz)
                    mu_n = mu_xyz(n,:)
                    do mx=1,lattice_dimx
                        do my=1,lattice_dimy
                            do mz=1,lattice_dimy
                                m = lattice_index_arr(mx,my,mz)
                                mu_m = mu_xyz(m,:)
                                crossproduct = cross(mu_n,mu_m)

                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine



