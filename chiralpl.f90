! chiralpl
! version 0
! 2024 Louis Minion
! This program uses the Frenkel-Holstein Hamiltonian to model a 3D lattice of organic chromophores
! It uses the multiparticle basis set of Philpott, truncating to one- and two-particle states only (two-particle approximation)
! The Hamiltonian is first constructed, then projected into its eigenbasis.
! Properties are then calculated from the eigenstates and vectors.
! Hartree units are used throughout.
module variables
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    
    !Constants
    real(wp), parameter :: eV = 8065.0_wp !1 eV = 8065cm^-1
    real(wp), parameter :: pi = 4.d0*datan(1.d0) ! good way of getting pi to highest machine precision available
    real(wp), parameter :: epsilon = 1.0_wp 
    real(wp), parameter :: Debye = 2.541765_wp ! Debye per dipole in AU
    real(wp), parameter :: c = 137.035999177 ! 1/alpha in Hartree atomic units
    real(wp), parameter :: kB = 0.6956925_wp ! boltzmann constant Eh/K
    ! INPUTS TO ENSURE A DEFAULT
    character*256 :: INPUT_NAME
    integer:: lattice_dimx = 1
    integer:: lattice_dimy = 1
    integer:: lattice_dimz = 1
    integer:: configs = 1
    logical :: bool_one_particle_states = .true.
    logical :: bool_two_particle_states = .true.
    logical :: H_out = .false.
    logical :: save_evals = .false.
    logical :: save_evecs = .false.
    logical :: manual_coupling = .true. ! if true, then doesn't calculate JCouplings individually, just uses inputs JCoulx,JCouly,JCoulz
    integer :: max_vibs
    integer :: n_nearest_neighbour = 1
    real(wp) :: lambda_neutral  = 1.0_wp
    real(wp) :: w00 = 14000.0_wp
    real(wp) :: hw = 1400.0_wp
    real(wp) :: mu_0 = 1.0_wp*Debye
    real(wp) :: temp = 0.0000001_wp
    real(wp) :: JCoulx = -0.2_wp
    real(wp) :: JCouly = -0.15_wp
    real(wp) :: JCoulz = 0.5_wp
    real(wp) :: te_x = 700.0_wp
    real(wp) :: te_y = 700.0_wp    
    real(wp) :: th_x = 700.0_wp
    real(wp) :: th_y = 700.0_wp  
    real(wp) :: lw  = 250.0_wp
    real(wp) :: phi = 0.0_wp ! twist angle for chiral aggregates
    real(wp) :: x_spacing = 1.0_wp
    real(wp) :: y_spacing = 1.0_wp
    real(wp) :: z_spacing = 1.0_wp ! in case stacking is further spaced than spacing along polymer
    real(wp) :: k = 1.0_wp ! k=w00/c
    real(wp) :: sigma=1
    real(wp),dimension(3) :: l0 = 0.0_wp
    ! disorder
    real(wp),allocatable :: A_covar(:,:)
    real(wp),allocatable :: diagonal_disorder_offsets(:)
    ! state counters
    integer :: general_counter = 0
    integer :: lattice_count = 0
    integer :: one_particle_counter = 0
    integer :: two_particle_counter = 0


    ! index arrays
    integer, allocatable :: lattice_index_arr(:,:,:)
    integer, allocatable :: one_particle_index_arr(:,:)
    integer, allocatable :: two_particle_index_arr(:,:,:,:)

    ! franck-condon table
    real(wp), allocatable :: fc_ground_to_neutral(:,:)

    ! array of dipole moment vectors
    real(wp), allocatable :: mu_xyz(:,:)

    ! hamiltonian, H
    ! The hamiltonian is a 2d array
    real(wp), allocatable :: H(:,:)
    real(wp), allocatable :: EVAL(:)
    integer, parameter :: empty = -1

    integer :: IU = 1
    integer :: EVAL_COUNT


    ! File out names
    character*256 :: eval_out_f


    ! For property calculations
    complex(kind=wp), parameter :: complex_zero = ( 0_wp, 0_wp )
    complex(kind=wp), allocatable:: abs_osc_strengths_x(:)
    complex(kind=wp), allocatable:: abs_osc_strengths_y(:)
    complex(kind=wp), allocatable:: abs_osc_strengths_x_configavg(:)
    complex(kind=wp), allocatable:: abs_osc_strengths_y_configavg(:)
    real(wp), allocatable :: xpl_osc(:,:)
    real(wp), allocatable :: ypl_osc(:,:)
    real(wp), allocatable :: xpl_osc_configavg(:,:)
    real(wp), allocatable :: ypl_osc_configavg(:,:)  
    integer, parameter  :: spec_steps = 2600
    real(wp), allocatable :: pl_specx(:)
    real(wp), allocatable :: pl_specy(:)
    real(wp), allocatable :: abs_specx(:)
    real(wp), allocatable :: abs_specy(:)   
    real(wp), dimension(2) :: rot_strengths
    real(wp), allocatable :: cpl_spec(:)

    contains

        function cross(a,b) result(crs)
            implicit none
            real(wp), dimension(3) :: crs
            real(wp), intent(in),dimension(3) :: a,b
            crs(1) = a(2)*b(3) - a(3)*b(2)
            crs(2) = a(3)*b(1) - a(1)*b(3)
            crs(3) = a(1)*b(2) - a(2)*b(1)
            return
        end function

end module

module random_normal_distr
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: box_muller, box_muller_vec
    contains
        function box_muller() result(Z)
            implicit none
            integer(wp) i
            real(wp) :: U(2)
            real(wp) :: Z(2)
            real(wp),parameter :: pi=4.0_wp*datan(1.0_wp)
            ! while loop to exclude log(0)
            do
                call RANDOM_NUMBER(U)
                if (U(1) > 0.0_wp) exit
            end do
            Z(1)=dsqrt(-2*log(U(1)))*cos(2.0_wp*pi*U(2))
            Z(2)=dsqrt(-2*log(U(1)))*cos(2.0_wp*pi*U(2))
        end function

        function box_muller_vec(N,mu,sigma) result(Z)
            implicit none
            integer(wp), intent(in) :: N
            real(wp), intent(in) :: mu, sigma
            real(wp) :: D(2)
            real(wp),allocatable :: Z(:)
            integer(wp) :: i,j

            allocate(Z(N))
            do i=1,N/2 ! two random numbers are returned for each
                j = 2*i-1
                Z(j:j+1) = box_muller()
            end do
            ! if N odd, then last element will be empty
            if (mod(N,2) .ne. 0) then
                D = box_muller()
                Z(N) = D(1)
            end if
            Z = Z*sigma
            Z = Z + mu
            end function
end module





program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    use omp_lib
    implicit none
    integer :: j, o
    real(wp) :: estimated_RAM,begin_disorder_avgtime,end_disorder_avgtime
    character*8  :: date_now
    character*10 :: time_now
    integer :: threads
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
    threads = omp_get_num_threads()
    print*, 'Running on',threads,'threads'
    call cpu_time(begin_disorder_avgtime)
    !$OMP PARALLEL PRIVATE(abs_osc_strengths_x,abs_osc_strengths_y,xpl_osc,ypl_osc) SHARED(abs_osc_strengths_x_configavg,abs_osc_strengths_y_configavg,xpl_osc_configavg,ypl_osc_configavg)
    !$OMP DO
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
    !$OMP END DO
    !$OMP END PARALLEL
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
    integer :: errstat
    logical :: exists
    integer :: io_stat, file_no, line_no, pos_space
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



subroutine LatticeIndex()
    use variables
    implicit none
    integer :: i_x, i_y, i_z

    ! Check if arr allocated, if not allocate as 3D array of size lattice_dimx,lattice_dimy,lattice_dimz
    if ( .not. allocated( lattice_index_arr ) ) then 
        allocate( lattice_index_arr( lattice_dimx, lattice_dimy, lattice_dimz) )
    end if
    lattice_index_arr = empty
    ! Write chromophore indices to index_arr elements
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                lattice_count = lattice_count + 1
                lattice_index_arr( i_x, i_y, i_z ) = lattice_count
            end do
        end do
    end do
print*, '*****************************************'
print*, lattice_count,'Lattice Sites'
end subroutine


subroutine oneParticleIndex()
    use variables
    implicit none
    integer :: i_x, i_y, i_z, vib, indx_xyz
    ! Check if arr allocated, if not then allocate as 2D array of size (no_sites_lattice, max_vibs+1)
    if ( .not. allocated( one_particle_index_arr ) ) then 
        allocate( one_particle_index_arr( lattice_count, 0:max_vibs) )
    end if
    one_particle_index_arr = empty
    ! iterate over all sites in 3d lattice, each of which is a row in the matrix, with columns as vib state inds
    do i_x = 1, lattice_dimx
        do i_y = 1, lattice_dimy
            do i_z = 1, lattice_dimz
                do vib = 0, max_vibs
                    general_counter = general_counter + 1
                    indx_xyz = lattice_index_arr(i_x, i_y, i_z)
                    one_particle_index_arr( indx_xyz, vib) = general_counter
                    one_particle_counter = one_particle_counter + 1
                end do
            end do 
        end do
    end do
    print*, one_particle_counter, 'One particle states'
end subroutine

subroutine twoParticleIndex()
    use variables
    ! index two-particle states; two particle states consist of a vibronic excitation on one site and a vibrational excitation on another
    ! therefore each combination of site with vibronic exc and vibrational exc on other site need an index number
    implicit none
    integer i_x, i_y, i_z, i_xyz, vib ! indices for vibronic excitations
    integer i_xv, i_yv, i_zv, i_xyzv, vibv ! indices for vibrational excitations


    if (.not. allocated(two_particle_index_arr)) then
        allocate( two_particle_index_arr(lattice_count, 0:max_vibs, lattice_count, 1:max_vibs))
    
    end if
    two_particle_index_arr = empty
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                do vib=0,max_vibs ! for each vibronic state, loop over all other lattice sites for vibrational excitations
                    do i_xv=1,lattice_dimx
                        do i_yv=1,lattice_dimy
                            do i_zv=1,lattice_dimz
                                do vibv=1,max_vibs ! excited vibrational states in ground electronic states cannot be v=0
                                    i_xyz= lattice_index_arr(i_x,i_y,i_z)
                                    i_xyzv=lattice_index_arr(i_xv,i_yv,i_zv)
                                    if ( i_xyz == i_xyzv ) cycle ! no states with vibrational and vibronic exc on same site
                                    if (vibv + vib > max_vibs) cycle 
                                    general_counter = general_counter + 1
                                    two_particle_counter = two_particle_counter + 1
                                    two_particle_index_arr(i_xyz, vib, i_xyzv, vibv) = general_counter
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    print*, two_particle_counter, 'Two particle states'
    print*, general_counter, 'Total basis states'
    print*, '*****************************************'                 
end subroutine


subroutine calcFranckCondonTables()
    ! precalculates vibrational overlap integrals into table
    ! starting with just ground-state to frenkel exciton type <m|n> factors
    use variables
    implicit none
    integer :: ground_vib, exc_vib
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
    integer , intent (in) :: n, m
    integer :: j
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
    integer, intent(in) :: j
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

subroutine dipole_moment
    use variables
    implicit none
    integer :: ix,iy,iz,ixyz
    allocate(mu_xyz(lattice_count,3))
    mu_xyz = 0.0_wp
    do ix=1,lattice_dimx
        do iy=1,lattice_dimy
            do iz=1,lattice_dimz
                ixyz = lattice_index_arr(ix,iy,iz)
                mu_xyz(ixyz,1) = mu_0*cos((iz*1.0_wp)*phi) ! H-stack of dipoles
                mu_xyz(ixyz,2) = mu_0*sin((iz*1.0_wp)*phi) ! H-stack of dipoles
            end do
        end do
    end do
end subroutine

!pure
pure real(wp) function coupling(x1,y1,z1,x2,y2,z2)
    use variables
    integer, intent(in) :: x1,y1,z1,x2,y2,z2
    integer :: d_x, d_y, d_z, i_xyz1, i_xyz2
    real(wp) :: dx,dy,dz
    real(wp), dimension(3) :: R, R_norm
    real(wp) :: R_mag
    real(wp) :: CT_couple_alongchain
    real(wp) :: dipole_dot_product
    real(wp) :: mu1_dotR
    real(wp) :: mu2_dotR
    coupling  = 0.0_wp

    d_x = abs(x2-x1)
    d_y = abs(y2-y1)
    d_z = abs(z2-z1)
    if (manual_coupling .eq. .true.) then
        go to 78
    end if
    CT_couple_alongchain = 0.0_wp
    dipole_dot_product = 0.0_wp
    mu1_dotR = 0.0_wp
    mu2_dotR = 0.0_wp

    if (x1 .eq. x2) then
        if (z1 .eq. z2) then
            if (y1 .eq. y2) then
                return ! don't calc coupling if same lattice site
            else
                ! maybe add CT contribution to coupling
                continue ! for now skip
            end if 
        end if
    end if
    ! Calculate coupling using point dipole approximation, Jc = (mu1.mu2 - 3(mu1.R_norm)(mu2.R_norm))/(4*pi*epsilon*R_mag)
    dx = (d_x*1.0_wp)*x_spacing
    dy = (d_y*1.0_wp)*y_spacing
    dz = (d_z*1.0_wp)*z_spacing
    R(1) = dx
    R(2) = dy
    R(3) = dz
    R_mag = sqrt((dx**2)+(dy**2)+(dz**2))
    R_norm = (1/R_mag)*R
    i_xyz1 = lattice_index_arr(x1,y1,z1)
    i_xyz2 = lattice_index_arr(x2,y2,z2)
    dipole_dot_product = mu_xyz(i_xyz1,1)*mu_xyz(i_xyz2,1) + mu_xyz(i_xyz1,2)*mu_xyz(i_xyz2,2) + mu_xyz(i_xyz1,3)*mu_xyz(i_xyz2,3)
    mu1_dotR = mu_xyz(i_xyz1,1)*R_norm(1) + mu_xyz(i_xyz1,2)*R_norm(2) + mu_xyz(i_xyz1,2)*R_norm(2)
    mu2_dotR = mu_xyz(i_xyz2,1)*R_norm(1) + mu_xyz(i_xyz2,2)*R_norm(2) + mu_xyz(i_xyz2,2)*R_norm(2)
    coupling = (1/(4*pi*epsilon*(R_mag**3)))*(dipole_dot_product-(3*mu1_dotR*mu2_dotR))
    return
    78 if (d_y .eq. 0 .and. d_z .eq. 0) then
        if (d_x .eq. 1) then
            coupling = JCoulx !+ te_x + th_x
        end if
        return

    end if

    if (d_x .eq. 0 .and. d_z .eq. 0) then
        if (d_y .eq. 1) then
            coupling = JCouly !+ te_x + th_x
        end if
        return
    end if
    if (d_x .eq. 0 .and. d_y .eq. 0) then
        if (d_z .eq. 1) then
            coupling = JCoulz !+ te_x + th_x
        end if
        return

    end if
    return

end function

subroutine build1particleHamiltonian()
    use variables
    implicit none
    real(wp), external :: coupling
    integer :: i_x1, i_y1, i_z1, vib_i1, i_xyz1, h_i ! x,y,z coords for vibronic exc on 1
    integer :: i_x2, i_y2, i_z2, vib_i2, i_xyz2, h_j ! x,y,z coords for vibronic exc on 2
    ! h_i and h_j are hamiltonian indices
    do i_x1 = 1, lattice_dimx
        do i_y1 = 1, lattice_dimy
            do i_z1=1,lattice_dimz
                do vib_i1=0,max_vibs
                    i_xyz1 = lattice_index_arr( i_x1, i_y1, i_z1 )
                    h_i = one_particle_index_arr( i_xyz1, vib_i1 )
                    if ( h_i == empty ) cycle
                    H(h_i, h_i) = vib_i1*1.0_wp + w00 + diagonal_disorder_offsets(i_xyz1)
                    do i_x2=1, lattice_dimx
                        do i_y2=1, lattice_dimy
                            do i_z2=1, lattice_dimz
                                do vib_i2 = 0, max_vibs
                                    i_xyz2 = lattice_index_arr(i_x2, i_y2, i_z2)
                                    h_j = one_particle_index_arr(i_xyz2, vib_i2)
                                    if ( h_j == empty ) cycle
                                    if (h_j .eq. h_i) cycle
                                    H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2)*fc_ground_to_neutral(0,vib_i1)*fc_ground_to_neutral(0,vib_i2)
                                    H(h_j, h_i) = H(h_i, h_j)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

end subroutine

subroutine build2particleHamiltonian()
    use variables
    implicit none
    real(wp), external :: coupling
    integer :: i_x1, i_y1, i_z1, vib_i1, i_xyz1, h_i
    integer :: i_x1v, i_y1v, i_z1v, vib_i1v, i_xyz1v ! indices for vibrationally only excited site
    integer :: i_x2, i_y2, i_z2, vib_i2, i_xyz2, h_j
    integer :: i_xyz2v, vib_i2v
    do i_x1 = 1, lattice_dimx
        do i_y1 = 1, lattice_dimy
            do i_z1=1,lattice_dimz
                do vib_i1=0,max_vibs
                    i_xyz1 = lattice_index_arr(i_x1,i_y1,i_z1)
                    do i_x1v=1, lattice_dimx
                        do i_y1v=1, lattice_dimy
                            do i_z1v=1, lattice_dimz
                                do vib_i1v=1, max_vibs
                                    i_xyz1v = lattice_index_arr(i_x1v,i_y1v,i_z1v)
                                    h_i = two_particle_index_arr(i_xyz1,vib_i1, i_xyz1v, vib_i1v)
                                    if (h_i .eq. empty) cycle
                                    H(h_i,h_i) = (vib_i1+vib_i1v)*1.0_wp + w00
                                    do i_x2=1,lattice_dimx
                                        do i_y2=1,lattice_dimy
                                            do i_z2=1,lattice_dimz
                                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                                do vib_i2=0,max_vibs
                                                    i_xyz2v = i_xyz1v
                                                    vib_i2v = vib_i1v
                                                    h_j = two_particle_index_arr(i_xyz2, vib_i2, i_xyz2v, vib_i2v )
                                                    if (h_j .eq. empty .or. h_j .eq. h_i) then
                                                        continue
                                                    else
                                                        H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2)* &
                                                        fc_ground_to_neutral(0, vib_i1)*fc_ground_to_neutral(0,vib_i2)
                                                        H(h_j,h_i) = H(h_i,h_j)
                                                    end if
                                                    if (i_xyz1v .eq. i_xyz2) then ! if vibrat exc in state 1 at same place as vibronic on state 2 then exchange type 
                                                        i_xyz2v = i_xyz1
                                                        do vib_i2v=1, max_vibs
                                                            h_j = two_particle_index_arr(i_xyz2, vib_i2, i_xyz2v, vib_i2v)
                                                            if ( h_j .eq. empty .or. h_j .eq. h_i ) cycle
                                                            H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2)*fc_ground_to_neutral(vib_i2v, vib_i1)*fc_ground_to_neutral(vib_i1v,vib_i2)
                                                            H(h_j, h_i) = H(h_i, h_j)
                                                        end do
                                                    end if
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
    end do
end subroutine


subroutine build1particle2particleHamiltonian()
    use variables
    implicit none
    real(wp), external :: coupling
    integer i_x1, i_y1, i_z1, i_xyz1, vib_i1, h_i ! coords of the vibronic exc in 1 in the one particle state
    integer i_x2, i_y2, i_z2, i_xyz2, vib_i2, h_j ! coords of the vibronic (electronic+vib) excitation in 2 in the two particle state
    integer i_x2v, i_y2v, i_z2v, i_xyz2v, vib_i2v ! coords of the pure vibrational excitation in 2 in the two particle state
    do i_x1=1,lattice_dimx
        do i_y1=1,lattice_dimy
            do i_z1=1,lattice_dimz
                do vib_i1=0,max_vibs
                    i_xyz1 = lattice_index_arr(i_x1,i_y1,i_z1) !get lattice index of state 1 1-particle exc
                    h_i = one_particle_index_arr(i_xyz1, vib_i1)
                    if (h_i .eq. empty) cycle
                    do i_x2=1,lattice_dimx
                        do i_y2=1,lattice_dimy
                            do i_z2=1,lattice_dimz
                                do vib_i2=0,max_vibs
                                    i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                    ! all matrix elements are zero except those where the pure vibrational excitation in 2 is in the same site as electronic in 1
                                    i_x2v = i_x1 
                                    i_y2v = i_y1
                                    i_z2v = i_z1
                                    i_xyz2v = lattice_index_arr( i_x2v, i_y2v, i_z2v )
                                    do vib_i2v=1,max_vibs
                                        h_j = two_particle_index_arr(i_xyz2,vib_i2,i_xyz2v,vib_i2v)
                                        if (h_j .eq. empty) cycle
                                        H(h_i,h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2)*fc_ground_to_neutral(vib_i2v,vib_i1)*fc_ground_to_neutral(0,vib_i2)
                                        H(h_j,h_i) = H(h_i,h_j)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine


subroutine Diagonalize(A,RANGE,N,W,M,I_U)
    use variables
    implicit none
    character*1, intent(in) :: RANGE ! RANGE is which eigenvalues are calculated, RANGE='A' is all, 'V' and 'I' allow selection
    integer, intent(in) :: N,I_U !N is the order of the matrix, IU controls largest eigenvalue returned if RANGE='I'. IU not ref'd if RANGE='V' or 'A'
    real(wp), intent(inout) :: A(n,n) ! A is the Hamiltonian matrix, on exit it is assigned to the eigenvecs
    real(wp), intent(out) :: W(n) ! Array of eigenvalues in ascending order
    integer, intent(out)  :: M ! number of eigenvalues found.
    character*1 :: JOBV, UPLO ! JOBV controls whether eigenvecs and vals are calc'd (='V'), or just eigenvals (='N'). UPLO controls whether on exit from DSYEVR A stores upper or lower triangular matrix
    parameter  (JOBV='V', UPLO='U')

    integer :: LDA ! array leading dimension
    real(wp) :: VL = -3.0_wp ! lower and upper of interval to search for eigenvalues if RANGE='V'. Not accessed
    real(wp) :: VU = 20.0_wp
    integer :: IL = 1 ! lower bound of IL->IU if only selected eigenvals requested.
    real(wp) :: ABSTOL ! absolute error tolerance for eigenvals, 

    real(wp), external :: dlamchm ! utility func to determine machine parameters for minimum error tolerance without overflow
    real(wp), allocatable  :: Z(:,:) ! eigenvectors of A in columns
    integer :: LDZ ! leading dimension of z
    integer :: ISUPPZ(2*max(1,N))
    real(wp), allocatable :: WORK(:)
    real(wp) :: WORK_DIM(1)
    integer :: LWORK
    integer, allocatable :: IWORK(:)
    integer :: IWORK_DIM(1)
    integer :: LIWORK
    integer :: INFO
    real(wp) :: DIAG_START, DIAG_END
    real(wp), external :: dlamch
    external :: dsyevr



    ABSTOL = dlamch('Safe minimum')
    LDA=N
    LDZ=N
    allocate( Z( N, N ) )
    ! First, query workspace with LWORK and LIWORK set to -1. This calcs the optimal size of the WORK and IWORK arrays
    LWORK = -1
    LIWORK = -1
    call dsyevr(JOBV,RANGE,UPLO,N,A,LDA,VL,VU,IL,I_U,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORK_DIM,LWORK,IWORK_DIM,LIWORK,INFO)
    LWORK = WORK_DIM(1)
    LIWORK = IWORK_DIM(1)
    ! allocate WORK and IWORK arrays
    allocate(WORK(LWORK))
    allocate(IWORK(LIWORK))
    ! print*, ' Begin Hamiltonian diagonalization'
    ! print*, '*****************************************'
    call cpu_time(DIAG_START)
    call dsyevr(JOBV,RANGE,UPLO,N,A,LDA,VL,VU,IL,I_U,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    call cpu_time(DIAG_END)
    ! print*, '        Done'
    ! print*, M, 'eigenvalues found'
    ! print*, ' Diagonalization done in',(DIAG_END-DIAG_START),'seconds'
    ! print*, '*****************************************'
    A(:,1:M) = Z(:,1:M) ! Assign A to Z, (less of Z if M != N)
    deallocate( Z, WORK, IWORK )
end subroutine

subroutine absorption()
    use variables
    implicit none
    integer i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer j
    integer :: h_i
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
    integer :: spec_point, j
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
    integer :: spec_point
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
    integer i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer i_x2,i_y2,i_z2,i_xyz2,vib2 ! indices for ground state vibrational excitation (for two particle states)
    integer :: h_i, j
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
    integer :: vt
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
    integer :: spec_point, vt, j
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
    integer :: spec_point
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
    integer i_x, i_y, i_z, n, vib ! indices for vibronic excitations
    integer i_x2,i_y2,i_z2,m,vib2, n1 ! indices for ground state vibrational excitation (for two particle states)
    integer i_x3,i_y3,i_z3,n2,vib3 ! indices for third iteration
    integer :: h_i, j
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
    integer :: spec_point, vt, j
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
    integer :: spec_point
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
    integer nx,ny,nz,n,mx,my,mz,m
    integer j
    integer :: h_i
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



! Select vector delta from probability distribution P of possible sets of offsets.
! delta has an entry for each chromophore site in the lattice.
! P(delta) = 1/((2pi)^(N/2)*sqrt(detA)) *exp(-Sum(1/2*A^-1_nm*delta_n*delta_m))
! Entries in A are given by sigma^2/2*exp(-(n-m)/l0) <- this accounts for spatial correlation. We can introduce 3d spatial correlation by replacing
! |n-m| with r_mag/l0 with l0 no longer dimensionless.
! First step; populate A
! Second step; cholesky decomposition to get a matrix B such that BB^T=A
! Third step; generate a vector of random (normally distributed) numbers Z
! Fourth step; delta = mu + B*Z, where mu is the mean vector of the prob distribution.
! Computational effort can be saved by only doing the first two steps once, rather than once per configuration.
subroutine construct_covariance_matrix()
    use variables
    implicit none
    integer :: nx,ny,nz,mx,my,mz, nxyz, mxyz
    integer :: lx,ly,lz
    real(wp) :: dx,dy,dz,argc

    allocate(A_covar(lattice_count,lattice_count))
    do nx = 1, lattice_dimx
        do ny = 1, lattice_dimy
            do nz = 1, lattice_dimz
                nxyz = lattice_index_arr(nx, ny, nz)
                do mx=1,lattice_dimx
                    do my=1,lattice_dimy
                        do mz=1,lattice_dimz
                            mxyz = lattice_index_arr(mx, my, mz)
                            lx = (nx-mx)
                            ly = (ny-my)
                            lz = (nz-mz)
                            dx = (lx*1.0_wp)*x_spacing
                            dy = (ly*1.0_wp)*y_spacing
                            dz = (lz*1.0_wp)*z_spacing
                            if (ANY(l0 .eq. 0.0_wp)) then
                                ! write(*,*) '0 l0'
                                argc = 0.0_wp
                                if (nxyz .eq. mxyz) then
                                    A_covar(nxyz,mxyz) = 1.0_wp
                                else 
                                    A_covar(nxyz,mxyz) = 0.0_wp
                                end if
                            else
                                argc = dsqrt((dx/l0(1))**2 + (dy/l0(2))**2 + (dz/l0(3))**2)
                                A_covar(nxyz,mxyz) = ((sigma**2)/2)*dexp(-1.0_wp*argc)
                            end if
                        end do
                    end do
                end do
            end do 
        end do
    end do
end subroutine

subroutine cholesky_decomp(ACOV,N)
    use variables
    implicit none
    external :: DPOTRF ! lapack routine to calculate the cholesky decomp of a positive semi-definite matrix
    external :: DPSTRF
    character*1 :: UPLO
    integer :: INFO
    integer :: LDA
    integer, intent(in) :: N
    real(wp), intent(inout) :: ACOV(n,n)
    real(wp) :: CHOLESK_START, CHOLESK_END

    UPLO = 'L'
    LDA = N
    call cpu_time(CHOLESK_START)
    call DPOTRF(UPLO,N,ACOV,LDA,INFO)
    call cpu_time(CHOLESK_END)
    print*, ' CHOLESKY decomp done in',(CHOLESK_END-CHOLESK_START),'seconds'
    if (INFO .ne. 0) then
        if (INFO > 0) then
            write(*,*) 'Matrix not positive definite, implement positive semi-definite'
        end if
    end if


end subroutine

subroutine draw_multivar_distr(N_samples,MU,SIGMA,A_FAC)
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use random_normal_distr, only: box_muller_vec
    implicit none
    external DGEMV
    real(wp) :: SIGMA
    integer(wp),intent(in) :: N_samples
    real(wp), allocatable :: DELTA(:)
    real(wp), intent(inout) :: MU(N_samples)! array of means (these will most often be the same (0) in my program)
    real(wp), intent(in) :: A_FAC(N_samples,N_samples)
    character*1 :: TRANS
    integer :: N,M,LDA,INCX,INCY
    real(wp) :: ALPHA,BETA
    allocate(DELTA(N_samples))
    ! write(*,*) MU(1)
    DELTA = box_muller_vec(N_samples,MU(1),SIGMA) ! this line assumes all means are the same. This is designed to work with this distribution but should be changed if using for general multivariate normal.
    ! write(*,*) DELTA
    TRANS='N'
    M=N_samples
    N=N_samples
    LDA=N_samples
    INCX=1
    INCY=1
    ALPHA=1.0_wp
    BETA=1.0_wp
    call DGEMV(TRANS,M,N,ALPHA,A_FAC,LDA,DELTA,INCX,BETA,MU,INCY)
    deallocate(DELTA)
end subroutine



