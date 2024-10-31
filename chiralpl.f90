! Program to calculate the PL/CPL spectra of chiral aggregate using the Holstein Hamiltonian
! 2D coupling will be incorporated via the methods of R Ghosh
! Problem is to build the Hamiltonian of a system representing a chiral stack of conjugated polymers
! Use the theory developed by F Spano to calculate TDs of each chromophore
! - Two particle basis set

! Starting by putting in subroutine calls as basic logic breakdown
! The majority of this code was heavily inspired by the code of N J Hestand and R Ghosh in exciton_1d and polaron_cmsf90
! No code is reproduced here, but algorithms are reproduced under the MIT license those codes were released with.
module variables
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    
    !Constants
    real(wp), parameter :: eV = 8065.0_wp !1 eV = 8065cm^-1
    real(wp), parameter :: pi = 4.d0*datan(1.d0) ! good way of getting pi to highest machine precision available
    real(wp), parameter :: epsilon = 1.0_wp 

    ! INPUTS TO ENSURE A DEFAULT
    character*256 :: INPUT_NAME
    integer:: lattice_dimx = 1
    integer:: lattice_dimy = 1
    integer:: lattice_dimz = 1
    logical :: bool_one_particle_states = .true.
    logical :: bool_two_particle_states = .true.
    logical :: H_out = .false.
    logical :: save_evals = .false.
    logical :: save_evecs = .false.
    logical :: twist_z = .false.
    logical :: manual_coupling = .true. ! if true, then doesn't calculate JCouplings individually, just uses inputs JCoulx,JCouly,JCoulz
    integer :: max_vibs
    integer :: n_nearest_neighbour = 1
    real(wp) :: lambda_neutral  = 1.0_wp
    real(wp) :: w00 = 14000.0_wp
    real(wp) :: hw = 1400.0_wp
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
    real(wp), allocatable :: pl_osc(:)
    real(wp), allocatable :: xpl_osc(:)
    real(wp), allocatable :: ypl_osc(:)
    integer, parameter  :: spec_steps = 2600
    real(wp), allocatable :: pl_spec(:)
    real(wp), allocatable :: pl_specx(:)
    real(wp), allocatable :: pl_specy(:) 
end module

program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    integer :: j
    real(wp) :: estimated_RAM
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
    EVAL = 0.0_wp
    H = 0.0_wp

    if ( bool_one_particle_states ) call build1particleHamiltonian()
    if ( bool_two_particle_states ) call build2particleHamiltonian()
    if ( bool_one_particle_states .and. bool_two_particle_states ) call build1particle2particleHamiltonian()
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
    call Diagonalize(H,'A',general_counter, EVAL, EVAL_COUNT, IU)

    ! Save eigenvalue logic
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
    print*,'Calculating PL oscillator strengths'
    write(*,'(A,F10.5,A)') '0K: EMISSION FROM LOWEST EXCITON WITH ENERGY:',EVAL(1),'cm^-1'
    call pl()
    call pl_output()
    write(*,*)'Calculating PL spectrum'
    call calc_pl_spec()
    call write_pl()
    ! call absorption

    ! ! Need to setup dipoles in chiral manner
    ! call photoluminescence

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
            case( 'TWIST_Z' )
                read(buffer, *, iostat=io_stat) twist_z
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
                mu_xyz(ixyz,1) = cos(iz*phi) ! H-stack of dipoles
                mu_xyz(ixyz,2) = sin(iz*phi) ! H-stack of dipoles
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
                    H(h_i, h_i) = vib_i1*1.0_wp + w00
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
    print*, ' Begin Hamiltonian diagonalization'
    print*, '*****************************************'
    call cpu_time(DIAG_START)
    call dsyevr(JOBV,RANGE,UPLO,N,A,LDA,VL,VU,IL,I_U,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    call cpu_time(DIAG_END)
    print*, '        Done'
    print*, M, 'eigenvalues found'
    print*, ' Diagonalization done in',(DIAG_END-DIAG_START),'seconds'
    print*, '*****************************************'
    A(:,1:M) = Z(:,1:M) ! Assign A to Z, (less of Z if M != N)
    deallocate( Z, WORK, IWORK )
end subroutine




subroutine pl()
    use variables
    implicit none
    integer i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer i_x2,i_y2,i_z2,i_xyz2,vib2 ! indices for ground state vibrational excitation (for two particle states)
    integer :: h_i
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

    allocate(pl_osc(max_vibs+1)) ! need to account for 0-0 as well
    pl_osc = 0.0_wp
    allocate(xpl_osc(max_vibs+1))
    xpl_osc = 0.0_wp
    allocate(ypl_osc(max_vibs+1))
    ypl_osc = 0.0_wp
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
                    c_nv = H(h_i,1) ! get coefficient of one-particle state in eigenbasis
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
    pl_osc(1) = I_00%re
    xpl_osc(1) = I_00x%re
    ypl_osc(1) = I_00y%re
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
                    c_nv = H(h_i,1)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(1,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(1,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(1,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,1) !
                                c_nvmv = H(h_i,1)
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
    pl_osc(2) = I_01%re
    xpl_osc(2) = I_01x%re
    ypl_osc(2) = I_01y%re
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
                    c_nv = H(h_i,1)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(2,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(2,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(2,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,2) !
                                c_nvmv = H(h_i,1)
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
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,1)
                                c_nvmv = H(h_i,1)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,1)
                                c_mvnv = H(H_i,1)
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
    pl_osc(3) = I_02%re
    xpl_osc(3) = I_02x%re
    ypl_osc(3) = I_02y%re
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
                    c_nv = H(h_i,1)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(3,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(3,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(3,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,3) !
                                c_nvmv = H(h_i,1)
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
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,1)
                                c_nvmv = H(h_i,1)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,2)
                                c_mvnv = H(H_i,1)
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
    pl_osc(4) = I_03%re
    xpl_osc(4) = I_03x%re
    ypl_osc(4) = I_03y%re
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
                    c_nv = H(h_i,1)
                    I_from_n = I_from_n + c_nv*fc_ground_to_neutral(4,vib) ! *dipole moment
                    I_from_nx = I_from_nx + c_nv*fc_ground_to_neutral(4,vib)*mu_xyz(n,1)
                    I_from_ny = I_from_ny + c_nv*fc_ground_to_neutral(4,vib)*mu_xyz(n,2)
                end do
                ! second terms aren't summed over all one-site vibrations
                do i_x2=1,lattice_dimx   ! sum over all over sites
                    do i_y2=1,lattice_dimy
                        do i_z2=1,lattice_dimz
                            do vib2=0,max_vibs
                                i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                                if (i_xyz2 .eq. n) cycle
                                h_i = two_particle_index_arr(i_xyz2,vib2,n,4) !
                                c_nvmv = H(h_i,1)
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
                            i_xyz2 = lattice_index_arr(i_x2,i_y2,i_z2)
                            I_from_n = complex_zero
                            I_from_nx = complex_zero
                            I_from_ny = complex_zero
                            do vib=0, max_vibs
                                h_i = two_particle_index_arr(n,vib,i_xyz2,2)
                                c_nvmv = H(h_i,1)
                                h_i = two_particle_index_arr(i_xyz2,vib,n,2)
                                c_mvnv = H(H_i,1)
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
    pl_osc(5) = I_04%re
    xpl_osc(5) = I_04x%re
    ypl_osc(5) = I_04y%re

end subroutine


subroutine pl_output()
    use variables
    implicit none
    integer :: vt
    character*256 :: peak
    print*, '*****************************************'
    print*, '         PL Oscillator Strengths'
    print*, '*****************************************'
    write(*,'(A5,A14,A14,A14)') 'PEAKS','UNPOL','X','Y'
    do vt=1,max_vibs+1
        write(peak,'(A4,I1)'), 'I_0-',(vt-1)
        write(*,'(A5,F14.7,F14.7,F14.7)'), peak, pl_osc(vt), xpl_osc(vt),ypl_osc(vt)
    end do

end subroutine


subroutine calc_pl_spec()
    use variables
    implicit none
    integer :: spec_point, vt
    real(wp) :: spectrum_start, spectrum_end, energy
    real(wp) :: lineshape, sum_w,sum_wx, sum_wy
    real(wp) :: exciton_energy
    allocate(pl_spec(spec_steps))
    allocate(pl_specx(spec_steps))
    allocate(pl_specy(spec_steps))
    pl_spec = 0.0_wp
    pl_specx = 0.0_wp
    pl_specy = 0.0_wp
    spectrum_start = 0.0_wp
    spectrum_end = 20000.0_wp/hw
    exciton_energy = EVAL(1)   
    do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        sum_w = 0.0_wp
        sum_wx = 0.0_wp
        sum_wy = 0.0_wp
        do vt=0,max_vibs
            lineshape = dexp(-(energy - exciton_energy + (vt*1.0_wp))**2/(2.0_wp*(lw**2)))/dsqrt(2.0_wp*lw**2*pi)
            sum_w = sum_w + lineshape*pl_osc((vt+1))*((exciton_energy-(vt*1.0_wp))**3)
            sum_wx = sum_wx + lineshape*xpl_osc((vt+1))*((exciton_energy-(vt*1.0_wp))**3)
            sum_wy = sum_wy + lineshape*ypl_osc((vt+1))*((exciton_energy-(vt*1.0_wp))**3)
        end do
        pl_spec(spec_point) = sum_w
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
    101 format(*(F14.7, :, ","))
    spectrum_start = 0.0_wp   
    spectrum_end = spectrum_start + 20000.0_wp/hw 
    open(unit=668, file=eval_out_f, action='write')
    102 format(*(A14, :, ","))
    write(668,*) trim(INPUT_NAME)
    write(668,*) 'ENERGY OF EMITTING EXCITON',EVAL(1)*hw
    write(668,*) 'hw',hw
    write(668,*) 'MAX_VIBS',max_vibs
    write(668,102) 'Energy','PL','PLX','PLY'
    do spec_point=1,spec_steps
        energy = spectrum_start+ ((spectrum_end-spectrum_start)*(1.0_wp*spec_point))/(spec_steps*1.0_wp)
        write( 668, 101 ) energy*hw ,pl_spec(spec_point), pl_specx(spec_point),pl_specy(spec_point)
    end do
    close(668)
end subroutine



subroutine cpl()
    use variables
    implicit none
    integer i_x, i_y, i_z,i_xyz, n, vib ! indices for vibronic excitations
    integer i_x2,i_y2,i_z2,i_xyz2,vib2 ! indices for ground state vibrational excitation (for two particle states)
    integer :: h_i
    real(wp) :: c_nv
    real(wp) :: c_nvmv, c_mvnv
    write(*,*) 'CPL calculation'
end subroutine