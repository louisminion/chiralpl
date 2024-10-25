! Program to calculate the PL/CPL spectra of chiral aggregate using the Holstein Hamiltonian
! 2D coupling will be incorporated via the methods of R Ghosh
! Problem is to build the Hamiltonian of a system representing a chiral stack of conjugated polymers
! Use the theory developed by F Spano to calculate TDs of each chromophore
! - Two particle basis set

! Starting by putting in subroutine calls as basic logic breakdown
! The majority of this code was heavily inspired by the code of N J Hestand and R Ghosh in exciton_1d and polaron_cmsf90
! No code is reproduced here, but algorithms are reproduced under the MIT license those codes were released with.
module variables
    ! INPUTS TO ENSURE A DEFAULT
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    character*256 :: INPUT_NAME
    integer:: lattice_dimx = 1
    integer:: lattice_dimy = 1
    integer:: lattice_dimz = 1
    logical :: bool_one_particle_states = .true.
    logical :: bool_two_particle_states = .true.
    logical :: H_out = .false.
    integer :: max_vibs
    real(wp) :: lambda_neutral  = 1.0_wp
    real(wp) :: w00 = 10000.0_wp
    real(wp) :: hw = 1400.0_wp
    real(wp) :: JCoulx = -0.2_wp
    real(wp) :: JCouly = -0.15_wp
    real(wp) :: JCoulz = 0.5_wp
    real(wp) :: te_x = 700.0_wp
    real(wp) :: te_y = 700.0_wp    
    real(wp) :: th_x = 700.0_wp
    real(wp) :: th_y = 700.0_wp    

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


    ! hamiltonian, H
    ! The hamiltonian is a 2d array
    real(wp), allocatable :: H(:,:)
    real(wp), allocatable :: EVAL(:)
    integer, parameter :: empty = -1

    integer :: IU = 1
    integer :: EVAL_COUNT

end module

program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    integer :: j
    real(wp) :: estimated_RAM
    external :: dsyevr, dlamch
    !declare variables


    ! do parameter setting

    call readInput()
    call LatticeIndex()
    if ( bool_one_particle_states ) call oneParticleIndex()
    if ( bool_two_particle_states ) call twoParticleIndex()
    
    call calcFranckCondonTables
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
        open(unit=10, file='h.csv')
        do j = 1, size(H,dim=1)
            write(10, '(*(F12.8 : ", "))') H(:, j)
        end do
        close(10)
    else
        write(*,*) 'Skipping saving Hamiltonian' 
    end if 
    call Diagonalize(H,'A',general_counter, EVAL, EVAL_COUNT, IU)

    write(8,*) EVAL

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
            case ('WRITE_HAMILTONIAN_OUT')
                read(buffer, *, iostat=io_stat) H_out
            case default
                write(*,"(*(a))") 'Unable to assign input variable: ', label
            end select

        end if
    end do

    close(unit =  file_no)

    ! Normalise units to hw
    w00 = w00/hw
    hw = hw/hw
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
        allocate(fc_ground_to_neutral(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j
    end if

    do ground_vib=0,max_vibs
        do exc_vib=0,max_vibs
            call calcVibrationalOverlap(0.0_wp,ground_vib,lambda_neutral,exc_vib,vibolap)
            fc_ground_to_neutral(ground_vib, exc_vib) = vibolap
        end do
    end do


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

        vibolap = vibolap + ((-1.0_wp)**(m-j)/(factorial(n-j)*factorial(j)*factorial(m-j)))*   &
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


pure real(wp) function coupling(x1,y1,z1,x2,y2,z2)
    use variables
    integer, intent(in) :: x1,y1,z1,x2,y2,z2
    integer :: d_x, d_y, d_z
    d_x = abs(x2-x1)
    d_y = abs(y2-y1)
    d_z = abs(z2-z1)
    coupling  = 0.0_wp
    if (d_y .eq. 0 .and. d_z .eq. 0) then
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
    integer :: i_x1, i_y1, i_z1, vib_i1, i_xyz1, h_i
    integer :: i_x2, i_y2, i_z2, vib_i2, i_xyz2, h_j
    ! h_i and h_j are hamiltonian indices
    do i_x1 = 1, lattice_dimx
        do i_y1 = 1, lattice_dimy
            do i_z1=1,lattice_dimz
                do vib_i1=0,max_vibs
                    i_xyz1 = lattice_index_arr( i_x1, i_y1, i_z1 )
                    h_i = one_particle_index_arr( i_xyz1, vib_i1 )
                    if ( h_i == empty ) cycle
                    H(h_i, h_i) = vib_i1*1.0_wp + w00
                    ! write(*,*) H(h_i, h_i) 
                    do i_x2=1, lattice_dimx
                        do i_y2=1, lattice_dimy
                            do i_z2=1, lattice_dimz
                                do vib_i2 = 0, max_vibs
                                    i_xyz2 = lattice_index_arr(i_x2, i_y2, i_z2)
                                    h_j = one_particle_index_arr(i_xyz2, vib_i2)
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
                                    ! write(*,*) H(h_i, h_i)
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
    ! write(21, *) H
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
                    ! write(*,*) 'h_i', h_i,'at i_xyz1', i_xyz1, 'vib_i1', vib_i1
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
                                        H(h_i,h_j) = coupling(i_x1,i_y2,i_z2,i_y1,i_y2,i_z2)*fc_ground_to_neutral(vib_i2v,vib_i1)*fc_ground_to_neutral(0,vib_i2)
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
    real(wp) :: VL = 10.0_wp ! lower and upper of interval to search for eigenvalues if RANGE='V'. Not accessed
    real(wp) :: VU = 40.0_wp
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