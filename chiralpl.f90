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
    integer :: max_vibs
    real(wp) :: lambda_neutral  = 1.0_wp
    real(wp) :: w00 = 10000.0_wp
    real(wp) :: hw = 1000.0_wp
    real(wp) :: JCoulx = 700.0_wp
    real(wp) :: JCouly = 700.0_wp
    real(wp) :: JCoulz = 700.0_wp
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
    integer, parameter :: empty = -1

end module

program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    !declare variables


    ! do parameter setting

    call readInput()
    call LatticeIndex()
    if ( bool_one_particle_states ) call oneParticleIndex()
    if ( bool_two_particle_states ) call twoParticleIndex()
    
    call calcFranckCondonTables
    ! write(7,*) fc_ground_to_neutral
    print*,'Hamiltonian size()','(',general_counter,general_counter,')'
    if ( .not. allocated(H)) then
         allocate(H(general_counter,general_counter))
    end if
    H = 0.0_wp

    call build1particleHamiltonian()

    ! call build2particleHamiltonian

    ! call build1particle2particleHamiltonian


    ! call Diagonalize

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
        STOP
    end if
    inquire( file=trim(fname), exist=exists)
    if ( .not. exists ) then
        write(*,*) 'INPUT FILE DOES NOT EXIST'
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

            case default
                write(*,"(*(a))") 'Unable to assign input variable: ', label
            end select

        end if
    end do

    close(unit =  file_no)

    ! Normalise units to hw
    hw = hw/hw
    w00 = w00/hw
end subroutine



subroutine LatticeIndex()
    use variables
    implicit none
    integer :: i_x, i_y, i_z

    ! Check if arr allocated, if not allocate as 3D array of size lattice_dimx,lattice_dimy,lattice_dimz
    if ( .not. allocated( lattice_index_arr ) ) then 
        allocate( lattice_index_arr( lattice_dimx, lattice_dimy, lattice_dimz) )
    end if
    ! Write chromophore indices to index_arr elements
    do i_x=1,lattice_dimx
        do i_y=1,lattice_dimy
            do i_z=1,lattice_dimz
                lattice_count = lattice_count + 1
                lattice_index_arr( i_x, i_y, i_z ) = lattice_count
            end do
        end do
    end do

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
    print*, '*****************************************'
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
        do i_y1 = 1, lattice_dimz
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
    write(21, *) H

end subroutine
