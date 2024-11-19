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
    real(wp), allocatable :: xpl_osc(:,:)
    real(wp), allocatable :: ypl_osc(:,:)
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
subroutine construct_covariance_matrix
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




program test

    use variables
    use random_normal_distr
    implicit none

    integer(wp) N_samples 
    real(wp),allocatable :: MU(:)
    integer l,liter
    real(wp),allocatable :: SAMPLES(:,:)
    real(wp) :: new_time,old_time,Xmean
    real(wp),allocatable :: X(:)
    lattice_dimx = 5
    N_samples = lattice_dimx
    liter=100000
    allocate(MU(N_samples))
    allocate(SAMPLES(N_samples,liter))
    allocate(X(N_samples))
    MU=0.0_wp
    l0=0.0_wp
    call LatticeIndex()
    call construct_covariance_matrix()
    ! call draw_multivar_distr(N_samples,MU,sigma,A_covar)
    ! write(*,*) MU
    call cpu_time(old_time)
    do l=1,liter
        MU = 0.0_wp
        call draw_multivar_distr(N_samples,MU,sigma,A_covar)
        SAMPLES(:,l) = MU
    end do
    call cpu_time(new_time)
    write(*,*) 'Took',new_time-old_time,'s to generate',N_samples,'numbers',liter,'times'
    ! write(*,*) SAMPLES
    print "('central moments',/,a10,*(i10))","mean",1,2,3,4
    do l=1,N_samples
        X = SAMPLES(l,:)
        Xmean = sum(X)/liter
        X= X - Xmean
        print "(*(f10.4))",Xmean,sum(x**1)/liter,sum(x**2)/liter,sum(x**3)/liter,sum(x**4)/liter
    end do
    print*,"theoretical"
    print "(*(f10.4))",0.0_wp,0.0_wp,sigma**2,0.0_wp,3*sigma**4
end program