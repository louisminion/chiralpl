module disorder
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: construct_covariance_matrix, draw_multivar_distr, cholesky_decomp
    contains
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
            ! Defines the random distribution of diagonal disorder and how correlated close together sites are
            use variables
            implicit none
            integer(wp) :: nx,ny,nz,mx,my,mz, nxyz, mxyz
            integer(wp) :: lx,ly,lz
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
            ! Calls LAPACKs DPOTRF to cholesky decomp the covariance matrix A.
            use variables
            implicit none
            external :: DPOTRF ! lapack routine to calculate the cholesky decomp of a positive semi-definite matrix
            external :: DPSTRF
            character*1 :: UPLO
            integer(wp) :: INFO
            integer(wp) :: LDA
            integer(wp), intent(in) :: N
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
            ! Used to actually draw samples from the distribution defined by the covariance matrix ACOV
            ! N_samples is the length of the vector of random numbers.
            ! MU is returned as the vector of random numbers, on input it defines the mean of the distribution at each site.
            ! We define this to be zero (no-disorder), so MU is on entry an array of zeros in this program.
            ! SIGMA is the std deviation.
            ! A_FAC is the Cholesky decomposed B from the covariance matrix A_COV.
            ! We proceed by drawing a vector of size N_samples from a random normally distributed system (via Box Muller)
            ! Then we use LAPACK to change these into the random distribution we want via delta = mu + B*Z
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
            integer(wp) :: N,M,LDA,INCX,INCY
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

end module