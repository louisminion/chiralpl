module hamiltonian
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    private
    public :: dipole_moment, coupling, build1particleHamiltonian, build2particleHamiltonian, build1particle2particleHamiltonian, Diagonalize

    contains
        subroutine dipole_moment
            use variables
            implicit none
            integer(wp) :: ix,iy,iz,ixyz
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

        real(wp) function coupling(x1,y1,z1,x2,y2,z2)
            integer(wp), intent(in) :: x1,y1,z1,x2,y2,z2
            integer(wp) :: d_x, d_y, d_z, i_xyz1, i_xyz2
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
            integer(wp) :: i_x1, i_y1, i_z1, vib_i1, i_xyz1, h_i ! x,y,z coords for vibronic exc on 1
            integer(wp) :: i_x2, i_y2, i_z2, vib_i2, i_xyz2, h_j ! x,y,z coords for vibronic exc on 2
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
            integer(wp) :: i_x1, i_y1, i_z1, vib_i1, i_xyz1, h_i
            integer(wp) :: i_x1v, i_y1v, i_z1v, vib_i1v, i_xyz1v ! indices for vibrationally only excited site
            integer(wp) :: i_x2, i_y2, i_z2, vib_i2, i_xyz2, h_j
            integer(wp) :: i_xyz2v, vib_i2v
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
            integer(wp) i_x1, i_y1, i_z1, i_xyz1, vib_i1, h_i ! coords of the vibronic exc in 1 in the one particle state
            integer(wp) i_x2, i_y2, i_z2, i_xyz2, vib_i2, h_j ! coords of the vibronic (electronic+vib) excitation in 2 in the two particle state
            integer(wp) i_x2v, i_y2v, i_z2v, i_xyz2v, vib_i2v ! coords of the pure vibrational excitation in 2 in the two particle state
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
            integer(wp), intent(in) :: N,I_U !N is the order of the matrix, IU controls largest eigenvalue returned if RANGE='I'. IU not ref'd if RANGE='V' or 'A'
            real(wp), intent(inout) :: A(n,n) ! A is the Hamiltonian matrix, on exit it is assigned to the eigenvecs
            real(wp), intent(out) :: W(n) ! Array of eigenvalues in ascending order
            integer(wp), intent(out)  :: M ! number of eigenvalues found.
            character*1 :: JOBV, UPLO ! JOBV controls whether eigenvecs and vals are calc'd (='V'), or just eigenvals (='N'). UPLO controls whether on exit from DSYEVR A stores upper or lower triangular matrix
            parameter  (JOBV='V', UPLO='U')

            integer(wp) :: LDA ! array leading dimension
            real(wp) :: VL = -3.0_wp ! lower and upper of interval to search for eigenvalues if RANGE='V'. Not accessed
            real(wp) :: VU = 20.0_wp
            integer(wp) :: IL = 1 ! lower bound of IL->IU if only selected eigenvals requested.
            real(wp) :: ABSTOL ! absolute error tolerance for eigenvals, 

            real(wp), external :: dlamchm ! utility func to determine machine parameters for minimum error tolerance without overflow
            real(wp), allocatable  :: Z(:,:) ! eigenvectors of A in columns
            integer(wp) :: LDZ ! leading dimension of z
            integer(wp) :: ISUPPZ(2*max(1,N))
            real(wp), allocatable :: WORK(:)
            real(wp) :: WORK_DIM(1)
            integer(wp) :: LWORK
            integer(wp), allocatable :: IWORK(:)
            integer(wp) :: IWORK_DIM(1)
            integer(wp) :: LIWORK
            integer(wp) :: INFO
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


end module
