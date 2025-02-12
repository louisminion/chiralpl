module hamiltonian
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: dipole_moment, coupling, build1particleHamiltonian, build2particleHamiltonian, build1particle2particleHamiltonian, Diagonalize

    contains
        subroutine dipole_moment
            ! Calculates the dipole moments for each chromophore in an aggregrate depending on the chosen geometry
            ! Can use global variables as not called in openmp loop
            use variables
            implicit none
            integer(wp) :: ix,iy,iz,ixyz, polymer_chain_index
            ! real(wp) :: counter
            allocate(mu_xyz(lattice_count,3))
            mu_xyz = 0.0_wp
            if (GEOMETRY_TYPE .eq. 0) then
                do ix=1,lattice_dimx
                    do iy=1,lattice_dimy
                        do iz=1,lattice_dimz
                            ixyz = lattice_index_arr(ix,iy,iz)
                            mu_xyz(ixyz,1) = mu_0*cos((iz*1.0_wp-1.0_wp)*phi)
                            mu_xyz(ixyz,2) = mu_0*sin((iz*1.0_wp-1.0_wp)*phi)
                        end do
                    end do
                end do

            else if (GEOMETRY_TYPE .eq. 1) then
                do iz=1,lattice_dimz
                    polymer_chain_index = 0
                    do iy=1,lattice_dimy
                        do ix=1,lattice_dimx
                            polymer_chain_index = polymer_chain_index + 1
                            ixyz = lattice_index_arr(ix,iy,iz)
                            mu_xyz(ixyz,1) = mu_0*(cos(twist_angle*polymer_chain_index)*cos(theta+(iz*1.0_wp-1.0_wp)*phi)+cos((iz*1.0_wp-1.0_wp)*phi)*cos(theta)-cos(twist_angle*polymer_chain_index)*cos((iz*1.0_wp-1.0_wp)*phi)*cos(theta))
                            mu_xyz(ixyz,2) = mu_0*(cos(twist_angle*polymer_chain_index)*sin(theta+(iz*1.0_wp-1.0_wp)*phi)+sin((iz*1.0_wp-1.0_wp)*phi)*cos(theta)-cos(twist_angle*polymer_chain_index)*sin((iz*1.0_wp-1.0_wp)*phi)*cos(theta))
                            mu_xyz(ixyz,3) = mu_0*(sin(twist_angle*polymer_chain_index)*sin(theta)*cos(2*((iz*1.0_wp-1.0_wp)*phi)))
                            ! mu_xyz(ixyz,1) = mu_0*(cos(theta)*cos((iz*1.0_wp-1.0_wp)*phi)-sin(theta)*sin((iz*1.0_wp-1.0_wp)*phi))
                            ! mu_xyz(ixyz,2) = mu_0*(sin(theta)*cos((iz*1.0_wp-1.0_wp)*phi)+cos(theta)*sin((iz*1.0_wp-1.0_wp)*phi))
                            ! mu_xyz(ixyz,3) = 0
                        end do
                    end do
                end do

            else if (GEOMETRY_TYPE .eq. 2) then ! helically coiled polymer, x dimension is the polymer length.
                ! counter = 0.0_wp
                do ix=1,lattice_dimx
                    ixyz = lattice_index_arr(ix,1,1)
                    ! sum_xpos = sum_xpos + x_spacing*cos(((ix*1.0_wp)-2.0_wp)*phi)
                    ! sum_ypos = sum_ypos + x_spacing*sin(((ix*1.0_wp)-2.0_wp)*phi)
                    mu_xyz(ixyz,1)= mu_0*(cos(((ix*1.0_wp)-2.0_wp)*phi))
                    ! write(*,*) 'MUX', cos(((ix*1.0_wp)-2.0_wp)*phi)
                    mu_xyz(ixyz,2)= mu_0*(sin(((ix*1.0_wp)-2.0_wp)*phi))
                    mu_xyz(ixyz,3) = mu_0*(cos(theta))
                end do

            end if 
        end subroutine

        real(wp) function coupling(x1,y1,z1,x2,y2,z2, mu_xyz,r_xyz,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon)
            integer(wp), intent(in) :: x1,y1,z1,x2,y2,z2
            integer(wp) :: d_x, d_y, d_z, i_xyz1, i_xyz2
            real(wp) :: dx,dy,dz
            real(wp), dimension(3) :: R, R_norm
            real(wp), intent(in) :: mu_xyz(:,:), r_xyz(:,:)
            integer(wp), intent(in) :: lattice_index_arr(:,:,:)
            real(wp), intent(in) :: x_spacing,y_spacing,z_spacing
            real(wp) :: R_mag
            real(wp) :: CT_couple_alongchain
            real(wp) :: dipole_dot_product
            real(wp) :: mu1_dotR
            real(wp) :: mu2_dotR
            real(wp) :: pi, epsilon
            logical, intent(in) :: manual_coupling
            coupling  = 0.0_wp
            i_xyz1 = lattice_index_arr(x1,y1,z1)
            i_xyz2 = lattice_index_arr(x2,y2,z2)


            d_x = abs(r_xyz(i_xyz2,1) - r_xyz(i_xyz1,1))
            d_y = abs(r_xyz(i_xyz2,2) - r_xyz(i_xyz1,2))
            d_z = abs(r_xyz(i_xyz2,3) - r_xyz(i_xyz1,3))
            ! if (manual_coupling .eq. .true.) then
            !     go to 78
            ! end if
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
            dipole_dot_product = mu_xyz(i_xyz1,1)*mu_xyz(i_xyz2,1) + mu_xyz(i_xyz1,2)*mu_xyz(i_xyz2,2) + mu_xyz(i_xyz1,3)*mu_xyz(i_xyz2,3)
            mu1_dotR = mu_xyz(i_xyz1,1)*R_norm(1) + mu_xyz(i_xyz1,2)*R_norm(2) + mu_xyz(i_xyz1,2)*R_norm(2)
            mu2_dotR = mu_xyz(i_xyz2,1)*R_norm(1) + mu_xyz(i_xyz2,2)*R_norm(2) + mu_xyz(i_xyz2,2)*R_norm(2)
            coupling = (1/(4*pi*epsilon*(R_mag**3)))*(dipole_dot_product-(3*mu1_dotR*mu2_dotR))
            return
            ! 78 if (d_y .eq. 0 .and. d_z .eq. 0) then
            !     if (d_x .eq. 1) then
            !         coupling = JCoulx !+ te_x + th_x
            !     end if
            !     return

            ! end if

            ! if (d_x .eq. 0 .and. d_z .eq. 0) then
            !     if (d_y .eq. 1) then
            !         coupling = JCouly !+ te_x + th_x
            !     end if
            !     return
            ! end if
            ! if (d_x .eq. 0 .and. d_y .eq. 0) then
            !     if (d_z .eq. 1) then
            !         coupling = JCoulz !+ te_x + th_x
            !     end if
            !     return

            ! end if
            return

        end function

        subroutine build1particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz,r_xyz,fc_ground_to_neutral,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,one_particle_index_arr,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
            implicit none
            real(wp), intent(inout) :: H(:,:)
            real(wp), intent(in) :: diagonal_disorder_offsets(:)
            integer(wp), intent(in) :: lattice_dimx,lattice_dimy,lattice_dimz,max_vibs
            real(wp), intent(in) :: x_spacing, y_spacing, z_spacing
            integer(wp), intent(in) :: one_particle_index_arr(:,0:)
            integer(wp), intent(in) :: lattice_index_arr(:,:,:)
            logical, intent(in) :: manual_coupling
            real(wp), intent(in) :: pi,epsilon,w00
            real(wp), intent(in) :: mu_xyz(:,:), r_xyz(:,:)
            real(wp), intent(in) :: fc_ground_to_neutral(0:,0:)
            integer(wp) :: empty = -1
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
                                            H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2, mu_xyz, r_xyz,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon)*fc_ground_to_neutral(0,vib_i1)*fc_ground_to_neutral(0,vib_i2)
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

        subroutine build2particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz,r_xyz,fc_ground_to_neutral,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,two_particle_index_arr,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
            implicit none
            real(wp), intent(inout) :: H(:,:)
            real(wp), intent(in) :: diagonal_disorder_offsets(:)
            integer(wp), intent(in) :: lattice_dimx,lattice_dimy,lattice_dimz,max_vibs
            real(wp), intent(in) :: x_spacing, y_spacing, z_spacing
            integer(wp), intent(in) :: two_particle_index_arr(:,0:,:,1:)
            integer(wp), intent(in) :: lattice_index_arr(:,:,:)
            logical, intent(in) :: manual_coupling
            real(wp), intent(in) :: pi,epsilon,w00
            real(wp), intent(in) :: mu_xyz(:,:), r_xyz(:,:)
            real(wp), intent(in) :: fc_ground_to_neutral(0:,0:)
            integer(wp) :: empty = -1
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
                                                                H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2, mu_xyz, r_xyz, manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon)* &
                                                                fc_ground_to_neutral(0, vib_i1)*fc_ground_to_neutral(0,vib_i2)
                                                                H(h_j,h_i) = H(h_i,h_j)
                                                            end if
                                                            if (i_xyz1v .eq. i_xyz2) then ! if vibrat exc in state 1 at same place as vibronic on state 2 then exchange type 
                                                                i_xyz2v = i_xyz1
                                                                do vib_i2v=1, max_vibs
                                                                    h_j = two_particle_index_arr(i_xyz2, vib_i2, i_xyz2v, vib_i2v)
                                                                    if ( h_j .eq. empty .or. h_j .eq. h_i ) cycle
                                                                    H(h_i, h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2, mu_xyz,r_xyz,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon)*fc_ground_to_neutral(vib_i2v, vib_i1)*fc_ground_to_neutral(vib_i1v,vib_i2)
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


        subroutine build1particle2particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz,r_xyz,fc_ground_to_neutral,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,one_particle_index_arr,two_particle_index_arr,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
            implicit none
            real(wp), intent(inout) :: H(:,:)
            real(wp), intent(in) :: diagonal_disorder_offsets(:)
            integer(wp), intent(in) :: lattice_dimx,lattice_dimy,lattice_dimz,max_vibs
            real(wp), intent(in) :: x_spacing, y_spacing, z_spacing
            integer(wp), intent(in) :: one_particle_index_arr(:,0:)
            integer(wp), intent(in) :: two_particle_index_arr(:,0:,:,1:)
            integer(wp), intent(in) :: lattice_index_arr(:,:,:)
            logical, intent(in) :: manual_coupling
            real(wp), intent(in) :: pi,epsilon,w00
            real(wp), intent(in) :: mu_xyz(:,:), r_xyz(:,:)
            real(wp), intent(in) :: fc_ground_to_neutral(0:,0:)
            integer(wp) :: empty = -1
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
                                                H(h_i,h_j) = coupling(i_x1,i_y1,i_z1,i_x2,i_y2,i_z2, mu_xyz, r_xyz,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon)*fc_ground_to_neutral(vib_i2v,vib_i1)*fc_ground_to_neutral(0,vib_i2)
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
            integer, intent(in) :: N !N is the order of the matrix 
            integer, intent(in) :: I_U !IU controls largest eigenvalue returned if RANGE='I'. IU not ref'd if RANGE='V' or 'A'
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
            ! write(*,*) SIZE(H), 'SIZE(H)'
            call dsyevr(JOBV,RANGE,UPLO,N,A,LDA,VL,VU,IL,I_U,ABSTOL,M,W,Z,LDZ,ISUPPZ,WORK_DIM,LWORK,IWORK_DIM,LIWORK,INFO)
            LWORK = WORK_DIM(1)
            LIWORK = IWORK_DIM(1)
            ! allocate WORK and IWORK arrays
            allocate(WORK(LWORK))
            ! write(*,*) LIWORK, 'LIWORK' 
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
