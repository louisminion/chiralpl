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
    integer(wp):: lattice_dimx = 1
    integer(wp):: lattice_dimy = 1
    integer(wp):: lattice_dimz = 1
    integer(wp):: configs = 1
    logical :: bool_one_particle_states = .true.
    logical :: bool_two_particle_states = .true.
    logical :: H_out = .false.
    logical :: save_evals = .false.
    logical :: save_evecs = .false.
    logical :: manual_coupling = .true. ! if true, then doesn't calculate JCouplings individually, just uses inputs JCoulx,JCouly,JCoulz
    integer(wp) :: max_vibs
    integer(wp) :: n_nearest_neighbour = 1
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
    integer(wp) :: general_counter = 0
    integer(wp) :: lattice_count = 0
    integer(wp) :: one_particle_counter = 0
    integer(wp) :: two_particle_counter = 0


    ! index arrays
    integer(wp), allocatable :: lattice_index_arr(:,:,:)
    integer(wp), allocatable :: one_particle_index_arr(:,:)
    integer(wp), allocatable :: two_particle_index_arr(:,:,:,:)

    ! franck-condon table
    real(wp), allocatable :: fc_ground_to_neutral(:,:)

    ! array of dipole moment vectors
    real(wp), allocatable :: mu_xyz(:,:)

    ! hamiltonian, H
    ! The hamiltonian is a 2d array
    real(wp), allocatable :: H(:,:)
    real(wp), allocatable :: EVAL(:)
    integer(wp), parameter :: empty = -1

    integer(wp) :: IU = 1
    integer(wp) :: EVAL_COUNT


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
    !$OMP THREADPRIVATE(abs_osc_strengths_x,abs_osc_strengths_y,xpl_osc,ypl_osc,H,EVAL,diagonal_disorder_offsets)
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