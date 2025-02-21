! Contains global variables.
! Minimise the use of these objects as much as possible and only use for reading into the program.
! Where possible, use variables declared here as arguments for subroutines in other modules, rather than accessing directly.
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
    integer(wp) :: GEOMETRY_TYPE = 0 ! 0=straight wires along x, twist between consecutive chains controlled by twist angle, 1=twisted wire with dipoles pointing out of the wire,can be H aggregated as well, 2=helical coil
    integer(wp):: lattice_dimx = 1
    integer(wp):: lattice_dimy = 1
    integer(wp):: lattice_dimz = 1
    integer(wp):: configs = 1
    integer :: num_threads = 1
    logical :: bool_one_particle_states = .true.
    logical :: bool_two_particle_states = .true.
    logical :: bool_charge_transfer_states = .false.
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
    real(wp) :: theta = 0.0_wp ! the angle of precession of dipoles in helical twisted straight chains (type 1)
    real(wp) :: twist_angle = 0.0_wp ! angle between successive chromophores in successively twisted chains
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
    integer(wp) :: lattice_count = 0
    integer(wp) :: one_particle_counter = 0
    integer(wp) :: two_particle_counter = 0
    integer(wp) :: chargetransfer_counter = 0

    ! index arrays
    integer(wp), allocatable :: lattice_index_arr(:,:,:)
    integer(wp), allocatable :: one_particle_index_arr(:,:)
    integer(wp), allocatable :: two_particle_index_arr(:,:,:,:)
    integer(wp), allocatable :: chargetransfer_index_arr(:,:,:,:)

    ! franck-condon table
    real(wp), allocatable :: fc_ground_to_neutral(:,:)

    ! array of dipole moment vectors
    real(wp), allocatable :: mu_xyz(:,:)
    ! array of position vectors
    real(wp), allocatable :: r_xyz(:,:)

    integer(wp), allocatable :: bonding_matrix(:,:)
    ! hamiltonian, H
    ! The hamiltonian is a 2d array
    real(wp), allocatable :: H(:,:)
    real(wp), allocatable :: EVAL(:)
    integer(wp), parameter :: empty = -1

    integer :: IU = 1
    integer :: EVAL_COUNT


    ! File out names
    character*256 :: eval_out_f


    ! For property calculations
    complex(kind=wp), parameter :: complex_zero = ( 0_wp, 0_wp )
    complex(kind=wp), allocatable:: abs_osc_strengths_x(:)
    complex(kind=wp), allocatable:: abs_osc_strengths_y(:)
    complex(kind=wp), allocatable:: abs_osc_strengths_z(:)
    ! complex(kind=wp), allocatable:: abs_osc_strengths_x_configavg(:)
    ! complex(kind=wp), allocatable:: abs_osc_strengths_y_configavg(:)
    ! complex(kind=wp), allocatable:: abs_osc_strengths_z_configavg(:)
    real(wp), allocatable :: xpl_osc(:,:)
    real(wp), allocatable :: ypl_osc(:,:)
    real(wp), allocatable :: zpl_osc(:,:)
    ! real(wp), allocatable :: xpl_osc_configavg(:,:)
    ! real(wp), allocatable :: ypl_osc_configavg(:,:)  
    integer, parameter  :: spec_steps = 2600
    real(wp), allocatable :: pl_specx(:)
    real(wp), allocatable :: pl_specy(:)
    real(wp), allocatable :: pl_specz(:)
    real(wp), allocatable :: pl_specx_by_v(:,:)
    real(wp), allocatable :: pl_specy_by_v(:,:)
    real(wp), allocatable :: pl_specz_by_v(:,:)
    real(wp), allocatable :: pl_specx_by_v_perconfig(:,:)
    real(wp), allocatable :: pl_specy_by_v_perconfig(:,:)
    real(wp), allocatable :: pl_specz_by_v_perconfig(:,:)
    real(wp), allocatable :: abs_specx(:)
    real(wp), allocatable :: abs_specy(:)
    real(wp), allocatable :: abs_specz(:)
    real(wp), allocatable :: abs_specx_configavg(:)
    real(wp), allocatable :: abs_specy_configavg(:)
    real(wp), allocatable :: abs_specz_configavg(:)
    real(wp), allocatable :: cd_rot_strengths(:)
    real(wp), allocatable :: cd_spec(:)
    real(wp), allocatable :: cd_spec_avg(:)

    real(wp), dimension(2) :: rot_strengths
    real(wp), dimension(2) :: rot_strengths_tmp
    real(wp), allocatable :: cpl_spec(:)

end module
