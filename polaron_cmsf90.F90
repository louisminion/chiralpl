!***************************************************************************************
                       !Module containing all the Variables for the program
!*************************************************************************************** 

module common_variables 
implicit none
 integer conf_max, conf_process
        real*8 dwidth, ycorrlen, xcorrlen, time_start, time_finish,    &
               delta_da,sigma
        real*8, allocatable :: offset(:,:)
          !lattice numbers
    integer :: lattice_x = 1
    integer :: lattice_y = 1
    integer :: lattice_z = 1
    integer :: n_proc = 1
    !mpi
    integer :: dim2 = 10
  
    !vibration
    real*8 :: hw = 1400.d0
    real*8 :: s  = 1.d0
    real*8 :: s_cat = 0.d0
    real*8 :: s_ani = 0.d0
    integer:: vibmax = 0

    
    !coupling
    real*8 :: jo_x = 0.d0
    real*8 :: jo_y = 0.d0
    real*8 :: jo_z = 0.d0
    real*8 :: ct_u = 0.d0
    real*8 :: ct_v = 0.d0
    real*8 :: te_x = 0.d0
    real*8 :: te_y = 0.d0
    real*8 :: te_z = 0.d0
    real*8 :: th_x = 0.d0
    real*8 :: th_y = 0.d0
    real*8 :: th_z = 0.d0

    !coupling
    logical :: extended_cpl = .false.
    real*8, allocatable :: jo_ex(:,:,:)
    character*256 extended_cpl_file

    !task title
    character*256 task_title

    !hamiltonian counters
    integer :: kount = 0
    integer :: icount = 1
    integer :: kount_1p = 0
    integer :: kount_2p = 0
    integer :: kount_3p = 0
    integer :: kount_ct = 0
    integer :: kount_ct2 = 0
    integer :: kount_lattice = 0
    integer :: nnz = 0
    integer :: sum1=0
     !doping
    real*8 :: dintra = 0.d0
    real*8 :: dinter = 0.d0
    real*8 :: r1 = 0.d0
    real*8 :: danion = 0.d0
    real*8 :: pos_anion = 0.d0
    real*8 :: binstart = 0.d0
    real*8 :: binend = 0.d0


    !indexes
    integer, allocatable :: nx_lattice(:,:,:)
    integer, allocatable :: nx_1p(:,:)
    integer, allocatable :: index_excited(:)
    integer, allocatable :: indx(:)
    integer, allocatable :: jndx(:)
    integer, allocatable :: nx_2p(:,:,:,:)
    real *8, allocatable :: disorder(:,:,:)  
    real, allocatable :: pdisorder(:,:,:,:,:)
    real *8, allocatable :: Force_anion(:,:,:) 
    real *8, allocatable :: Force_cation(:,:,:) 
    real *8, allocatable :: D_anion1(:,:,:)
    real *8, allocatable :: D_anion2(:,:,:)
    real *8, allocatable :: D_anion3(:,:,:)
    real *8, allocatable :: D_anion4(:,:,:)
    real *8, allocatable :: D_cation(:,:,:)    

    !parameterize the hamiltonian
    real*8, allocatable :: h(:,:)
    real*8, allocatable :: vl(:)
    real*8, allocatable :: eval(:)
    
    !franck condon factors
    real*8, allocatable :: fc_gf(:,:)
    real*8, allocatable :: fc_gc(:,:)
    real*8, allocatable :: fc_ga(:,:)
    real*8, allocatable :: fc_af(:,:)
    real*8, allocatable :: fc_cf(:,:)

    !constants
    real*8, parameter :: pi = 4.d0*datan(1.d0)
    real*8, parameter :: ev = 8065.d0
    real*8, parameter :: hc = 1.23984193d3 * ev !(plancks constant times the speed of light in nm*hw )
    real*8, parameter :: kb = 0.6956925d0        !units of cm-1 k
    real*8, parameter :: beta = 2.35
    !empty
    integer, parameter :: empty = -1

    !multiparticle states
    logical :: one_state    =.true.
    logical :: two_state    =.false.
    logical :: LL    =.true.
    logical :: SS    =.false.
    logical :: LS = .false.
    logical :: offdiagonal = .true.
    logical :: diagonal = .false.
    logical :: three_state  =.false.
    logical :: ct_state     =.false.
    logical :: ct2_state    =.false.
    logical :: dopant = .false.
    
    !periodic
    logical periodic
    
    !oscillator strength
    real*8, allocatable :: osc_x(:), osc_y(:), osc_z(:)
    real*8, allocatable :: plosc_x(:,:), plosc_y(:,:), plosc_z(:)
    real*8, allocatable :: ux(:), uy(:), uz(:)
    real*8 :: abs_lw = 100.d0
    real*8 :: mon_tr_e = 0.d0
    integer, parameter  :: spec_step = 2600        !for 10 cm-1 resolution

    !number of threads
    integer :: numthreads = 1
    
    !dielectric
    real*8 :: dielectric=1.d0
    
    !ct states trunctation
    logical :: ct_truncate = .false.
    logical :: cty_truncate = .false.
    integer :: ct_range = 1000
    integer :: cty_range = 1000
    integer :: pem_e = 10
    integer :: pem_g = 10
    logical :: two_truncate = .false.
    integer :: two_range = 1000
    
    logical :: lorentzian = .false.
    logical :: linewidthprogress = .false.
    real*8 :: delta_sigma = 0.d0
    real*8, allocatable :: vibsum_evec(:), vibsum_state(:)
    logical :: abs_freq_dep = .false.
    
    
    integer :: kount2
    integer :: iu=1
    character(1) :: rrange='A'

!***************************************************************************************
          !Petsc and Slepc variables and headers required for the program
!*************************************************************************************** 

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>        
#include <petsc/finclude/petscmat.h>
#include <slepc/finclude/slepcsys.h>
#include <slepc/finclude/slepceps.h>

!*************************************************************************************** 

#define xx_ag(ib)  xx_vg(xx_ig + (ib))
#define xx_ae(ib)  xx_ve(xx_ie + (ib))

 
       Mat            A                          ! Matrix A
       EPS            eps                        ! Eigen Problem Solver
       EPSType        tname                      ! Type of Eigen Problem Solver  
       PetscReal      error                       
       PetscRandom    rnd                        ! Petsc Random no. generator 
       PetscInt       i, Istart, Iend, j
       Vec            xr_g,xi_g
       Vec            xr_e,xi_e
       Vec            xr_g_0, xr_e_0
       PetscScalar    kr_g,ki_g,kr_gg,kr_ground,ki_ground
       PetscScalar    kr_excited,ki_excited
       PetscScalar    kr_e,ki_e,value,ki_gg
       PetscInt       nev, maxit, its, nconv,n,mpd
       PetscMPIInt    rank,color,row_rank,row_size,new_comm
       PetscErrorCode ierr
       PetscBool      flg,terse
       PetscOffset    xx_ig,xx_ie
       double precision xx_vg(1),xx_ve(1)

    end module

        
!***************************************************************************************
                        !main program starts 
!***************************************************************************************


program cms
use common_variables
implicit none
integer ista,ien,nproc
integer config,complete


    integer :: bin1 = 0
    integer :: mbin1 = 0
    integer :: bin2 = 0
    integer :: mbin2 = 0
    integer :: bin3 = 0
    integer :: mbin3 = 0
    integer :: bin4 = 0
    integer :: mbin4 = 0
    integer :: bin5 = 0
    integer :: mbin5 = 0
  
  
    real*8  start, finish
    real*8 :: ab_x(spec_step) = 0.D0
    real*8 :: mab_x(spec_step) = 0.D0
    real*8 :: ab_y(spec_step) = 0.D0
    real*8 :: mab_y(spec_step) = 0.D0
    real*8 :: pl_x(spec_step) = 0.D0
    real*8 :: mpl_x(spec_step) = 0.D0
    real*8 :: pl_y(spec_step) = 0.D0
    real*8 :: mpl_y(spec_step) = 0.D0
    real*8 :: sumoscx = 0.d0, sumoscy = 0.d0
    real*8 :: msumoscx = 0.d0, msumoscy = 0.d0
    real*8, allocatable :: cohfxn(:,:,:), cohtmp(:,:,:)
    real*8, allocatable :: mcohfxn(:,:,:)    
    character*128 sstring
    
  
!***************************************************************************************
                       !set the parameters (local)
!***************************************************************************************    
    call set_para()
    
   
!***************************************************************************************
                       !Indexing the lattice
!***************************************************************************************    
   

    call index_lattice()
    if ( one_state )   call index_1p()
    if ( two_state )   call index_2p()    
    
    !output the kount
    print*, '***************************'
    print*, 'kount   :', kount                      !total basis set size
    print*, 'kount 1p:', kount_1p                   ! Count the no. of one particle states
    print*, 'kount 2p:', kount_2p                   ! Count the no. of two particle states
    print*, '***************************'
    
    !Frank Condon tables
    call set_fctable()                         ! Sets the Frank Condon Tables
        
    
    !Output parameters
    call para_out_cms()                       !outputs the parameters
    
    ! Disorder values (diagonal and offdiagonal)
    
    ! if (diagonal) call disorder_elements()
   ! if (offdiagonal) call paracrystalline_SS()

    
!***************************************************************************************
                       !Programming with Petsc
!***************************************************************************************       

  call SlepcInitialize(PETSC_NULL_CHARACTER, ierr)      ! Initializes the slepc library, SlepcInitialize() calls the PetscInitialize(), should be called in the beginning of the program
  call MPI_COMM_SIZE(PETSC_COMM_WORLD, nproc, ierr)       ! calls the no. of processors 
  call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, ierr)      ! calls the rank of the processor
  color = rank/1
  call MPI_Comm_Split(PETSC_COMM_WORLD,color,rank,new_comm,ierr)
  call MPI_Comm_rank(new_comm,row_rank,ierr)               ! calls the rank of the new communicator
  call MPI_Comm_size(new_comm,row_size,ierr)               ! calls the size of each row of the new communicator
  call para_range(1,conf_max,n_proc,color,ista,ien)       ! This subroutine is called at the end of the program, I have explained everything there
  !n=16
  !call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-n', n, flg, ierr)      
      
! print*, 'ista', ista, 'ien', ien, 'conf_max', conf_max, 'color',color, &
 !        'rank',rank, 'n_proc', n_proc     
! create and set up the matrix
  
 
    if (diagonal) call disorder_elements()               ! Adds static diagonal disorder
    if (dopant) call doping_disorder()                    ! Adds the coulomb potential between hole and anion
    if (offdiagonal) call paracrystalline_SS()            ! Adds off-diagonal disorder
   
        do config = ista, ien       ! ista is the Start of iterations for rank iproc; ien is the end of iterations for rank iproc
       
        
      call MatCreate(new_comm,A,ierr)                                ! creates the matrix
      call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,kount,kount,ierr)         ! Sets the sizes
      call MatSetFromOptions(A,ierr)
      call MatSetUp(A,ierr)                             !Sets up the internal matrix data structures for the later use. 
     
           
!***************************************************************************************
                       !Build the hamiltonian
        
        if ( one_state ) call hamiltonian_1p_sparse(config)            ! creates the one particle sparse hamiltonian
        if ( two_state ) call hamiltonian_2p_sparse(config)             ! creates the two particle sparse hamiltonian
        if ( one_state .and. two_state ) call hamiltonian_1p2p_sparse(config) ! creates the one particle two particle sparse hamiltonian
        
        
!***************************************************************************************   
      
       
      
       if ( SS ) call short_short(config)
     
    ! Should be called after setting the values. Assembles all the matrix values.
  
      call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY,ierr)                    !  Begins assembling the matrix. This routine should be called after completing all calls to MatSetValues(). 
      call MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY,ierr)                       !  Completes assembling the matrix. This routine should be called after MatAssemblyBegin(). 
       
     ! call MatView(A, PETSC_VIEWER_STDOUT_WORLD, ierr)        !optional
      

!***************************************************************************************
                       !Programming with SLEPC
!***************************************************************************************  

!   ** Create eigensolver context
      call EPSCreate(new_comm,eps,ierr)                  ! create the eigen problem solver

 !    ** Set operators. In this case, it is a standard eigenvalue problem 
      call EPSSetOperators(eps,A,PETSC_NULL_OBJECT,ierr)          
      call EPSSetProblemType(eps,EPS_HEP,ierr)                   ! My matrix is hermitian hence EPS_HEP
      nev = 1000                                               ! it is the total no. of eigen values   
      nconv = 2*nev                                              ! No. of converged eigen values
      mpd = nconv                                                
      call EPSSetDimensions(eps,nev,nconv,mpd,ierr)              ! We can set how many eigen values we need, if this command is not there, the program by default finds only one eigen pair           

 !     ** Set solver parameters at runtime
      call EPSSetFromOptions(eps,ierr)

 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 !     Solve the eigensystem
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       call EPSSetWhichEigenpairs(eps,eps_smallest_real,ierr)     ! Sort the eigen values according to priority. If not done it store the eigen values randomly, I need from smallest to highest
     !  call EPSSetType(eps,EPSLAPACK,ierr)
       call cpu_time(start)
       call EPSSolve(eps,ierr)                                    ! Solve
       call cpu_time(finish)

       print*, 'Diagonalization time:', finish-start 
       call EPSGetConverged (eps,nconv,ierr)                      ! Find the no. of converged eigenpairs
      if (rank .eq. 0) then
       write(*,150) nconv
     endif 
      150  format (' Number of converged eigenpairs:',I2) 
         
       

!    ** Optional: Get some information from the solver and display it
       call EPSGetType(eps,tname,ierr)
      if (rank .eq. 0) then
       write(*,120) tname
     endif
  120  format (' Solution method: ',A)
       call EPSGetDimensions(eps,nev,PETSC_NULL_INTEGER,                 &
     &                      PETSC_NULL_INTEGER,ierr)
      if (rank .eq. 0) then
      write(*,130) nev
       endif
  130  format (' Number of requested eigenvalues:',I4)
    
    
     !call VecGetArray(xr,xx_v,xx_i,ierr)
     
     
     !do j=0,nconv-1
     call EPSGetEigenpair(eps,0,kr_g,ki_g,xr_g,xi_g,ierr)         ! Get only the ground state eigen pair
     call density_of_groundstates(bin1,bin2,bin3,bin4,bin5)
     
     call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xr_g,ierr)   ! creates parallel vectors
     call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xi_g,ierr)
     call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xr_e,ierr)
     call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xi_e,ierr)

    if (row_rank == 0) then  
    call VecGetArray( xr_g,xx_vg,xx_ig,ierr) !real
   end if
 !***************************************************************************************
                       !Calculating the absorption spectra
                       
       call cms_osc()
       call cms_osc_sum( sumoscx, sumoscy, config)
       call cms_spec( ab_x, ab_y )
       
 !***************************************************************************************  
 
 !***************************************************************************************
                       !Calculating the photoluminiscence spectra
                       
      ! call pl_osc()
       call pl_osc1(config)
       call pl_spec( pl_x, pl_y )
       
 !***************************************************************************************  

 !***************************************************************************************
                       !Calculating coherence function
      
         if (.not. allocated( cohfxn ) ) then 
            allocate( cohfxn( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
            cohfxn = 0.d0
        end if 

        if (.not. allocated( mcohfxn ) ) then 
            allocate( mcohfxn( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
            mcohfxn = 0.d0
        end if       
        
       !allocate the temp function and initialize to zero every time!
        if (.not. allocated( cohtmp ) )                 &
            allocate( cohtmp( -lattice_x+1:lattice_x-1, &
                              -lattice_y+1:lattice_y-1, &
                              -lattice_z+1:lattice_z-1 ) )
        cohtmp = 0.d0
       
        call cms_cohfxn( cohtmp )
        cohfxn = cohfxn + cohtmp
!***************************************************************************************   


call EPSDestroy(eps,ierr)
call VecDestroy(xr_g,ierr)
call VecDestroy(xi_g,ierr)
call VecDestroy(xr_e,ierr)
call VecDestroy(xi_e,ierr)
call MatDestroy(A,ierr)
end do
  
     if (rank == 0) then
     print*, '***************************'
     print*, ' Finishing Up...'
     print*, '***************************'
     end if
    
  ! Please look into how mpi reduce works
  
     Call MPI_REDUCE(bin1, mbin1, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     bin1 = mbin1
     
     Call MPI_REDUCE(bin2, mbin2, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     bin2 = mbin2
     
      Call MPI_REDUCE(bin3, mbin3, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD,ierr)
     bin3 = mbin3
     
     Call MPI_REDUCE(bin4, mbin4, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     bin4 = mbin4
  
     Call MPI_REDUCE(bin5, mbin5, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     bin5 = mbin5
      
     Call MPI_REDUCE(ab_x, mab_x, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     ab_x = mab_x/( conf_max * 1.d0)


     Call MPI_REDUCE(ab_y, mab_y, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     ab_y = mab_y/( conf_max * 1.d0)
     
     Call MPI_REDUCE(pl_x, mpl_x, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     pl_x = mpl_x/( conf_max * 1.d0)


     Call MPI_REDUCE(pl_y, mpl_y, spec_step, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     pl_y = mpl_y/( conf_max * 1.d0)

     Call MPI_REDUCE(sumoscx, msumoscx, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     sumoscx = msumoscx/(conf_max * 1.d0 )*hw
          
     Call MPI_REDUCE(sumoscy, msumoscy, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     sumoscy = msumoscy/( conf_max * 1.d0 )*hw

     Call MPI_REDUCE(cohfxn, mcohfxn, ((lattice_x*2-1)*(lattice_y*2-1)), MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
     cohfxn = mcohfxn/( conf_max * 1.d0 )
   !Normalize
 
    if (rank == 0) then 
   ! ab_x = ab_x/( conf_max * 1.d0)
   ! ab_y = ab_y/( conf_max * 1.d0)
   ! cohfxn = cohfxn/( conf_max * 1.d0 )
   ! sumoscx = sumoscx/( conf_max * 1.d0 )*hw
   ! sumoscy = sumoscy/( conf_max * 1.d0 )*hw    !multiply by hw to keep in units of (cm-1)
    print*,'             row_rank  : ', row_rank
    print*,'             color  : ', color
    print*,            'Absorption Sums:'
    print*,'             x  : ', sumoscx
    print*,'             y  : ', sumoscy
    print*,'             sum: ', sumoscx + sumoscy
    write( sstring, * ) 'Absorption Sum (xysum),', sumoscx, ',', sumoscy, ',', sumoscx + sumoscy
    print*, 'bin1:', bin1
    print*, 'bin2:', bin2
    print*, 'bin3:', bin3
    print*, 'bin4:', bin4
    print*, 'bin5:', bin5
    !call para_out_append( sstring )

   ! !write the spectra (local)  
    call cms_out( ab_x, ab_y )
    call pl_out( pl_x, pl_y )
    call cohfxn_out( cohfxn )
    
    print*, ' Done'
    
 end if


      
       
       
     
     call MPI_Comm_free(new_comm,ierr)
     call SlepcFinalize(ierr)
   
    
end program

!***************************************************************************************
                        !main program ends
!***************************************************************************************



!***************************************************************************************
                        !Subroutines for the program
!***************************************************************************************

!***************************************************************************************!
!    Set the parameters
!***************************************************************************************!

subroutine set_para()
     use common_variables
    implicit none

    logical         exists
    character*100   buffer, label, fname
    integer         fno, ios, line, pos, errstat
    parameter       ( fno = 90 )

    !Set the default parameters
    print*, 'Setting the default parameters.'
    lattice_x       = 9                     
    lattice_y       = 1 
    vibmax          = 4
    hw              = 1400.d0
    s               = 1.d0
    jo_x            = -0.40d0     !These are really te or th, but using exciton subroutines
    jo_y            = -0.15d0     !I have to use the jo variables
    delta_da        = 0.d0        !donor - acceptor energy difference. Set to 0 for homopolymer
    task_title      = 'test'
    abs_lw          = 350.d0
    lorentzian      = .false.
    rrange          = 'A'         !By default, find all of the eigenstates
    iu              = 1           
   ! rrange          = 'I'
    !iu              = 10
    conf_max        = 1 !disorder conf parameters
    dwidth          = 0.d0
    sigma           = 0.20d0
    xcorrlen        = 0.d0
    ycorrlen        = 0.d0
    one_state       = .true.
    two_state       = .true.
    LL = .true.
    LS = .false.
    SS = .true.
    dintra = 0.4
    dinter = 0.4
    danion = 0.40
    pos_anion = 1.6
    binstart = 0
    binend = 2
    r1 = 0.20
    pem_e = 405
    pem_g = 404
    n_proc = 20
    dopant = .true.
    diagonal = .false.
    offdiagonal = .true.
    two_truncate    = .false.
    two_range       = .false.
    goto 1010
    !Read the name of the input file
    call get_command_argument(1, fname, status=errstat )
    if ( errstat .ne. 0 ) goto 1010
    inquire( file = trim(fname), exist = exists)
    if ( .not. exists ) then
        print*, 'Input file not found...aborting'
        stop
    end if
    
    !Open the file to read in
    open( unit = fno, file = fname, status = 'old', action = 'read' )
    !-----------------------------------!
    !   This reads the control file
    !   The ios changes if end of record
    !   or end of file is reached
    !-----------------------------------!
    ios = 0
    line = 0
    print*, 'Reading the input file...'
    do while ( ios == 0 )
        read( fno, '(a)', iostat=ios ) buffer
        if( ios == 0 ) then
            line = line + 1

            !Find first instance of whitespace
            pos = scan( buffer, ' ' )
            label = buffer( 1:pos )
            buffer = buffer( pos + 1:)

            select case ( label )
            case('lattice_x')
                read( buffer, *, iostat=ios ) lattice_x
                print*, '    Setting lattice_x to: ', lattice_x
            case('lattice_y')
                read( buffer, *, iostat=ios ) lattice_y
                print*, '    Setting lattice_y to: ', lattice_y
            case('vibmax')
                read( buffer, *, iostat=ios ) vibmax
                print*, '    Setting vibmax to: ', vibmax
            case('hw')
                read( buffer, *, iostat=ios ) hw
                print*, '    Setting vibrational energy to (cm-1): ', hw
            case('s')
                read( buffer, *, iostat=ios ) s
                print*, '    Setting s to : ', s
            case('te_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting te_x to (eV): ', jo_x
            case('th_x')
                read( buffer, *, iostat=ios) jo_x
                print*, '    Setting th_x to (eV): ', jo_x    
            case('delta_da')
                read( buffer, *, iostat=ios) delta_da
                print*, '    Setting delta_da to (cm-1): ', delta_da           
            case('te_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting te_y to (eV): ', jo_y
            case('th_y')
                read( buffer, *, iostat=ios) jo_y
                print*, '    Setting th_y to (eV): ', jo_y
            case('task_title')
                read( buffer, *, iostat=ios) task_title
                print*, '    Setting task_title to: ', trim(task_title)
            case('abs_lw')
                read( buffer, *, iostat=ios) abs_lw
                print*, '    Setting the linewidth to (cm-1): ', abs_lw
            case('two_range')
                read( buffer, *, iostat=ios) two_range
                print*, '    Setting the vibrational radius to: ', two_range
            case('two_truncate')
                read( buffer, *, iostat=ios) two_truncate
                if ( two_truncate ) print*, '    Vibrational radius is on'
                if ( .not. two_truncate ) print*, '    Vibrational radius is off'
            case('lorentzian')
                read( buffer, *, iostat=ios) lorentzian
                if ( lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' lorentzian'
                if ( .not.lorentzian ) print*, '    Setting the lineshape to:', &
                                          ' gaussian'
            case('one_state')
                read( buffer, *, iostat=ios) one_state
                if ( one_state ) print*, '    One particle states are turned on.'
                if ( .not.one_state ) print*, '    One particle states are turned off.'
            case('two_state')
                read( buffer, *, iostat=ios) two_state
                if ( two_state ) print*, '    Two particle states are turned on.'
                if ( .not.two_state ) print*, '    Two particle states are turned off.'
            case('rrange')
                read( buffer, *, iostat=ios) rrange
                if ( rrange .eq. 'A' ) print*, '    dsyevr will find all  eigenvectors'
                if ( rrange .eq. 'I' ) print*, '    dsyevr will find a range of eigenvectors'
            case('iu')
                read( buffer, *, iostat=ios) iu
                print*, '    desyevr will look for the lowest', iu, ' eigenvectors'           
            case('conf_max')
                read( buffer, *, iostat=ios) conf_max
                print*, '    Setting conf_max to: ', conf_max
            case('dwidth')
                read( buffer, *, iostat=ios) dwidth
                print*, '    Setting disorder width to (cm-1): ', dwidth
            case('xcorrlen')
                read( buffer, *, iostat=ios) xcorrlen
                print*, '    Setting xcorrlen length to: ', xcorrlen            
            case('ycorrlen')
                read( buffer, *, iostat=ios) ycorrlen
                print*, '    Setting ycorrlen length to: ', ycorrlen
            case default
                print*, '    invalid label at line, ', line
            end select
        end if

    end do

    close(fno)

1010 continue
    print*, 'Calculating derived parameters.'
    !Normalize parameters to hamiltonian units of vib quanta
    jo_x = jo_x * eV / hw
    jo_y = jo_y * eV / hw
    abs_lw = abs_lw / hw
    dwidth = dwidth / hw
    delta_da = delta_da / hw
    if ( vibmax == 0 ) s = 0.d0
    
!    call omp_set_num_threads(1)

end subroutine
!***************************************************************************************!
!    Write the parameters to a file
!***************************************************************************************!
subroutine para_out_cms()
    use common_variables
        implicit none
    
    character*64 fname
    integer fno
    parameter( fno = 16 )
    character*8 ddd
    character*10 ttt
    
    !create the file
    fname = trim(task_title)//'_para.csv'
    open( unit = fno, file = fname, action='write' )
    call date_and_time( ddd, ttt )
    
    write( fno, * ) 'parameter file for, cms_local_disorder.f95'
    write( fno, * ) 'task title, ', trim(task_title)
    write( fno, * ) 'date and time, ', ddd, ' ', ttt
    write( fno, * ) 'x dimension, ', lattice_x
    write( fno, * ) 'y dimension, ', lattice_y
    write( fno, * ) 'vibmax, ', vibmax
    write( fno, * ) 'vib energy (cm-1), ', hw 
    write( fno, * ) 'huang-rhys, ', s
    write( fno, * ) 'th x (eV), ', jo_x*hw/eV
    write( fno, * ) 'th y (eV), ', jo_y*hw/eV    
    write( fno, * ) 'abs lw (cm-1), ', abs_lw * hw
    if ( lorentzian ) write( fno, * ) 'lineshape, lorentzian'
    if ( .not. lorentzian ) write( fno, * ) 'lineshape, gaussian'    
    write( fno, * ) 'conformations, ', conf_max
    write( fno, * ) 'disorder width (cm-1), ', dwidth*hw
    write( fno, * ) 'x correlation length, ', xcorrlen
    write( fno, * ) 'y correlation length, ', ycorrlen
    write( fno, * ) 'one state, ', one_state
    write( fno, * ) 'two state, ', two_state
    write( fno, * ) 'two truncate, ', two_truncate
    write( fno, * ) 'two range, ', two_range   
    
    close( fno )
    
end subroutine

!***************************************************************************************!
!    Calculate the CMS oscillator strength
!***************************************************************************************!
subroutine cms_osc()
    use common_variables
    implicit none
    
    
     
    integer ground
    parameter ( ground = 0 )
    integer lx, ly, lz, lxyz, vib, hx,  &
            lxv, lyv, lzv, lxyzv, vibv, &
            excited
    real*8  tmp_x, tmp_y,mtmp_x, mtmp_y
    VecScatter ctxr_e, ctxr_g
    
    ! call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xr_g,ierr)   ! creates parallel vectors
    ! call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xi_g,ierr)
    ! call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xr_e,ierr)
    ! call VecCreateMPI(new_comm,PETSC_DECIDE,KOUNT,xi_e,ierr)
 if (row_rank == 0) then 
     call EPSGetEigenpair(eps,0,kr_g,ki_g,xr_g,xi_g,ierr)         ! Get only the ground state eigen pair
end if

    If ( .not. allocated( osc_x ) ) allocate( osc_x(1:nev), osc_y(1:nev) )
    
   
   !  if (row_rank == 0) then  
   ! call VecGetArray( xr_g,xx_vg,xx_ig,ierr) !real
  ! end if
 !    call EPSGetEigenvalue(eps, ground, kr_g, ki_g, ierr)
    do excited = 1, nev-1
        call EPSGetEigenpair(eps,excited,kr_e,ki_e,xr_e,xi_e,ierr)
                  
     if (row_rank == 0) then 
   !  if (kr_g > -5 .and. kr_g < -4) then 
    call VecGetArray(xr_e,xx_ve,xx_ie,ierr)

        !initialize the tmp variables
        tmp_x = 0.d0
        tmp_y = 0.d0
    
        !go over all 1p basis states
        if ( one_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            hx = nx_1p( lxyz, vib )
            if ( hx == empty ) cycle
           
            !set the transition dipole moments
            tmp_x = tmp_x + lx * xx_ag(hx) * xx_ae(hx )
            tmp_y = tmp_y + ly * xx_ag(hx) * xx_ae(hx )
    
        end do
        end do
        end do
        end do
        end if
        
        !go over all 2p basis states
        if ( two_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
        do vibv = 1, vibmax
            lxyzv = nx_lattice( lxv, lyv, lzv )
            hx = nx_2p( lxyz, vib, lxyzv, vibv )
            if ( hx == empty ) cycle
           
            !set the transition dipole moments
            tmp_x = tmp_x + lx * xx_ag(hx) * xx_ae(hx )
            tmp_y = tmp_y + ly * xx_ag(hx) * xx_ae(hx )
     
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end if

       
        !square to get the oscillator strength
        osc_x( excited ) = tmp_x * tmp_x
        osc_y( excited ) = tmp_y * tmp_y
!end if
    call VecRestoreArray( xr_e,xx_ve,xx_ie,ierr)
end if

    end do
  !  print*, 'osc_x', osc_x(excited)
 !   print*, 'tmp_x_osc', tmp_x

  if (row_rank == 0) then 
   call VecRestoreArray( xr_g,xx_vg,xx_ig,ierr )
  end if
    
end subroutine
!***************************************************************************************!
!    Calculate the CMS oscillator strength sums
!***************************************************************************************!
subroutine cms_osc_sum( sumx, sumy,config )
    use common_variables
    implicit none
    
    real*8, intent(inout) :: sumx, sumy
    integer excited, ground,config
    parameter ( ground = 0 )
    real*8 trans_e
    
   
! if (rank == 0) then   
  call EPSGetEigenvalue(eps, 0, kr_g, ki_g, ierr) 
     
      do excited =1,nev-1

        call EPSGetEigenvalue(eps, excited, kr_e, ki_e, ierr)
   

        

        trans_e = kr_e - kr_g
     !   if (trans_e < 2) then
     !   sum1 = sum1 + 1
        write(15,*) 'trans_e', trans_e, 'excited', excited
      !  end if
       !  if (trans_e > 2) then
      !  write(17,*) 'trans_e', trans_e, 'config', config
      !  end if

       !      write (15,*) 'trans_e', trans_e, 'excited', excited
 
     end do
        
    

!end if




!print*, sum1
end subroutine
!***************************************************************************************!
!    Calculate the CMS spectrum
!***************************************************************************************!
subroutine cms_spec( ab_x, ab_y )
    use common_variables
    implicit none
    
    real*8, intent(out) :: ab_x( spec_step ), ab_y( spec_step )
    integer ground
    parameter ( ground = 0 )
    integer spec_point, excited
    real*8 spec_start, spec_end, energy, trans_e, &
           lineshape, gamma
               
    spec_start = 0.d0                              !start
    spec_end = spec_start + 20000.D0/hw            !end
      
    do spec_point = 1, spec_step 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
      
        call EPSGetEigenvalue(eps, ground, kr_g, ki_g, ierr)            
          
        do excited = 1, nev-1
        call EPSGetEigenpair(eps,excited,kr_e,ki_e,xr_e,xi_e,ierr)
 if (row_rank == 0) then 
      
          trans_e = kr_e - kr_g
!print*, 'trans_e', trans_e
                                            
            gamma = abs_lw
            
            if ( lorentzian ) then
                lineshape = gamma/( (energy - trans_E)**2 + gamma**2 )/pi
!                lineshape = lineshape/(1.D0*kount_lattice)
            else                !Gaussian
                lineshape = dexp( - ( ( energy - trans_E ) / gamma ) ** 2 )
                lineshape = lineshape / ( dsqrt(pi) * gamma )            !Keep gaussian normalized
!                lineshape = lineshape / ( dsqrt(pi) * gamma * kount_lattice)
            end if
 
            !cms is frequency dependent, so must multiply by the transition energy
            lineshape = lineshape*trans_e
            
 
            ab_x( spec_point ) = ab_x( spec_point ) + lineshape * osc_x( excited )
            ab_y( spec_point ) = ab_y( spec_point ) + lineshape * osc_y( excited )
  end if

        end do        
    end do
  
end subroutine
!***************************************************************************************!
!    Write the CMS spectrum
!***************************************************************************************!
subroutine cms_out( ab_x, ab_y )
    use common_variables
    implicit none
    
    character(256) f_name
    integer f_no, spec_point
    real*8 spec_start, spec_end, energy,                                 &
           ab_x( spec_step ), ab_y( spec_step )
        
    !File info
    f_name = trim(task_title)//'_cms1.dat'
    f_no = 2
        
    !Open the file and put headers
    open( unit = f_no, file = f_name )
    write(f_no, *) 'Task Title: ', trim(task_title)
    write(f_no, *) 'CMS Data'
    write(f_no, *) 'Energy CMS CMS '
    write(f_no, *) ' cm\+(-1) a.u. a.u. a.u. '
    write(f_no, *) ' n.a. sum x y'
        
    spec_start = 0.d0                              !start
    spec_end = spec_start + 20000.D0/hw            !end
    
    do spec_point = 1, spec_step 
 if (row_rank == 0) then 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
        !multiply all by hw to get in units of cm-1
        write( f_no, '(5f14.7)' ) energy,                                           &
                                  (ab_x( spec_point ) + ab_y( spec_point )), &
                                  ab_x( spec_point ) , ab_y( spec_point )

        ! write( f_no, '(5f14.7)' ) energy,
                                  

end if   
 end do
    close( f_no )

 End Subroutine


!***************************************************************************************!
!    Calculate the photoluminiscence oscillator strength
!***************************************************************************************!
subroutine pl_osc1(config)
    use common_variables
    implicit none
    
    
     
    integer ground,config
    parameter ( ground = 0 )
    integer lx, ly, lz, lxyz, vib,vib1, hx,  &
            lxv, lyv, lzv, lxyzv, vibv, &
            excited,emitting,state
    real*8  tmp_x, tmp_y,mtmp_x, mtmp_y, energy, energy2
    VecScatter ctxr_e, ctxr_g
  

If ( .not. allocated( plosc_x ) ) allocate( plosc_x(0:nev,0:nev), plosc_y(0:nev,0:nev) )

     call EPSGetEigenvalue(eps, 0, kr_ground, ki_ground, ierr)     

do excited = 85,85
     call EPSGetEigenvalue(eps, excited, kr_excited, ki_excited, ierr)
     energy = (kr_excited - kr_ground)/5.76

!print*, 'energy', energy, 'excited', excited, 'config', config
if (energy > binstart .and. energy < binend) then
!print*, 'energy', energy, 'excited', excited, 'config', config


 call EPSGetEigenpair(eps,excited,kr_e,ki_e,xr_e,xi_e,ierr)
     do state=0,excited-1
 call EPSGetEigenpair(eps,state,kr_g,ki_g,xr_g,xi_g,ierr)
     


    call VecGetArray(xr_e,xx_ve,xx_ie,ierr)
    call VecGetArray( xr_g,xx_vg,xx_ig,ierr)

        !initialize the tmp variables
        tmp_x = 0.d0
        tmp_y = 0.d0
    
        !go over all 1p basis states
        if ( one_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            hx = nx_1p( lxyz, vib )
            if ( hx == empty ) cycle
           

           !set the transition dipole moments
            tmp_x = tmp_x + lx * xx_ag(hx) * xx_ae(hx )
            tmp_y = tmp_y + ly * xx_ag(hx) * xx_ae(hx )
        
        end do
        end do
        end do
        end do
        end if                            !if statement for 1p states end
        
        !go over all 2p basis states
        if ( two_state ) then
        do lx = 1, lattice_x
        do ly = 1, lattice_y
        do lz = 1, lattice_z
        do vib = 0, vibmax
            lxyz = nx_lattice( lx, ly, lz )
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
        do vibv = 1, vibmax
            lxyzv = nx_lattice( lxv, lyv, lzv )
            hx = nx_2p( lxyz, vib, lxyzv, vibv )
            if ( hx == empty ) cycle
           
            !set the transition dipole moments
            tmp_x = tmp_x + lx * xx_ag(hx) * xx_ae(hx )
            tmp_y = tmp_y + ly * xx_ag(hx) * xx_ae(hx )
     
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end if                       !if statement for 2p states end

plosc_x( excited, state ) = tmp_x * tmp_x
plosc_y( excited, state ) = tmp_y * tmp_y

!print *, 'plosc_x(',excited, state, ')', plosc_x( excited, state)


    end do
end if
end do


  if (row_rank == 0) then 
   call VecRestoreArray( xr_e,xx_ve,xx_ie,ierr)
   call VecRestoreArray( xr_g,xx_vg,xx_ig,ierr )
  end if
    
end subroutine

!***************************************************************************************!
!    Calculate the photoluminiscence spectrum
!***************************************************************************************!
subroutine pl_spec( pl_x, pl_y )
    use common_variables
    implicit none
    
    real*8, intent(out) :: pl_x( spec_step ), pl_y( spec_step )
    integer ground
    parameter ( ground = 0 )
    integer spec_point, excited,vib,state
    real*8 spec_start, spec_end, energy, trans_e, &
           lineshape, gamma, energy1
     
    spec_start = 0.d0                              !start
    spec_end = spec_start + 20000.D0/hw            !end
      
    do spec_point = 1, spec_step 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
 
      
do excited = 85, 85
 call EPSGetEigenvalue(eps, excited, kr_ground, ki_ground, ierr)   

do state = 0, excited-1
     call EPSGetEigenvalue(eps, state, kr_excited, ki_excited, ierr)
     energy1 = (kr_ground - kr_excited)/5.76
    ! trans_e = kr_excited - kr_ground
    trans_e = kr_ground - kr_excited

if (energy1 > binstart .and. energy1 < binend) then

     
  
               gamma = abs_lw
            
            if ( lorentzian ) then
                lineshape = gamma/( (energy - trans_e)**2 + gamma**2 )/pi
!                lineshape = lineshape/(1.D0*kount_lattice)
            else                !Gaussian
                lineshape = dexp( - ( ( energy - trans_e ) / gamma ) ** 2 )
                lineshape = lineshape / ( dsqrt(pi) * gamma )            !Keep gaussian normalized
!                lineshape = lineshape / ( dsqrt(pi) * gamma * kount_lattice)
            end if
 
            !cms is frequency dependent, so must multiply by the transition energy
            lineshape = lineshape*trans_e*trans_e*trans_e
            pl_x( spec_point ) = pl_x( spec_point ) + lineshape * plosc_x( excited,state)
            pl_y( spec_point ) = pl_y( spec_point ) + lineshape * plosc_y( excited,state)
end if


 end do
     
    end do
end do
   
end subroutine

!***************************************************************************************!
!    Write the photoluminiscence spectrum
!***************************************************************************************!
subroutine pl_out( pl_x, pl_y )
    use common_variables
    implicit none
    
    character(256) f_name
    integer f_no, spec_point
    real*8 spec_start, spec_end, energy,                                 &
           pl_x( spec_step ), pl_y( spec_step )
        
    !File info
    f_name = trim(task_title)//'_pl.dat'
    f_no = 3
      
    !Open the file and put headers
    open( unit = f_no, file = f_name )
    write(f_no, *) 'Task Title: ', trim(task_title)
    write(f_no, *) 'CMS Data'
    write(f_no, *) 'Energy CMS CMS '
    write(f_no, *) ' cm\+(-1) a.u. a.u. a.u. '
    write(f_no, *) ' n.a. sum x y'
        
    spec_start = 0.d0                              !start
    spec_end = spec_start + 20000.D0/hw            !end
    
    do spec_point = 1, spec_step 
 if (row_rank == 0) then 
        energy = spec_start + spec_point/( spec_step * 1.D0 )*( spec_end - spec_start )
        !multiply all by hw to get in units of cm-1
        write( f_no, '(5f14.7)' ) energy,                                           &
                                  (pl_x( spec_point ) + pl_y( spec_point )), &
                                  pl_x( spec_point ) , pl_y( spec_point )
end if   
 end do
    close( f_no )

 End Subroutine

!***************************************************************************************!
!Coherence function
!***************************************************************************************!
subroutine cms_cohfxn( cohfxn )
    use common_variables
    implicit none
    
    real*8, intent(inout) :: cohfxn( -lattice_x+1:lattice_x-1, &
                                     -lattice_y+1:lattice_y-1, &
                                     -lattice_z+1:lattice_z-1 )
    integer lx, ly, lz, lxyz,     &
            rx, ry, rz,           &
            lx2, ly2, lz2, lxyz2, &
            vib, vib2, ground,    &
            hx, hx2,              &
            lxv, lyv, lzv, lxyzv, &
            vibv, vibg,           &
            lxyz2v, vib2v
    parameter ( ground = 1, vibg = 0 )
  
    !choose r
    do rx = -lattice_x + 1, lattice_x - 1
    do ry = -lattice_y + 1, lattice_y - 1
    do rz = -lattice_z + 1, lattice_z - 1
    
        !one-particle one-particle
        if ( one_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
         !   write(14,*)'rx:',rx,'ry:',ry
         !   write(14,*) 'lx:', lx, 'ly:', ly, 'lx2:',lx2, 'ly2:', ly2
         !   write(14,*) 'lxyz:', lxyz, 'lxyz2:', lxyz2
           
            !if on the same molecule, vib must be the same on both
            if ( lxyz == lxyz2 ) then
                do vib = 0, vibmax
                    hx = nx_1p( lxyz, vib )
                 !   write(14,*) 'hx:', hx
                    if ( hx == empty ) cycle
                    cohfxn( rx, ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                           xx_ag(hx) * xx_ag(hx )
             !   write(14,*) 'cohfxn(', rx, ry, rz ,')', cohfxn(rx,ry,rz)
                 end do
            !if not multiply by the FC factors
            else
                do vib = 0, vibmax
                do vib2 = 0, vibmax
                    hx = nx_1p( lxyz, vib )
                    if ( hx == empty ) cycle
                    hx2 = nx_1p( lxyz2, vib2 )
                    if ( hx2 == empty ) cycle
                    cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                          xx_ag(hx) * xx_ag(hx2) *   &
                                          fc_gf( 0, vib ) * fc_gf( 0, vib2 )
                end do
                end do
            end if
        end do
        end do
        end do
        end if
        
        !one-particle two-particle
        if ( one_state .and. two_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            if ( lxyz == lxyz2 ) cycle            !if on the same molecule, this is zero
            lxyz2v = lxyz                         !free vibration must be on hole site in 1p state
            do vib = 0, vibmax
            do vib2 = 0, vibmax
            do vib2v = 1, vibmax
                hx = nx_1p( lxyz, vib )
                if ( hx == empty ) cycle
                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty ) cycle
                cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                      xx_ag(hx) * xx_ag(hx2) *   &
                                      fc_gf( vib2v, vib ) * fc_gf( 0, vib2 ) &
                                      * 2.d0                                !the 2 takes into account 2p-1p interactions (will be the same)
            end do
            end do
            end do
        end do
        end do
        end do
        end if
        
        !two-particle two-particle
        if ( two_state ) then
        do lx = max( 1-rx, 1 ), min( lattice_x - rx, lattice_x )
        do ly = max( 1-ry, 1 ), min( lattice_y - ry, lattice_y )
        do lz = max( 1-rz, 1 ), min( lattice_z - rz, lattice_z )
            lx2 = lx + rx
            ly2 = ly + ry
            lz2 = lz + rz
            lxyz = nx_lattice( lx, ly, lz )
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
        
            !linker type (just like 1p-1p)
            do lxv = 1, lattice_x
            do lyv = 1, lattice_y
            do lzv = 1, lattice_z
                lxyzv = nx_lattice( lxv, lyv, lzv )
                lxyz2v = lxyzv
            do vibv = 1, vibmax
               vib2v = vibv
               if ( lxyz == lxyz2 ) then    !vib must be the same or else fc factor goes to zero
                    do vib = 0, vibmax
                        vib2 = vib
                        hx = nx_2p( lxyz, vib, lxyzv, vibv )
                        if ( hx == empty ) cycle
                        hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty ) cycle
                        cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                              xx_ag(hx) * xx_ag(hx2)
                    end do
                else                        !vib can be different, now we have fc factors
                    do vib = 0, vibmax
                    do vib2 = 0, vibmax
                        hx = nx_2p( lxyz, vib, lxyzv, vibv )
                        if ( hx == empty ) cycle
                        hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty ) cycle
                        cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                 &
                                              xx_ag(hx) * xx_ag(hx2) *   &
                                              fc_gf( 0, vib ) * fc_gf( 0, vib2 )
                    end do
                    end do
                end if
            end do
            end do
            end do
            end do
            
            ! This is the exchange type in the coherence function
            if ( lxyz == lxyz2 ) cycle    !can't be exhange type
            lxyzv = lxyz2
            lxyz2v = lxyz
            do vib = 0, vibmax
            do vib2 = 0, vibmax
            do vibv = 1, vibmax
            do vib2v = 1, vibmax
                    hx = nx_2p( lxyz, vib, lxyzv, vibv )
                    if ( hx == empty ) cycle
                    hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                    if ( hx2 == empty ) cycle
                    cohfxn( rx,ry, rz ) = cohfxn( rx, ry, rz ) +                     &
                                          xx_ag(hx) * xx_ag(hx2) *       &
                                          fc_gf( vib2v, vib ) * fc_gf( vibv, vib2 )
            end do
            end do
            end do
            end do
        end do
        end do
        end do
        end if
        
    end do
    end do
    end do

end subroutine
!***************************************************************************************!
!    Write the coherence function
!***************************************************************************************!
subroutine cohfxn_out( cohfxn )
     use common_variables
   implicit none 

    character(256) f_name
    integer f_no, rx, ry, rz  
    real*8 cohfxn( -lattice_x+1:lattice_x-1, &
                   -lattice_y+1:lattice_y-1, &
                   -lattice_z+1:lattice_z-1 )
    
    !File info
    f_name = trim(task_title)//'_coh.dat'
    f_no = 2
        
    !Open the file and put headers
    open( unit = f_no, file = f_name )
    write(f_no, *) 'Task Title: ', trim(task_title)
    write(f_no, *) 'Polaron Coherence Function'
    write(f_no, *) 'x\y ',(ry, ry = -lattice_y + 1, lattice_y - 1 )
    
    rz = 0
    do rx = -lattice_x+1, lattice_x-1
        write( f_no, * ) rx, ( cohfxn( rx, ry, rz ) , ry = -lattice_y + 1, lattice_y - 1 )
    end do
    
    close( f_no )

 End Subroutine


!**********************************************!
!                        index the lattice                                           !
!**********************************************!
subroutine index_lattice()
    use common_variables
    implicit none 

    integer lx, ly, lz
    
    if ( .not. allocated( nx_lattice ) ) then 
        allocate( nx_lattice( lattice_x, lattice_y, lattice_z ) )
    end if
    nx_lattice = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
        kount_lattice = kount_lattice + 1
        nx_lattice( lx, ly, lz ) = kount_lattice
    end do
    end do
    end do
        
end subroutine

!**********************************************!
!                        index the 1 p states                                   !
!**********************************************!
subroutine index_1p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib

    !allocate the index
    if ( .not. allocated( nx_1p ) ) then 
            allocate( nx_1p ( kount_lattice, 0:vibmax ) )
    end if
    nx_1p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        kount = kount + 1
        kount_1p = kount_1p + 1
        lxyz = nx_lattice( lx, ly, lz )
        nx_1p( lxyz, vib ) = kount
        write (9,*) 'nx_1p(', lxyz, vib,' )', kount
!                print*, lxyz, vib, kount
    end do
    end do
    end do
    end do
        
end subroutine
!**********************************************!
!                        index the 2 p states                                   !
!**********************************************!

subroutine index_2p()
    use common_variables
    implicit none 
    
    integer lx, ly, lz, lxyz, vib
    integer lxv, lyv, lzv, lxyzv, vibv

    !allocate the index
    if ( .not. allocated( nx_2p ) ) then 
            allocate( nx_2p ( kount_lattice, 0:vibmax, kount_lattice, 1:vibmax ) )
    end if
    nx_2p = empty
    
    do lx = 1, lattice_x
    do ly = 1, lattice_y
    do lz = 1, lattice_z
    do vib = 0, vibmax
        do lxv = 1, lattice_x
        do lyv = 1, lattice_y
        do lzv = 1, lattice_z
!                do vib = 0, vibmax
        do vibv = 1, vibmax
            lxyz = nx_lattice( lx, ly, lz )
            lxyzv = nx_lattice( lxv, lyv, lzv ) 
            if ( lxyz == lxyzv ) cycle                !not a 2 p state
            if (vibv + vib > vibmax) cycle
            if ( two_truncate .and. two_range*1.d0 < dsqrt( 1.d0*((lx-lxv)**2+(ly-lyv)**2 + (lzv-lzv)**2)) ) cycle      !truncate at two_range
            kount = kount + 1
            kount_2p = kount_2p + 1
            nx_2p( lxyz, vib, lxyzv, vibv ) = kount
            write (9,*) 'nx_2p(', lxyz, vib, lxyzv, vibv,' )', kount
            !write (7,*) 'nx_2p(',lxyz,vib,lxyzv,vibv,')', kount
!                        print*, lxyz, vib,lxyzv, vibv, kount
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
        
end subroutine

!**********************************************************************************!
!                        create 1p sparse hamiltonian                                 !
!**********************************************************************************!
subroutine hamiltonian_1p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1
    integer lx2, ly2, lz2, lxyz2, vib2, hx2,config
    real*8 get_j,coupling,val,dis_order
    real*8 SS_disorder
    if ( .not. allocated( vibsum_state ) ) allocate( vibsum_state(kount) )
    vibsum_state = 0.d0
 call MatGetOwnershipRange(A,Istart,Iend,ierr)
    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
    lxyz1 = nx_lattice( lx1, ly1, lz1 )
      
     
      
       hx1 = nx_1p( lxyz1, vib1 )
       if ( hx1 == empty ) cycle

       hx1 = hx1-1
       if (hx1 >= Istart .and. hx1<Iend) then
       
        !diagonal
        val = vib1*1.d0
      !  if (diagonal) then
 !print*, 'vib1', vib1
 !print*, 'disorder (', config,lx1,ly1,')', disorder(config,lx1,ly1)
  !      SS_disorder = disorder(config,lx1,ly1)
 !print*, 'SS_disorder', SS_disorder
   !     dis_order = val + SS_disorder
!end if
      if (val .eq. 0) then
continue
else
 call MatSetValues ( A, 1, hx1, 1, hx1, val,INSERT_VALUES,ierr)        !Set the values
  end if
     
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
        lxyz2 = nx_lattice( lx2, ly2, lz2 )
        hx2 = nx_1p( lxyz2, vib2 )
            if ( hx2 == empty ) cycle
            hx2=hx2-1
            if ( hx2 == hx1 ) cycle
     
  
   
            !off diagonal
     coupling = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                            fc_gf( 0, vib1 ) * fc_gf( 0, vib2 )
 if (coupling .eq. 0)  then
continue
else
     call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
    end if         
        


       
        end do
        end do
        end do
        end do
end if
  
    end do
    end do
    end do
    end do

end subroutine
  
!**********************************************************************************!
!                        create 1p2p sparse hamiltonian                            !
!**********************************************************************************!

subroutine hamiltonian_1p2p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lx2v, ly2v, lz2v, lxyz2v, vib2v
    real*8 get_j,val,coupling
    
    call MatGetOwnershipRange(A,Istart,Iend,ierr)
    
    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
        hx1 = nx_1p( lxyz1, vib1 )
        if ( hx1 == empty ) cycle
        hx1 = hx1-1
        
        if (hx1 >= Istart .and. hx1<Iend) then  
      
        do lx2 = 1, lattice_x
        do ly2 = 1, lattice_y
        do lz2 = 1, lattice_z
        do vib2 = 0, vibmax
            lxyz2 = nx_lattice( lx2, ly2, lz2 )
            
            !these are the only ones that are non-zero
            lx2v = lx1 
            ly2v = ly1
            lz2v = lz1
            lxyz2v = nx_lattice( lx2v, ly2v, lz2v )
            do vib2v = 1, vibmax

                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                if ( hx2 == empty ) cycle
                hx2=hx2-1  
                
              coupling = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                                fc_gf( vib2v, vib1 ) * fc_gf( 0, vib2 )
           if (coupling .eq. 0)  then
             continue
             else
     call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
     call MatSetValues ( A, 1,hx2, 1, hx1, coupling ,INSERT_VALUES,ierr)
            end if   
               
               
           
            end do
        end do
        end do
        end do
        end do
    
    end if

    end do
    end do
    end do
    end do
   
end subroutine


!**********************************************************************************!
!                        create 2p sparse hamiltonian                              !
!**********************************************************************************!


subroutine hamiltonian_2p_sparse(config)
    use common_variables
    implicit none
    
    integer lx1, ly1, lz1, lxyz1, vib1, hx1,config
    integer lx1v, ly1v, lz1v, lxyz1v, vib1v
    integer lx2, ly2, lz2, lxyz2, vib2, hx2
    integer lxyz2v, vib2v
    real*8 get_j,coupling,val
    real*8 SS_disorder,dis_order
    call MatGetOwnershipRange(A,Istart,Iend,ierr) 

    do lx1 = 1, lattice_x
    do ly1 = 1, lattice_y
    do lz1 = 1, lattice_z
    do vib1 = 0, vibmax
        lxyz1 = nx_lattice( lx1, ly1, lz1 )
    
        do lx1v = 1, lattice_x
        do ly1v = 1, lattice_y
        do lz1v = 1, lattice_z
        do vib1v = 1, vibmax
            lxyz1v = nx_lattice( lx1v, ly1v, lz1v )
    
            
            hx1 = nx_2p( lxyz1, vib1, lxyz1v, vib1v )
            
           
            if ( hx1 == empty ) cycle
            hx1=hx1-1
            if (hx1 >= Istart .and. hx1<Iend) then
            val = (vib1 +vib1v)*1.d0
          ! if (diagonal) then
       ! SS_disorder = disorder(config,lx1,ly1)
! print*, 'SS_disorder', SS_disorder
       ! dis_order = val + SS_disorder
!end if
            if (val .eq. 0) then
            continue
            else
            call MatSetValues ( A, 1, hx1, 1, hx1, val,INSERT_VALUES,ierr)
            end if         
         
            do lx2 = 1, lattice_x
            do ly2 = 1, lattice_y
            do lz2 = 1, lattice_z
                lxyz2 = nx_lattice( lx2, ly2, lz2 )
            do vib2 = 0, vibmax
            
                !linker type
                lxyz2v = lxyz1v
                vib2v = vib1v
                hx2 = nx_2p( lxyz2, vib2, lxyz2v, vib2v )
                hx2=hx2-1
                if ( hx2 == empty .or. hx2 == hx1 ) then
                    continue
                else
               coupling = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config ) * &
                                    fc_gf( 0, vib1 ) * fc_gf( 0, vib2 )
              if (coupling .eq. 0)  then
              continue
              else
              call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
              end if
                   
                end if
        
                !exchange type
                if ( lxyz1v == lxyz2 ) then
                    lxyz2v = lxyz1
                    do vib2v = 1, vibmax
                        hx2 = nx_2p(  lxyz2, vib2, lxyz2v, vib2v )
                        if ( hx2 == empty) cycle
                        hx2=hx2-1
                        if (hx2 == hx1 ) cycle
                        
                        coupling = get_j( lx1, lx2, ly1, ly2, lz1, lz2,config) * &
                                        fc_gf( vib2v, vib1 ) * fc_gf( vib1v, vib2 )
                  if (coupling .eq. 0)  then
                      continue
                          else
                   call MatSetValues ( A, 1,hx1, 1, hx2, coupling ,INSERT_VALUES,ierr)
                      end if   
                     
                    end do
               end if
            end do
            end do
            end do
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
             
end subroutine

!***************************************************
!***************************************************
!***************************************************
                !  SS: short range both sides
!***************************************************
!***************************************************
!***************************************************
 subroutine short_short(config)
 use common_variables
 implicit none

integer vx,vx1,vy,vy1,vz,vz1,v1,v2,vxyz,vxyz1,hx1,hx2,config
real*8 val

call MatGetOwnershipRange(A,Istart,Iend,ierr) 
do vx = 1, lattice_x    
        do vy = 1, lattice_y
        do vz = 1, lattice_z
           do v1=0,vibmax  
          vxyz = nx_lattice( vx, vy, vz )
            hx1= nx_1p(vxyz,v1) 
            if (hx1==empty) cycle
            hx1=hx1-1
         if (hx1 >= Istart .and. hx1<Iend) then 
            val = v1*1.d0
            if (diagonal) then 
            val = val + disorder(config,vx,vy)
            end if
            if (dopant) then
            val = val + Force_anion(config,vx,vy)
            end if
        call MatSetValues ( A, 1,hx1, 1, hx1, val,INSERT_VALUES,ierr)   
       end if
if (two_state) then
 do vx1=1,lattice_x
            do vy1=1,lattice_y
            do vz1=1,lattice_z
            do v2=1,vibmax
             vxyz1=nx_lattice(vx1,vy1,vz1)
          hx1=nx_2p(vxyz,v1,vxyz1,v2)
          if (hx1==empty) cycle
          hx1= hx1-1
        if (hx1 >= Istart .and. hx1<Iend) then  
          val = (v1 +v2)*1.d0 
          if (diagonal) then
          val = val + disorder(config,vx,vy)
          end if
          if (dopant) then
          val = val + Force_anion(config,vx,vy)
          end if
          call MatSetValues ( A, 1,hx1, 1, hx1, val,INSERT_VALUES,ierr)    

         end if 
       
        end do
        end do
        end do
        end do
        end if

        end do
        end do
        end do
        end do     
      !  end if

       
       
 end subroutine
!**********************************************!
!                        get the coupling function                           !
!**********************************************!
real*8 function get_j( lx1, lx2, ly1, ly2, lz1, lz2,config )
    use common_variables
    implicit none
    
     
    integer, intent(in) :: lx1, lx2, ly1, ly2, lz1, lz2,config
    integer dx, dy, dz
    real*8 delta_jo_x
      real*8 z1,z2,M,theta,v1,v2
       integer, allocatable :: seed(:)
            integer :: p = 1
             call random_seed(size = p)
           allocate(seed(p))


    dx = abs( lx1 - lx2 )
    dy = abs( ly1 - ly2 )
    dz = abs( lz1 - lz2 )
    if (periodic) then
        dx = min( dx, lattice_x - dx )
        dy = min( dy, lattice_y - dy )
        dz = min( dz, lattice_z - dz )
    end if
    
    call RANDOM_NUMBER(v1)
      call RANDOM_NUMBER(v2)
        M = dsqrt(-2*log(v1))
        theta = 2*pi*v2
        z1 = M*cos(theta)
        z2 = M*sin(theta)
    
    !initialize
    get_j = 0.d0

    if ( .not. extended_cpl ) then
        if ( dy == 0 .and. dz == 0 ) then
  
     if (dx == 1  ) get_j = jo_x
  !  if (ly1==1 .and. lx1==4 .and. ly2==1 .and.  lx2==5) get_j =0
  !  if (ly1==1 .and. lx1==8 .and. ly2==1 .and.  lx2==9) get_j =0
  !   if (ly1==1 .and. lx1==15 .and. ly2==1 .and.  lx2==16) get_j =0
  !   if (ly1==1 .and. lx1==16 .and. ly2==1 .and.  lx2==17) get_j =0 
    return
        end if
      

  if ( dx == 0 .and. dz == 0 ) then
if (dy == 1 ) get_j = jo_y
!if (ly1==1 .and. ly2==2) get_j =0
!if (ly1==2 .and. ly2==1) get_j =0
return

 
!              if ( dy == 1 .and. offdiagonal ) get_j = jo_y + ((jo_y*exp(-beta*pdisorder(config,lx1,ly1,lx1,ly2)))-jo_y)
 !         return
         end if

        if ( dx == 0 .and. dy == 0 ) then
                if ( dz == 1 ) get_j = jo_z
                return
        end if
    else if ( extended_cpl ) then
        get_j = jo_ex( dx, dy, dz )
    end if
    
    return

end function
subroutine set_fctable()
        use common_variables
        implicit none

        integer vib_g, vib_exc, vib_ion, vib_n
        ! vib_g is for ground state, vib_exc for excited, vib_ion is for charged mol, and vib_n is for neutral mol.
        real*8 fc, s_diff

        !allocate matrices
        if ( .not. allocated( fc_gf ) ) then 
                allocate( fc_gf ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_gc ) ) then 
                allocate( fc_gc ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_ga ) ) then 
                allocate( fc_ga ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_cf ) ) then 
                allocate( fc_cf ( 0:vibmax, 0:vibmax ) )
        end if
        if ( .not. allocated( fc_af ) ) then 
                allocate( fc_af ( 0:vibmax, 0:vibmax ) )
        end if

        !----- generate fc table ----!
        do vib_g =0,vibmax
        do vib_exc =0,vibmax
                call fcfac( vib_g,vib_exc,s,fc)
                fc_gf(vib_g, vib_exc) = fc
        enddo
        enddo

        !ground to cation
        do vib_n  =0,vibmax
        do vib_ion=0,vibmax
                call fcfac(vib_n, vib_ion,s_cat,fc)
                fc_gc(vib_n, vib_ion) = fc
        enddo
        enddo
        
        !ground to anion
        do vib_n  =0,vibmax
        do vib_ion=0,vibmax
                call fcfac(vib_n, vib_ion,s_ani,fc)
                fc_ga(vib_n, vib_ion) = fc
        enddo
        enddo

        !---- generate fc table for ct ----!
        ! the procedure is designed in a way that the first potential well is inside
        ! the second. if this order is reversed, it may give wrong signe.
        ! for example, if s_frenkel > se, then the order should be se then s_frenkel.
        ! however, if s_frenkel < se, then the order becomes s_frenkel and se.
        ! 
        ! 1/27/2009
        ! up to this point, s(1) = 1 and nuclear displacement factor of
        ! ct states (lamda) are set such that the sum of the lamda is s(1), which is 1
        ! since lamda^2 is s, root sct_h + root sct_e = root s(1)=1.
        
        !cation to frenkel      
        do vib_ion=0,vibmax
        do vib_n  =0,vibmax
                s_diff = dsqrt(s) - dsqrt(s_cat)
                s_diff = s_diff**2
                if( s <= s_cat ) then
                        call fcfac(vib_n,vib_ion,dabs( s_diff ),fc)
                else
                        call fcfac(vib_ion,vib_n,dabs( s_diff ),fc)
                endif
                fc_cf(vib_ion, vib_n) = fc
        enddo
        enddo
        
        !anion to frenkel
        do vib_ion=0,vibmax
        do vib_n  =0,vibmax
                s_diff = dsqrt(s) - dsqrt(s_ani)
                s_diff = s_diff**2
                if( s <= s_ani ) then
                        call fcfac(vib_n,vib_ion,dabs( s_diff ),fc)
                else
                        call fcfac(vib_ion,vib_n,dabs( s_diff ),fc)
                endif
                fc_af(vib_ion, vib_n) = fc
        enddo
        enddo

end subroutine

!********************************************************************
! franc-condon calculations     !
!********************************************************************
subroutine fcfac(n,m,s,fc)
        implicit none
        integer n,m,k
        real*8 s,fc,f_m,f_n,f_k,f_nmk,f_mk,facin
        real*8 fact

        fc = 0.d0

        do k = 0,m
                if(n-m+k < 0) go to 100 ! if n-m+k is negative, factorial is not calculatable.

                f_mk  = fact(m-k)
                f_nmk = fact(n-m+k)
                f_k   = fact(k)
                facin = 1.d0/(1.d0*f_k*f_mk*f_nmk)

                fc = fc + facin*s**(k*0.5d0)*s**(1.0d0*(n-m+k)*0.5d0)*(-1)**(n-m+k)
100             continue
        enddo

        f_n = fact(n)
        f_m = fact(m)
!       print*,'f_n,f_m=',f_n,f_m
        fc = fc*dsqrt(1.d0*f_m*f_n)*dexp(-s/2.d0)

        return
end subroutine

!********************************************************************
!       calculating factorial
!       * i think, mathimatically numbers should be integers, but
!         here, i use real numbers for convenience.
!********************************************************************
real*8 function fact( n )
        implicit none
        integer n
        integer i

        fact = 1.0d0
        if( n .ge. 0 ) then
                do i = 2, n
                        fact = fact * i
                enddo
        endif
end function

       
        subroutine disorder_elements()
        use common_variables
        implicit none
        
        integer vx,vy,config
        real*8 z1,z2,R,theta,v1,v2
      
         integer, allocatable :: seed(:)
            integer :: p = 1
             call random_seed(size = p)
           allocate(seed(p))
            call random_seed(put=seed)
         !   write (*, *) seed
        
      !  call RANDOM_SEED()
    
         if ( .not. allocated( disorder ) ) then 
            allocate( disorder ( conf_max, lattice_x, lattice_y ))
    end if
       
        do config =1, conf_max
       
        do vx=1,lattice_x
        do vy=1,lattice_y
    
        
        
        call RANDOM_NUMBER(v1)
        call RANDOM_NUMBER(v2)
        R = dsqrt(-2*log(v1))
        theta = 2*pi*v2
        z1 = R*cos(theta)
        z2 = R*sin(theta)
        
        z2 = z2 * ((sigma*8065)/1400) 
        z1 = z1 * ((sigma*8065)/1400)

      !  do vy=1,lattice_y
        disorder (config,vx,vy) = z2
!  print*, 'disorder (', config,vx,vy,')', disorder(config,vx,vy)
            end do
end do
end do

end subroutine  

subroutine doping_disorder()
 use common_variables
        implicit none
        
        integer vx,vy,config
       real dx1,dy,dx2,dx
       
        if ( .not. allocated( D_anion1 ) ) then 
            allocate( D_anion1 (conf_max, lattice_x, lattice_y ))
    end if
     if ( .not. allocated( D_anion2 ) ) then 
            allocate( D_anion2 (conf_max, lattice_x, lattice_y ))
    end if
     if ( .not. allocated( D_anion3 ) ) then 
            allocate( D_anion3 (conf_max, lattice_x, lattice_y ))
    end if
        if ( .not. allocated( Force_anion ) ) then 
            allocate( Force_anion ( conf_max, lattice_x, lattice_y ))
    end if
     if ( .not. allocated( D_anion4 ) ) then 
            allocate( D_anion4 (conf_max, lattice_x, lattice_y ))
    end if
    
      if ( .not. allocated( D_cation ) ) then 
            allocate( D_cation (conf_max, lattice_x, lattice_y ))
    end if
    
        if ( .not. allocated( Force_cation ) ) then 
            allocate( Force_cation ( conf_max, lattice_x, lattice_y ))
    end if
!    Force = 0
       do config =1, conf_max
       
      !  D(config,3,3) = r 
        do vx=1,lattice_x
        do vy=1,lattice_y
        
    
      dx = abs(vx-5.5)
       dx1 = abs(vx-3)
       dy = abs(vy-2.5)
        dx2 = abs(vx-6.125)                    



!D_anion1(config, vx,vy) = sqrt( ((vx-1)*0.4-1.8)**2  + ((vy-1)*0.4-0.4)**2 + (0.7)**2)

   

D_anion1(config, vx,vy) = sqrt( ((vx-1)*0.4-pos_anion)**2  + ((vy-1)*0.4+danion)**2)

Force_anion(config,vx,vy) = (-8.3/1)/(D_anion1(config,vx,vy))


    
        end do
        end do
        end do
        end subroutine
       
  subroutine paracrystalline_SS()
        use common_variables
        implicit none
        
        integer lx1,lx2, ly1,ly2
        real *8 rand,g,v1,v2,z1,z2,R,theta
        integer config
     
          
         integer, allocatable :: seed(:)
            integer :: p = 1
             call random_seed(size = p)
           allocate(seed(p))
            call random_seed(put=seed)


      if ( .not. allocated( pdisorder ) ) then 
            allocate( pdisorder ( conf_max,lattice_x,lattice_y, lattice_x, lattice_y))
      end if
      rand = 0

  do config=1, conf_max
  do lx1=1, lattice_x
  ! lx2 = lx1
  do ly1=1, lattice_y-1
  ly2=ly1+1
        
        call RANDOM_NUMBER(v1)
        call RANDOM_NUMBER(v2)
        R = dsqrt(-2*log(v1))
        theta = 2*pi*v2
        z1 = R*cos(theta)
        z2 = R*sin(theta)
        z2 = z2 * 0.19
        g = 0.19/3.8


! pdisorder (config,lx1,1,lx1,2) = 0.1
! pdisorder (config,lx1,2,lx1,1) = 0.1
! pdisorder (config,lx1,2,lx1,3) = 0
! pdisorder (config,lx1,3,lx1,2) = 0
! pdisorder (config,lx1,3,lx1,4) = -0.1
! pdisorder (config,lx1,4,lx1,3) = -0.1

      
 pdisorder (config,lx1,ly1,lx1,ly2) = z2
 pdisorder (config,lx1,ly2,lx1,ly1) = z2
 
 ! print*, 'pdisorder (', config,lx1,ly1,lx1, ly2, ')', pdisorder (config,lx1,ly1,lx1, ly2)
 ! print*, 'pdisorder (', config,lx1,ly2,lx1, ly1, ')', pdisorder (config,lx1,ly2,lx1, ly1)

   
    
   end do
   end do
   end do 
      
  end subroutine

subroutine para_range(n1, n2, communicators, color, ista, ien)
integer(4) :: n1 ! Lowest value of iteration variable
integer(4) :: n2 ! Highest value of iteration variable
integer(4) :: communicators ! # of communicators
integer(4) :: color ! # cores
integer(4) :: row_rank
integer(4) :: ista ! Start of iterations for ith communicator 
integer(4) :: ien ! End of iterations for ith communicator 
ista = n1+(color*n2/communicators)
ien = ((color+1)*n2)/communicators
return
end subroutine para_range
        
subroutine density_of_groundstates(bn1,bn2,bn3,bn4,bn5)
  use common_variables
        implicit none
       
        
        integer bn1,bn2,bn3,bn4,bn5  
       
if (row_rank == 0) then
        if (kr_g > -6 .and. kr_g < -5) then 
        bn1=bn1+1
        end if
        if (kr_g > -5 .and. kr_g < -4) then 
        bn2=bn2+1
        end if
        if (kr_g > -4 .and. kr_g < -3) then 
        bn3=bn3+1
        end if
        if (kr_g > -3 .and. kr_g < -2) then 
        bn4=bn4+1
        end if
        if (kr_g > -2 .and. kr_g < -1) then 
        bn5=bn5+1
        end if
       
  !print*, 'bin1:', bn1
 ! print*, 'bin2:', bn2
 ! print*, 'bin3:', bn3
 ! print*, 'bin4:', bn4
 ! print*, 'bin5:', bn5
 end if
 end subroutine            
