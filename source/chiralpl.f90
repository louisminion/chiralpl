program chiralpl
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use input
    use output
    use variables
    use index
    use franckcondon
    use geometry
    use hamiltonian
    use disorder
    use spectra
    use omp_lib
    implicit none
    integer(wp) :: o, threads, thread, vt,x
    real(wp) :: start_time, end_time
    call printOutputHeader()
    call readInput()

    allocate( lattice_index_arr( lattice_dimx, lattice_dimy, lattice_dimz) )
    call LatticeIndex(lattice_dimx, lattice_dimy, lattice_dimz, lattice_index_arr, lattice_count, general_counter, empty)
    !implement accurate geometry
    call construct_geometry(r_xyz,GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,x_spacing,y_spacing,z_spacing,phi,theta,lattice_count,lattice_index_arr)
    call construct_bonding_matrix(bonding_matrix,GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,lattice_count,lattice_index_arr)

    if ( bool_one_particle_states ) then
        allocate( one_particle_index_arr( lattice_count, 0:max_vibs) )
        call oneParticleIndex(one_particle_index_arr,lattice_dimx,lattice_dimy,lattice_dimz, max_vibs,one_particle_counter, empty, general_counter, lattice_index_arr)
    end if
    if ( bool_two_particle_states ) then
        allocate( two_particle_index_arr(lattice_count, 0:max_vibs, lattice_count, 1:max_vibs))
        call twoParticleIndex(two_particle_index_arr,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,two_particle_counter,empty, general_counter,lattice_index_arr)
    end if
    if ( bool_charge_transfer_states ) then
        allocate( chargetransfer_index_arr(lattice_count, 0:max_vibs, lattice_count, 0:max_vibs))
        call chargeTransferIndex(chargetransfer_index_arr,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,chargetransfer_counter,empty, general_counter,lattice_index_arr,bonding_matrix)
    end if
    print*, general_counter, 'Total basis states'
    print*, '*****************************************'   
    print*, 'Precalculating vibrational overlap integrals'
    ! These can be done using global variables as don't need them in parallel environment
    call calcFranckCondonTables()
    ! open(unit=90, file='RARR')
    ! do x = 1, size(bonding_matrix,dim=1)
    !     write(90, '(*(I0 : ", "))') bonding_matrix(x,:)
    ! end do
    ! close(90)
    call writePositionArray(INPUT_NAME,lattice_dimx,lattice_dimy,lattice_dimz,r_xyz)
    ! open(unit=90, file='RARR')
    ! do x = 1, size(r_xyz,dim=1)
    !     write(90, '(*(F12.8 : ", "))') r_xyz(x, 1),r_xyz(x, 2),r_xyz(x, 3)
    ! end do
    ! close(90)
    call dipole_moment()
    call writeDipoleArray(INPUT_NAME,lattice_dimx,lattice_dimy,lattice_dimz,mu_xyz)

    print*,'Hamiltonian size()','(',general_counter,general_counter,')'
    call checkRAM(general_counter)
    call construct_covariance_matrix()
    call cholesky_decomp(A_covar,lattice_count)
    call omp_set_num_threads(num_threads)
    CALL RANDOM_SEED 
    write(*,*) general_counter, 'general_counter'

    allocate(abs_specx_configavg(spec_steps))
    allocate(abs_specy_configavg(spec_steps))
    allocate(abs_specz_configavg(spec_steps))

    abs_specx_configavg = 0.0_wp
    abs_specy_configavg = 0.0_wp
    abs_specz_configavg = 0.0_wp

    allocate(pl_specx_by_v(spec_steps,0:max_vibs))
    allocate(pl_specy_by_v(spec_steps,0:max_vibs))
    allocate(pl_specz_by_v(spec_steps,0:max_vibs))
    pl_specx_by_v = 0.0_wp
    pl_specy_by_v = 0.0_wp
    pl_specz_by_v = 0.0_wp
    allocate(cd_spec_avg(spec_steps))
    start_time = omp_get_wtime()
    !$OMP PARALLEL
    threads = omp_get_num_threads()
    !$OMP END PARALLEL
    write(*,'(a,I0,a)') 'Entering parallelised region with ',threads,' threads'
    !$OMP PARALLEL SHARED(A_covar,pl_specx_by_v,pl_specy_by_v,pl_specz_by_v,abs_specx_configavg,abs_specy_configavg,abs_specz_configavg, cd_spec_avg) PRIVATE(diagonal_disorder_offsets,H&
    !$OMP,EVAL,xpl_osc,ypl_osc,zpl_osc,pl_specx,pl_specy,pl_specz,pl_specx_by_v_perconfig,pl_specy_by_v_perconfig,pl_specz_by_v_perconfig,abs_osc_strengths_x,abs_osc_strengths_y,abs_osc_strengths_z,abs_specx,abs_specy,abs_specz,cd_rot_strengths,cd_spec)
    allocate(H(general_counter,general_counter))
    allocate( EVAL( general_counter ) )
    allocate(diagonal_disorder_offsets(lattice_count))
    allocate(xpl_osc(max_vibs+1,general_counter))
    allocate(ypl_osc(max_vibs+1,general_counter))
    allocate(zpl_osc(max_vibs+1,general_counter))
    allocate(abs_osc_strengths_x(general_counter))
    allocate(abs_osc_strengths_y(general_counter))
    allocate(abs_osc_strengths_z(general_counter))
    allocate(pl_specx(spec_steps))
    allocate(pl_specy(spec_steps))
    allocate(pl_specz(spec_steps))

    allocate(pl_specx_by_v_perconfig(spec_steps,0:max_vibs))
    allocate(pl_specy_by_v_perconfig(spec_steps,0:max_vibs))
    allocate(pl_specz_by_v_perconfig(spec_steps,0:max_vibs))

    allocate(abs_specx(spec_steps))
    allocate(abs_specy(spec_steps))
    allocate(abs_specz(spec_steps))
    allocate(cd_rot_strengths(general_counter))
    allocate(cd_spec(spec_steps))
    !$OMP DO PRIVATE(o) 
    do o=1,configs
        pl_specx_by_v_perconfig = 0.0_wp
        pl_specy_by_v_perconfig = 0.0_wp
        pl_specz_by_v_perconfig = 0.0_wp
        diagonal_disorder_offsets = 0.0_wp
        EVAL = 0.0_wp
        H = 0.0_wp
        call draw_multivar_distr(lattice_count,diagonal_disorder_offsets,sigma,A_covar)
        ! write(*,*) diagonal_disorder_offsets
        ! Build hamiltonian
        ! Lots of arguments! But this is mostly so we can avoid any global variables
        if ( bool_one_particle_states ) call build1particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz,r_xyz,fc_ground_to_neutral,lattice_dimx &
        ,lattice_dimy,lattice_dimz,max_vibs,one_particle_index_arr,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
        if ( bool_two_particle_states ) call build2particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz,r_xyz,fc_ground_to_neutral,lattice_dimx &
        ,lattice_dimy,lattice_dimz,max_vibs,two_particle_index_arr,manual_coupling,x_spacing,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
        if ( bool_one_particle_states .and. bool_two_particle_states ) call build1particle2particleHamiltonian(H,diagonal_disorder_offsets,mu_xyz, r_xyz &
        ,fc_ground_to_neutral,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,one_particle_index_arr,two_particle_index_arr,manual_coupling,x_spacing &
        ,y_spacing,z_spacing, lattice_index_arr,pi, epsilon,w00)
        if (bool_charge_transfer_states) call buildChargeTransferHamiltonian(H,lattice_dimx,lattice_dimy,lattice_dimz,lattice_index_arr,max_vibs &
        ,chargetransfer_index_arr,E_CT,te,th,fc_ground_to_neutral,fc_ground_to_anion,fc_ground_to_cation)
        if ( bool_one_particle_states .and. bool_charge_transfer_states) call buildOneParticleChargeTransferHamiltonian(H,lattice_dimx,lattice_dimy,lattice_dimz,lattice_index_arr,max_vibs,one_particle_index_arr &
        ,chargetransfer_index_arr,te,th,fc_ground_to_neutral,fc_ground_to_anion,fc_ground_to_cation,fc_cation_to_frenkel,fc_anion_to_frenkel)
        ! Hamiltonian diagonalization
        call Diagonalize(H,'A',general_counter, EVAL, EVAL_COUNT, IU)
        thread = omp_get_thread_num()
        call absorption(abs_osc_strengths_x,abs_osc_strengths_y,abs_osc_strengths_z,general_counter, lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,lattice_index_arr &
        ,one_particle_index_arr,mu_xyz,fc_ground_to_neutral,H)
        call pl(xpl_osc,ypl_osc,zpl_osc,general_counter,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,lattice_index_arr,one_particle_index_arr,two_particle_index_arr &
        ,mu_xyz,fc_ground_to_neutral,H, bool_two_particle_states)
        do vt=0,max_vibs
            call calc_vibpl_spec_per_config(vt,spec_steps,w00,lw,xpl_osc,ypl_osc,zpl_osc,pl_specx,pl_specy,pl_specz, EVAL,general_counter, temp,kB,pi)
            pl_specx_by_v_perconfig(:,vt) = pl_specx 
            pl_specy_by_v_perconfig(:,vt) = pl_specy 
            pl_specz_by_v_perconfig(:,vt) = pl_specz 
        end do
        call calc_abs_spec(abs_specx,abs_specy,abs_specz, EVAL,w00,lw,pi, abs_osc_strengths_x,abs_osc_strengths_y,abs_osc_strengths_z,lattice_count,spec_steps,general_counter)
        call cd(cd_rot_strengths,H,lattice_dimx,lattice_dimy,lattice_dimz,lattice_index_arr,one_particle_index_arr,general_counter,mu_xyz,max_vibs,fc_ground_to_neutral,k,mu_0,x_spacing,y_spacing,z_spacing,r_xyz)
        call calc_cd_spec(cd_rot_strengths, cd_spec, EVAL,w00,lw,pi,lattice_count,spec_steps,general_counter)
        ! DO ADDING AT THE END OF THE LOOP SO WAITING AT A MINIMUM
        !$omp critical(UPDATESPEC)
        pl_specx_by_v = pl_specx_by_v + pl_specx_by_v_perconfig
        pl_specy_by_v = pl_specy_by_v + pl_specy_by_v_perconfig
        pl_specz_by_v = pl_specz_by_v + pl_specz_by_v_perconfig
        abs_specx_configavg = abs_specx_configavg + abs_specx
        abs_specy_configavg = abs_specy_configavg + abs_specy
        abs_specz_configavg = abs_specz_configavg + abs_specz
        cd_spec_avg = cd_spec_avg + cd_spec
        !$omp end critical(UPDATESPEC)
    end do
    !$OMP END DO
    deallocate(H)
    deallocate( EVAL )
    deallocate(diagonal_disorder_offsets)
    deallocate(xpl_osc)
    deallocate(ypl_osc)
    deallocate(zpl_osc)
    deallocate(abs_osc_strengths_x)
    deallocate(abs_osc_strengths_y)
    deallocate(abs_osc_strengths_z)
    deallocate(abs_specx)
    deallocate(abs_specy)
    deallocate(abs_specz)
    deallocate(pl_specx_by_v_perconfig)
    deallocate(pl_specy_by_v_perconfig)
    deallocate(pl_specz_by_v_perconfig)
    deallocate(cd_spec)
    !$OMP END PARALLEL
    end_time = omp_get_wtime()
    write(*,'(a,F12.5,a)') 'Configurational average completed in',end_time-start_time,' seconds.'
    pl_specx_by_v = pl_specx_by_v/(configs*1.0_wp)
    pl_specy_by_v = pl_specy_by_v/(configs*1.0_wp)
    pl_specz_by_v = pl_specz_by_v/(configs*1.0_wp)
    cd_spec_avg = cd_spec_avg/(configs*1.0_wp)
    abs_specx_configavg = abs_specx_configavg/(configs*1.0_wp)
    abs_specy_configavg = abs_specy_configavg/(configs*1.0_wp)
    abs_specz_configavg = abs_specz_configavg/(configs*1.0_wp)
    ! Reallocate pl_specx etc as only local to those threads before
    allocate(pl_specx(spec_steps))
    allocate(pl_specy(spec_steps))
    allocate(pl_specz(spec_steps))

    pl_specx = 0.0_wp
    pl_specy = 0.0_wp
    pl_specz = 0.0_wp
    do vt=0,max_vibs
        pl_specx = pl_specx + pl_specx_by_v(:,vt)
        pl_specy = pl_specy + pl_specy_by_v(:,vt)
        pl_specz = pl_specz + pl_specz_by_v(:,vt)

    end do
    call write_pl(pl_specx,pl_specy,pl_specz,INPUT_NAME,spec_steps, hw, max_vibs,w00,lw)
    call write_abs(abs_specx_configavg,abs_specy_configavg,abs_specz_configavg,w00,hw,max_vibs,spec_steps,INPUT_NAME)
    call write_cd(cd_spec_avg,w00,hw,max_vibs,spec_steps,INPUT_NAME)
    ! abs_osc_strengths_x_configavg = abs_osc_strengths_x_configavg/(configs*1.0_wp)
    ! call calc_abs_spec()
    ! call write_abs()
    ! write(*,*) abs_osc_strengths_x_configavg

end program