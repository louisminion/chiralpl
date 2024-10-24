program ind1p
    integer, allocatable :: one_particle_index_arr(:,:)
    integer max_vibs
    integer lattice_count
    read(*,*) lattice_count
    read(*,*) max_vibs
    if ( .not. allocated( one_particle_index_arr ) ) then 
        allocate( one_particle_index_arr( lattice_count, 1:max_vibs) )
    end if
    write(*,*) size(one_particle_index_arr,dim=2)
end program