module index
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    private
    public :: LatticeIndex, oneParticleIndex, twoParticleIndex
    contains
        subroutine LatticeIndex()
            use variables
            implicit none
            integer(wp) :: i_x, i_y, i_z
        
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
        
        
        subroutine oneParticleIndex()
            use variables
            implicit none
            integer(wp) :: i_x, i_y, i_z, vib, indx_xyz
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
            print*, one_particle_counter, 'One particle states'
        end subroutine
        
        subroutine twoParticleIndex()
            use variables
            ! index two-particle states; two particle states consist of a vibronic excitation on one site and a vibrational excitation on another
            ! therefore each combination of site with vibronic exc and vibrational exc on other site need an index number
            implicit none
            integer(wp) i_x, i_y, i_z, i_xyz, vib ! indices for vibronic excitations
            integer(wp) i_xv, i_yv, i_zv, i_xyzv, vibv ! indices for vibrational excitations
        
        
            if (.not. allocated(two_particle_index_arr)) then
                allocate( two_particle_index_arr(lattice_count, 0:max_vibs, lattice_count, 1:max_vibs))
            
            end if
            two_particle_index_arr = empty
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

end module