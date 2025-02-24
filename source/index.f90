module index
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: LatticeIndex, oneParticleIndex, twoParticleIndex, chargeTransferIndex
    contains
        subroutine LatticeIndex(lattice_dimx, lattice_dimy, lattice_dimz, lattice_index_arr, lattice_count, general_counter, empty)
            implicit none
            integer(wp) :: i_x, i_y, i_z
            integer, intent(inout) :: general_counter
            integer(wp), intent(inout) ::  lattice_count
            integer(wp), intent(in) :: lattice_dimx, lattice_dimy, lattice_dimz, empty
            integer(wp), dimension(:,:,:), intent(inout) :: lattice_index_arr
        
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
        
        
        subroutine oneParticleIndex(one_particle_index_arr,lattice_dimx,lattice_dimy,lattice_dimz, max_vibs,one_particle_counter, empty, general_counter, lattice_index_arr)
            ! |n,ṽ>
            implicit none
            integer(wp), dimension(:,0:) :: one_particle_index_arr
            integer(wp), dimension(:,:,:), intent(in) :: lattice_index_arr
            integer(wp), intent(inout) :: one_particle_counter
            integer, intent(inout) :: general_counter
            integer(wp), intent(in) :: max_vibs, lattice_dimx, lattice_dimy, lattice_dimz, empty
            integer(wp) :: i_x, i_y, i_z, vib, indx_xyz
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
        
        subroutine twoParticleIndex(two_particle_index_arr,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,two_particle_counter,empty, general_counter,lattice_index_arr)
            ! |n,ṽ;n',v'>
            ! index two-particle states; two particle states consist of a vibronic excitation on one site and a vibrational excitation on another
            ! therefore each combination of site with vibronic exc and vibrational exc on other site need an index number
            implicit none
            integer(wp) i_x, i_y, i_z, i_xyz, vib ! indices for vibronic excitations
            integer(wp) i_xv, i_yv, i_zv, i_xyzv, vibv ! indices for vibrational excitations
            integer(wp), intent(in) :: max_vibs, lattice_dimx, lattice_dimy, lattice_dimz, empty
            integer(wp), dimension(:,:,:), intent(in) :: lattice_index_arr
            integer(wp), intent(inout) :: two_particle_index_arr(:,0:,:,1:)
            integer(wp), intent(inout) :: two_particle_counter
            integer, intent(inout) :: general_counter


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
            ! print*, general_counter, 'Total basis states'
            ! print*, '*****************************************'                 
        end subroutine



        subroutine chargeTransferIndex(chargetransfer_index_arr,lattice_dimx,lattice_dimy,lattice_dimz,max_vibs,chargetransfer_counter,empty, general_counter,lattice_index_arr,bonding_matrix)
            ! |n,v+;n+s,v->
            !LM index charge transfer states consisting of a vibrationally excited cation and an anion.
            !LM charge transfer states are tricky. Spano and Hestand use the index s for one-dimensional aggregates.
            !LM I think I can still use s when I unravel the whole aggregate into a 1d array with lattice_index_arr
            
            !LM 2025: restrict charge transfer sites to bonded sites with a nearest neighbour separation of 1 only.

            implicit none
            integer(wp) i_xc, i_yc, i_zc, i_xyzc, vibc ! indices for cations/ vibrations
            integer(wp) i_xa, i_ya, i_za, i_xyza, viba ! indices for anions/ vibrations
            integer(wp), intent(in) :: max_vibs, lattice_dimx, lattice_dimy, lattice_dimz, empty
            integer(wp), dimension(:,:,:), intent(in) :: lattice_index_arr
            integer(wp), intent(in) :: bonding_matrix(:,:)
            integer(wp), intent(inout) :: chargetransfer_index_arr(:,0:,:,0:)
            integer(wp), intent(inout) :: chargetransfer_counter
            integer, intent(inout) :: general_counter
            chargetransfer_index_arr = empty
            do i_xc=1,lattice_dimx
                do i_yc=1,lattice_dimy
                    do i_zc=1,lattice_dimz
                        do vibc=0,max_vibs ! for each cation state, loop over all other lattice sites for anion sites
                            do i_xa=1,lattice_dimx
                                do i_ya=1,lattice_dimy
                                    do i_za=1,lattice_dimz
                                        do viba=0,max_vibs ! excited vibrational states in ground electronic states cannot be v=0
                                            i_xyzc= lattice_index_arr(i_xc,i_yc,i_zc) ! n
                                            i_xyza=lattice_index_arr(i_xa,i_ya,i_za) ! n+s
                                            if ( i_xyzc .eq. i_xyza ) cycle ! |s|>=1 (s=0 would be Frenkel excitons)
                                            if (viba + vibc > max_vibs) cycle 
                                            if (bonding_matrix(i_xyza,i_xyzc) .ne. 1) cycle
                                            general_counter = general_counter + 1
                                            chargetransfer_counter = chargetransfer_counter + 1
                                            chargetransfer_index_arr(i_xyzc, vibc, i_xyza, viba) = general_counter
                                        end do
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
            print*, chargetransfer_counter, 'Charge-transfer states'
            ! print*, general_counter, 'Total basis states'
            ! print*, '*****************************************'      
        end subroutine

end module