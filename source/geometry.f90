module geometry
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: construct_geometry, construct_bonding_matrix
    contains
        subroutine construct_geometry(r_xyz,GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,x_spacing,y_spacing,z_spacing,phi,theta,lattice_count,lattice_index_arr)
            integer(wp), intent(in) :: GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,lattice_count
            real(wp), intent(in) :: phi, theta,x_spacing,y_spacing,z_spacing
            real(wp), allocatable, intent(inout) :: r_xyz(:,:)
            integer(wp), dimension(:,:,:), intent(in) :: lattice_index_arr
            integer(wp) :: ix,iy,iz,ixyz
            real(wp) :: sum_xpos,sum_ypos,counter
            real(wp) :: midpointx,midpointy
            allocate(r_xyz(lattice_count,3))
            r_xyz = 0.0_wp

            midpointx = (lattice_dimx*1.0_wp+1.0_wp)/2.0_wp
            midpointy = (lattice_dimy*1.0_wp+1.0_wp)/2.0_wp
            ! write(*,*) midpointx, 'GHHH'
            ! if (lattice_dimy .eq. 1) then
            !     midpointy=0
            ! end if
            if (GEOMETRY_TYPE .eq. 0) then
                do ix=1,lattice_dimx
                    do iy=1,lattice_dimy
                        do iz=1,lattice_dimz
                            ixyz = lattice_index_arr(ix,iy,iz)
                            r_xyz(ixyz,1) = ((ix*1.0_wp)-midpointx)*x_spacing*cos((((iz*1.0_wp)-1.0_wp))*phi)-((iy*1.0_wp)-midpointy)*y_spacing*sin(((iz*1.0_wp)-1.0_wp)*phi)
                            r_xyz(ixyz,2) = ((ix*1.0_wp)-midpointx)*x_spacing*sin(((iz*1.0_wp)-1.0_wp)*phi)+((iy*1.0_wp)-midpointy)*y_spacing*cos(((iz*1.0_wp)-1.0_wp)*phi)  !(iy*1.0_wp)*y_spacing
                            r_xyz(ixyz,3) = ((iz*1.0_wp)-1.0_wp)*z_spacing
                        end do
                    end do
                end do

            else if (GEOMETRY_TYPE .eq. 1) then
                do ix=1,lattice_dimx
                    do iy=1,lattice_dimy
                        do iz=1,lattice_dimz
                            ixyz = lattice_index_arr(ix,iy,iz)
                            r_xyz(ixyz,1) = ((ix*1.0_wp)-midpointx)*x_spacing*cos((((iz*1.0_wp)-1.0_wp))*phi)-((iy*1.0_wp)-midpointy)*y_spacing*sin(((iz*1.0_wp)-1.0_wp)*phi)
                            r_xyz(ixyz,2) = ((ix*1.0_wp)-midpointx)*x_spacing*sin(((iz*1.0_wp)-1.0_wp)*phi)+((iy*1.0_wp)-midpointy)*y_spacing*cos(((iz*1.0_wp)-1.0_wp)*phi)  !(iy*1.0_wp)*y_spacing
                            r_xyz(ixyz,3) = ((iz*1.0_wp)-1.0_wp)*z_spacing
                        end do
                    end do
                end do

            else if (GEOMETRY_TYPE .eq. 2) then ! helically coiled polymer, x dimension is the polymer length.
                counter = 0.0_wp
                sum_xpos = 0.0_wp
                sum_ypos = 0.0_wp
                do ix=1,lattice_dimx
                    ixyz = lattice_index_arr(ix,1,1)
                    sum_xpos = sum_xpos + x_spacing*cos(((ix*1.0_wp)-2.0_wp)*phi)
                    sum_ypos = sum_ypos + x_spacing*sin(((ix*1.0_wp)-2.0_wp)*phi)
                    r_xyz(ixyz,1)= sum_xpos
                    r_xyz(ixyz,2)= sum_ypos
                    r_xyz(ixyz,3) = x_spacing*((ix*1.0_wp)-1.0_wp)
                end do

            end if 
        end subroutine

        subroutine construct_bonding_matrix(bonding_matrix,GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,lattice_count,lattice_index_arr)
            ! The bonding matrix tells us which sites are bonded together by having a 1 in the matrix, versus which are non-bonded (having a 0)
            ! Used to assign different values of the transfer integrals t_e and t_h depending on whether sites are bonded.
            integer(wp), intent(in) :: GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,lattice_count
            integer(wp), allocatable, intent(inout) :: bonding_matrix(:,:)
            integer(wp), dimension(:,:,:), intent(in) :: lattice_index_arr
            integer(wp) :: ix,iy,iz,ixyz   
            integer(wp) :: ix1,iy1,iz1,ixyz1         
            allocate(bonding_matrix(lattice_count,lattice_count))
            bonding_matrix = 0
            if (GEOMETRY_TYPE .eq. 0 .or. GEOMETRY_TYPE .eq. 1) then
                do ix=1,lattice_dimx
                    do iy=1,lattice_dimy
                        do iz=1,lattice_dimz
                            ixyz = lattice_index_arr(ix,iy,iz)
                            do ix1=1,lattice_dimx
                                do iy1=1,lattice_dimy
                                    do iz1=1,lattice_dimz
                                        ixyz1 = lattice_index_arr(ix1,iy1,iz1)
                                        if (iz .eq. iz1) then
                                            if (abs(ix-ix1) .eq. 1) then
                                                bonding_matrix(ixyz,ixyz1) = 1
                                            end if
                                        else
                                            cycle
                                            
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            else if (GEOMETRY_TYPE .eq. 2) then
                do ix=1,lattice_dimx
                    do iy=1,lattice_dimy
                        do iz=1,lattice_dimz
                            ixyz = lattice_index_arr(ix,iy,iz)
                            do ix1=1,lattice_dimx
                                do iy1=1,lattice_dimy
                                    do iz1=1,lattice_dimz
                                        ixyz1 = lattice_index_arr(ix1,iy1,iz1)
                                        if (abs(ix-ix1) .eq. 1) then
                                            ! For coiled helices we are only dealing with one of them so x-1 and x+1 are bonded
                                            bonding_matrix(ixyz,ixyz1) = 1
                                        else
                                            cycle
                                        end if
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end if
        end subroutine
end module