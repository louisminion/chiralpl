module geometry
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: construct_geometry
    contains
        subroutine construct_geometry(r_xyz,GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz,x_spacing,y_spacing,z_spacing,phi,theta,period_Jtwist,lattice_count,lattice_index_arr)
            integer(wp), intent(in) :: GEOMETRY_TYPE,lattice_dimx,lattice_dimy,lattice_dimz, period_Jtwist,lattice_count
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
            write(*,*) midpointx, 'GHHH'
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
end module