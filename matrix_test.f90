program matrix_test
    implicit none
    real, allocatable :: mu_xyz(:,:)
    integer :: N

    read(*,*) N
    allocate(mu_xyz(N,3))
    mu_xyz = 0.D0
    mu_xyz(1,1) = 1
    mu_xyz(1,2) = 2
    mu_xyz(1,3) = 3
    write(*,*) mu_xyz(:,2)

end program