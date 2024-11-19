module random_normal_distr
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: box_muller, box_muller_vec
    contains
        function box_muller() result(Z)
            implicit none
            integer(wp) i
            real(wp) :: U(2)
            real(wp) :: Z(2)
            real(wp),parameter :: pi=4.0_wp*datan(1.0_wp)
            ! while loop to exclude log(0)
            do
                call RANDOM_NUMBER(U)
                if (U(1) > 0.0_wp) exit
            end do
            Z(1)=dsqrt(-2*log(U(1)))*cos(2.0_wp*pi*U(2))
            Z(2)=dsqrt(-2*log(U(1)))*cos(2.0_wp*pi*U(2))
        end function

        function box_muller_vec(N,mu,sigma) result(Z)
            implicit none
            integer(wp), intent(in) :: N
            real(wp), intent(in) :: mu, sigma
            real(wp) :: D(2)
            real(wp),allocatable :: Z(:)
            integer(wp) :: i,j

            allocate(Z(N))
            do i=1,N/2 ! two random numbers are returned for each
                j = 2*i-1
                Z(j:j+1) = box_muller()
            end do
            ! if N odd, then last element will be empty
            if (mod(N,2) .ne. 0) then
                D = box_muller()
                Z(N) = D(1)
            end if
            Z = Z*sigma
            Z = Z + mu
            end function
end module

program random_test
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use random_normal_distr, only: box_muller_vec
    implicit none
    integer(wp) :: j,N
    real(wp) :: old_time, new_time
    real(wp) :: mu, sigma, Xmean
    real(wp),allocatable :: X(:)
    read(*,*) N
    ! N=1000
    allocate(X(N))
    mu=0.0_wp
    sigma=0.0_wp
    call cpu_time(old_time)
    X = box_muller_vec(N,mu,sigma)
    call cpu_time(new_time)
    write(*,*) 'Took',new_time-old_time,'s to generate',N,'numbers'
    Xmean = sum(X)/N
    X= X - Xmean
    print "('central moments',/,a10,*(i10))","mean",1,2,3,4
    print "(*(f10.4))",Xmean,sum(x**1)/N,sum(x**2)/N,sum(x**3)/N,sum(x**4)/N
    print*,"theoretical"
    print "(*(f10.4))",mu,0.0_wp,sigma**2,0.0_wp,3*sigma**4

end program
