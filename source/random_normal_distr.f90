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