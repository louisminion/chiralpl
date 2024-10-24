

program vibolaptest
    use variables
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    integer :: n, m
    real(kind=wp) :: lambda_hr, vibolap
    real*8, external :: volap
    write(*,*) 'GET n'
    read(*,*) n
    write(*,*) 'GET m'
    read(*,*) m
    lambda_hr = 1.0_wp    
    call calcVibrationalOverlap(0.0_wp,n,lambda_hr,m,vibolap)
    write(*,*) 'vibolap', vibolap
    write(*,*) 'volap', volap(0.0_wp,n,lambda_hr,m)


end program


function factorial(j) result(fac)
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    integer, intent(in) :: j
    real(kind=wp) :: fac
    if (j < 0) then
        write(*,*)  '<ERROR> Factorial undefined for arg', j, '<=0'
        write(*,*) ' <STOP> Bad argument to factorial'
        stop
    end if
    fac = gamma(real(j+1, kind=wp))

end function factorial


function volap( lambda1, vib1, lambda2, vib2 )
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    real(kind=wp) :: volap
    integer, intent (in) :: vib1, vib2
    real(kind=wp), intent(in) :: lambda1, lambda2
    integer k
    real(kind=wp), external :: factorial
    real(kind=wp) lambda
    
    !calculate the displacement between the two potential wells
    lambda = lambda2 - lambda1

    volap = 0.d0
    !calculate the vibrational overlap
    !first calculate the summation
    do k = 0, min( vib1, vib2 )
        volap = volap+(-1.0_wp)**(vib2-k)/        &
                (factorial(vib1-k)*factorial(k)*&
                 factorial(vib2-k))*            &
                 lambda**(vib1+vib2-2*k)
    end do
    volap = volap*dsqrt(1.0_wp*factorial(vib1)*   &
                             factorial(vib2))*  &
                             dexp(-1.0_wp*        &
                             lambda**2/2.0_wp)

end function


subroutine calcVibrationalOverlap(lambda_1,n,lambda_2,m,vibolap)
    use variables
    implicit none
    real(kind=wp), external :: factorial
    integer , intent (in) :: n, m
    integer :: j
    real(kind=wp), intent(in):: lambda_1, lambda_2
    real(kind=wp) ::  vibolap
    real(kind=wp) :: lambda_hr

    lambda_hr = lambda_2 - lambda_1
    vibolap = 0.0_wp
    do j=0,min(n,m)

        vibolap = vibolap + ((-1.0_wp)**(m-j)/(factorial(n-j)*factorial(j)*factorial(m-j)))*   &
        (lambda_hr**(n+m-2*j))

    end do
    vibolap = vibolap*dsqrt(1.0_wp*factorial(n)*    &
    factorial(m))*dexp(-1.0_wp*lambda_hr**2/2.0_wp)

end subroutine