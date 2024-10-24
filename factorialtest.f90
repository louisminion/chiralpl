function factorial(j) result(fac)
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    integer, intent(in) :: j
    real(kind=wp) :: fac
    if (j <= 0) then
        write(*,*)  '<ERROR> Factorial undefined for arg', j, '<=0'
        write(*,*) ' <STOP> Bad argument to factorial'
        stop
    end if
    fac = gamma(real(j+1, kind=wp))

end function factorial

program factorialtest
    use iso_fortran_env, only: wp => real64, int64
    implicit none
    integer :: f
    real(kind=wp), external :: factorial
    write(*,*) 'FACTORIAL OF?'
    read(*,*) f
    write(*,*) factorial(f)
end program