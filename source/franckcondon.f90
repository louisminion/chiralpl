module franckcondon
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    private
    public :: calcFranckCondonTables, calcVibrationalOverlap,factorial
    
    contains
        subroutine calcFranckCondonTables()
            ! precalculates vibrational overlap integrals into table
            ! starting with just ground-state to frenkel exciton type <m|n> factors
            ! also calculate fc tables for charge-transfer excitons; vibrational overlap between neutral ground state and cation/anion vibrational manifold.
            implicit none
            integer(wp) :: ground_vib, exc_vib
            real(wp) :: vibolap
            if (.not. allocated(fc_ground_to_neutral)) then
                allocate(fc_ground_to_neutral(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
            end if
            if (.not. allocated(fc_ground_to_cation)) then
                allocate(fc_ground_to_cation(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
            end if
            if (.not. allocated(fc_ground_to_anion)) then
                allocate(fc_ground_to_anion(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
            end if
            if (.not. allocated(fc_cation_to_frenkel)) then
                allocate(fc_cation_to_frenkel(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
            end if
            if (.not. allocated(fc_anion_to_frenkel)) then
                allocate(fc_anion_to_frenkel(0:max_vibs, 0:max_vibs)) ! 2d array with each element vibrational overlap between state i and state j. Zero-indexed
            end if
            do ground_vib=0,max_vibs
                do exc_vib=0,max_vibs
                    call calcVibrationalOverlap(0.0_wp,ground_vib,lambda_neutral,exc_vib,vibolap)
                    fc_ground_to_neutral(ground_vib, exc_vib) = vibolap
                end do
            end do
            do ground_vib=0,max_vibs
                do exc_vib=0,max_vibs
                    call calcVibrationalOverlap(0.0_wp,ground_vib,lambda_plus,exc_vib,vibolap)
                    fc_ground_to_cation(ground_vib, exc_vib) = vibolap
                end do
            end do
            do ground_vib=0,max_vibs
                do exc_vib=0,max_vibs
                    call calcVibrationalOverlap(0.0_wp,ground_vib,lambda_minus,exc_vib,vibolap)
                    fc_ground_to_anion(ground_vib, exc_vib) = vibolap
                end do
            end do
            do ground_vib=0,max_vibs
                do exc_vib=0,max_vibs
                    call calcVibrationalOverlap(lambda_minus,ground_vib,lambda_neutral,exc_vib,vibolap)
                    fc_anion_to_frenkel(ground_vib, exc_vib) = vibolap
                end do
            end do
            do ground_vib=0,max_vibs
                do exc_vib=0,max_vibs
                    call calcVibrationalOverlap(lambda_plus,ground_vib,lambda_neutral,exc_vib,vibolap)
                    fc_cation_to_frenkel(ground_vib, exc_vib) = vibolap
                end do
            end do
            !FCWRITE


        end subroutine


        subroutine calcVibrationalOverlap(lambda_1,n,lambda_2,m,vibolap)
            ! calc <m|n> via the recursion formula
            implicit none
            integer(wp) , intent (in) :: n, m
            integer(wp) :: j
            real(kind=wp), intent(in):: lambda_1, lambda_2
            real(kind=wp) ::  vibolap
            real(kind=wp) :: lambda_hr

            lambda_hr = lambda_2 - lambda_1
            vibolap = 0.0_wp
            do j=0,min(n,m)

                vibolap = vibolap + (((-1.0_wp)**(m-j))/(factorial(n-j)*factorial(j)*factorial(m-j)))*   &
                (lambda_hr**(n+m-2*j))

            end do
            vibolap = vibolap*dsqrt(1.0_wp*factorial(n)*    &
            factorial(m))*dexp(-1.0_wp*lambda_hr**2/2.0_wp)
        end subroutine


        function factorial(j) result(fac)
            ! calc factorial(j) via gamma(j+1) = j!
            use iso_fortran_env, only: wp => real64, int64
            implicit none
            integer(wp), intent(in) :: j
            real(kind=wp) :: fac
            if (j < 0) then
                write(*,*)  '<ERROR> Factorial undefined for arg', j, '<=0'
                write(*,*) ' <STOP> Bad argument to factorial'
                stop
            end if
            fac = gamma(real(j+1, kind=wp))

        end function factorial

end module