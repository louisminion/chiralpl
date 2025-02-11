module output
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    implicit none
    private
    public :: print_line, printOutputHeader, writeDipoleArray,checkRAM
    contains

        subroutine print_line(string, INPUT_NAME)
            ! Prints line to stdlib and mirrors it to an output file.
            ! NOT USED!
            implicit none
            character (*) string
            character (*) INPUT_NAME
            character(:), allocatable :: OUTPUT_NAME
            allocate(character(len=(LEN(TRIM(INPUT_NAME))+4)) :: OUTPUT_NAME)
            write(OUTPUT_NAME,'(a,a)') trim(INPUT_NAME), trim('.out')
            print*, string
            open(19, file=OUTPUT_NAME, position='APPEND', action='WRITE')
            write(19,*) string
            close(19)


        end subroutine

    subroutine printOutputHeader()
        implicit none
        character*256 :: dummy
        character*8  :: date_now
        character*10 :: time_now
        character(len=15) :: start='Program Started'
        print*, '*****************************************'
        print*, '                 chiralpl'
        print*, '*****************************************'
        print*, '            Louis Minion 2024'
        call date_and_time(DATE=date_now,TIME=time_now)
        print*, iso_8601()
        dummy = trim(start) // ' ' // date_now
        print*, trim(dummy)
        print*, 'Beginning input read.'
        print*, 'Echoes input below'
    end subroutine

    subroutine writeDipoleArray(INPUT_NAME,lattice_dimx,lattice_dimy,lattice_dimz,mu_xyz)
        implicit none
        character(256) :: file_name
        character(*) :: INPUT_NAME
        integer(wp), intent(in) :: lattice_dimx,lattice_dimy,lattice_dimz
        real(wp), intent(in) :: mu_xyz(:,:)
        integer(wp) :: j

        file_name = trim(INPUT_NAME) // trim('_dipoles.dat')
        file_name = trim(file_name)
        open(unit=33, file=file_name)
        write(33,'(A,I0,A,I0,A,I0)') 'LATTICE DIMENSIONS (X,Y,Z) ',lattice_dimx,',',lattice_dimy,',',lattice_dimz
        do j = 1, size(mu_xyz,dim=1)
            write(33, '(*(F12.8 : ", "))') mu_xyz(j, 1),mu_xyz(j, 2),mu_xyz(j, 3)
        end do
        close(33)

    end subroutine

    subroutine checkRAM(general_counter)
        implicit none
        real(wp) :: estimated_RAM
        integer,intent(in) :: general_counter
        ! print*,'Hamiltonian size()','(',general_counter,general_counter,')'
        estimated_RAM = ((((1.0_wp*general_counter)**2)*64)/(8.0_wp*1.0E9_wp))
        if ( estimated_RAM> 8.0_wp) then
            print*,estimated_RAM,'GB virtual mem requested, higher than the available 8GB'
            write(*,*) '<STOP> Estimated need for RAM greater than that available, reduce size of basis.'
            STOP
        else
            print*,estimated_RAM, 'GB RAM requested'
        end if

    end subroutine

    function iso_8601()
        ! return date using ISO-8601 format at a resolution of seconds
        character(len=8)  :: dt
        character(len=10) :: tm
        character(len=5)  :: zone
        character(len=25) :: iso_8601
        call date_and_time(dt, tm, zone)
           ISO_8601 = dt(1:4)//'-'//dt(5:6)//'-'//dt(7:8) &
                    & //' '//                             &
                    & tm(1:2)//':'//tm(3:4)//':'//tm(5:6) &
                    & //zone(1:3)//':'//zone(4:5)
        end function iso_8601

end module