module input
    use, intrinsic :: iso_fortran_env, only: wp => real64, int64
    use variables
    implicit none
    private
    public :: readInput
    contains
        subroutine readInput()
            character*255 :: fname
            character*100 :: buffer, label
            integer(wp) :: errstat
            logical :: exists
            integer(wp) :: io_stat, file_no, line_no, pos_space
            parameter (file_no=87)
            call get_command_argument(1, fname, status=errstat )
            if (errstat .ne. 0) then
                write(*,*) 'NO INPUT FILE SPECIFIED'
                write(*,*) '<STOP> Input file required.'
                STOP
            end if
            inquire( file=trim(fname), exist=exists)
            if ( .not. exists ) then
                write(*,*) 'INPUT FILE DOES NOT EXIST'
                write(*,*) '<STOP> Non-existent input file.'
                STOP
            end if
            write(*,*) trim(fname)
            io_stat = 0
            line_no = 0
            open( unit = file_no, file = fname, status = 'old', action = 'read' )
            do while (io_stat .eq. 0)
                read( file_no, '(a)', iostat=io_stat ) buffer
                if (io_stat .eq. 0) then
                    line_no = line_no + 1
                    pos_space = scan( buffer, ' ' )
                    label = buffer( 1:pos_space )
                    buffer = buffer( pos_space + 1:)
                    if ( label(1:1) == '#' ) cycle ! Enable comments with line start #
                    write(*,"(*(a))") '->',trim(label), ' ', trim(buffer)
                    select case ( label )
                    case ( 'INPUT_NAME' )
                        read(buffer, *, iostat=io_stat) INPUT_NAME
                    case ( 'LATTICE_DIMX' )
                        read(buffer, *, iostat=io_stat) lattice_dimx
                    case ( 'LATTICE_DIMY' )
                        read(buffer, *, iostat=io_stat) lattice_dimy
                    case ( 'LATTICE_DIMZ' )
                        read(buffer, *, iostat=io_stat) lattice_dimz
                    case ( 'X_SPACING')
                        read(buffer, *, iostat=io_stat) x_spacing
                    case ( 'Y_SPACING')
                        read(buffer, *, iostat=io_stat) y_spacing
                    case ( 'Z_SPACING')
                        read(buffer, *, iostat=io_stat) z_spacing
                    case ( 'MANUAL_COUPLING' )
                        read(buffer, *, iostat=io_stat) manual_coupling
                    case ( 'NEAREST_NEIGHBOURS' )
                        read(buffer, *, iostat=io_stat) n_nearest_neighbour
                    case ( 'TWIST_ANGLE' )
                        read(buffer, *, iostat=io_stat) phi
                    case ( 'ONE_PARTICLE_STATES' )
                        read(buffer, *, iostat=io_stat) bool_one_particle_states
                    case ( 'TWO_PARTICLE_STATES' )
                        read(buffer, *, iostat=io_stat) bool_two_particle_states
                    case ( 'MAX_VIB')
                        read(buffer, *, iostat=io_stat) max_vibs
                    case ( 'HUANG_RHYS_NEUTRAL')
                        read(buffer, *, iostat=io_stat) lambda_neutral
                    case ( 'MONOMER_TRANSITION')
                        read(buffer, *, iostat=io_stat) w00
                    case ( 'TRANSITION_DIPOLE' )
                        read(buffer, *, iostat=io_stat) mu_0
                    case ( 'TEMPERATURE' )
                        read(buffer, *, iostat=io_stat) temp
                    case ( 'JCOULX' )
                        read(buffer, *, iostat=io_stat) JCoulx
                    case ( 'JCOULY' )
                        read(buffer, *, iostat=io_stat) JCouly
                    case ( 'JCOULZ' )
                        read(buffer, *, iostat=io_stat) JCoulz
                    case ( 'LINEWIDTH' )
                        read(buffer, *,iostat=io_stat) lw
                    case ('WRITE_HAMILTONIAN_OUT')
                        read(buffer, *, iostat=io_stat) H_out
                    case ('SAVE_EIGENVALS')
                        read(buffer, *, iostat=io_stat) save_evals
                    case ('SAVE_EIGENVECS')
                        read(buffer, *, iostat=io_stat) save_evecs
                    case ( 'CORRELATION_LENGTH' )
                        read(buffer, *, iostat=io_stat) l0
                    case ( 'NUMBER_CONFIGURATIONS' )
                        read(buffer, *, iostat=io_stat) configs
                    case ( 'DISORDER_WIDTH' )
                        read(buffer, *, iostat=io_stat) sigma
                    case ( 'THREADS' )
                        read(buffer, *, iostat=io_stat) num_threads
                    case default
                        write(*,"(*(a))") 'Unable to assign input variable: ', label
                    end select
        
                end if
            end do
        
            close(unit =  file_no)
        
            ! Normalise units to hw
            w00 = w00/hw
            lw = lw/hw
            ! hw = hw/hw
            JCoulx = JCoulx * eV / hw
            JCouly = JCouly * eV / hw
            JCoulz = JCoulz * eV / hw
        
        
            phi = pi*(phi/180.0_wp)
            mu_0 = mu_0/Debye
            k = w00/c
        end subroutine



end module