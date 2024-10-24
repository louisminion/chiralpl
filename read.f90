program read
    character*100 :: r, a
    integer :: ios

    read(*,*) r

    read(r,*,iostat=ios) a

    write(*,*) a
    write(*,*) ios, 'IOS'
end program