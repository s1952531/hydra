program split_and_zip_cgc_inputs
    implicit none
    ! ===== User parameters =====
    integer, parameter :: n_total = 23810       ! total number of records
    integer, parameter :: n_per_chunk = 2381    ! records per output file
    ! ============================

    integer,parameter:: ng=128, nt=2*ng !ng set in line 105 of flow-setup
    integer,parameter:: ngridp=ng*nt
    integer,parameter:: npm=200*ngridp


    integer :: i, chunk, ios
    integer:: next(0:npm),npt
    double precision:: x(npm),y(npm),z(npm)

    character(len=64) :: filename, cmd

    open(unit=10, file='cgc_inputs.dat', form='unformatted', access='stream', status='old')

    do i = 1, n_total
        read(10, iostat=ios) x, y, z, next, npt
        if (ios /= 0) exit

        chunk = (i - 1) / n_per_chunk + 1
        write(filename, '(A,I4.4,A)') 'cgc_inputs_', chunk, '.dat'

        open(unit=20, file=filename, form='unformatted', access='stream', &
             position='append', status='unknown')
        write(20) x, y, z, next, npt
        close(20)
    end do

    close(10)

    do chunk = 1, (n_total + n_per_chunk - 1) / n_per_chunk
        write(filename, '(A,I4.4,A)') 'cgc_inputs_', chunk, '.dat'
        write(cmd, '(A,A)') 'gzip -f ', trim(filename)
        call system(cmd)
    end do
end program split_and_zip_cgc_inputs
