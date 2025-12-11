program cgcDev
    use omp_lib
    implicit none

    !constants for the default test case
    !twopi, dlfi, pi, dq, ngf, ntf, hpi, clonf, slonf, f12, f14
    double precision,parameter:: one=1.d0, two=2.d0, four=4.d0
    double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
    double precision,parameter:: hpi=pi/two, twopi=two*pi
    double precision,parameter:: f12=one/two, f14=one/four
    double precision,parameter:: small=1.d-12
    double precision,parameter:: fpole=two*6.2831853071 !set in flow-setup as fpole=twopi usin shorter pi and 10 dps (see line 191)
    integer,parameter:: ncont=80 !default as 80 in line 430 of flow-setup
    double precision,parameter:: dq=two*fpole/dble(ncont)
    integer,parameter:: ng=128, nt=2*ng !ng set in line 105 of flow-setup
    integer,parameter:: mgf=4, ngf=ng*mgf, ntf=nt*mgf
    double precision:: clonf(ntf),slonf(ntf)
    integer:: i, j, iterCount
    double precision:: rlonf
    double precision:: dlf,dlfi
    integer,parameter:: ngridp=ng*nt
    integer,parameter:: npm=200*ngridp
    double precision,parameter:: dl=twopi/dble(nt)
    !end of vars/constants outside of contours

    !defined in contours
    double precision:: x(npm),y(npm),z(npm)
    integer:: next(0:npm),npt
    integer(kind=8):: callCount !as used in calc of pos for large file
    double precision:: fcor(ng)

    integer:: numInputs
    integer, parameter:: numIters=1 !238
    integer:: totReads

    double precision:: qc(ng,nt)
    double precision, allocatable:: qcDiffs(:) !of shape (totReads)
                                               !according to gprof there are 23810 calls for the default test case

    double precision, allocatable:: x_arr(:, :) !of shape (npm, numInputs)
    double precision, allocatable:: y_arr(:, :) !of shape (npm, numInputs)
    double precision, allocatable:: z_arr(:, :) !of shape (npm, numInputs)
    integer, allocatable:: next_arr(:, :) !of shape (0:npm, numInputs)
    integer, allocatable:: npt_arr(:) !of shape (numInputs)

    double precision, allocatable:: qc_arr(:,:,:) !of shape (ng,nt, numInputs)

    integer:: qcSize=size(qc) * storage_size(qc)/8 !div by 8 to get bytes
    integer(kind=8):: input_block_size= storage_size(x)/8 * size(x) + &
                                storage_size(y)/8 * size(y) + &
                                storage_size(z)/8 * size(z) + &
                                storage_size(next)/8 * size(next) + &
                                storage_size(npt)/8


    !Timers
    double precision:: con2gridToTTime=0.0d0, preAvgTotTime=0.0d0, avgTotTime=0.0d0
    double precision:: l1TotTime=0.0d0, l2TotTime=0.0d0, l3TotTime=0.0d0, l4TotTime=0.0d0, l5TotTime=0.0d0
    double precision:: l6TotTime=0.0d0
    double precision:: combineTotTime=0.0d0

    !hard coded repeat ks
    integer, allocatable :: callsRepeatKs(:)
    integer, allocatable :: nonRepeatKs(:)

    !for load balancing repeat Ks
    type :: group_t
        integer, allocatable :: values(:)
    end type group_t

    type(group_t), allocatable :: groups(:)
    integer :: ngroups = 0
    character(len=:), allocatable :: filename

    call init
    
    !main loop
    do iterCount=1,numIters
        do callCount = 1, numInputs
            x = x_arr(:, callCount)
            y = y_arr(:, callCount)
            z = z_arr(:, callCount)
            next = next_arr(:, callCount)
            npt = npt_arr(callCount)
            call readRepeatKs
            call getNonRepeatKs
            call getBalancedRepeatKs

            !call checkAllKIncluded

            call con2grid_balancedRepeatKs(qc)
            !call con2grid_masterRepeatKs(qc)
            !call con2grid_serial(qc)
            call compare_qcs
        end do
    end do

    !print the max of the qcDiffs
    print *, 'Max difference in qc: ', maxval(qcDiffs)

    call finish

    contains

    subroutine init
	print *, 'Initializing...'
        call openFiles
        call initVars
        
        do callCount = 1, numInputs
            call readInput
            !call readInputReversed
            call readOutputs
            !call readOutputsReversed
        end do
        call closeFiles
	print *, 'Initialized'
    end subroutine


    subroutine initVars  
        dlf =twopi/dble(ntf)

        do i=1,ntf
            rlonf=dlf*dble(i-1)-pi
            clonf(i)=cos(rlonf)
            slonf(i)=sin(rlonf)
        enddo

        do j=1,ng
            fcor(j)=fpole*sin((dble(j)-f12)*dl-hpi)
        enddo

        dlfi=dble(ntf)/(twopi+small)

        call getNumWrites

        totReads=numInputs*numIters

        call allocateVars

        return
    end subroutine

    subroutine allocateVars
        allocate(x_arr(npm, numInputs))
        allocate(y_arr(npm, numInputs))
        allocate(z_arr(npm, numInputs))
        allocate(next_arr(0:npm, numInputs))
        allocate(npt_arr(numInputs))
        allocate(qc_arr(ng,nt, numInputs))
        allocate(qcDiffs(totReads))
        return
    end subroutine

    subroutine getNumWrites
        !determine how many writes are in cgc_* files
        !inputs and outputs written same number of times, easier to calc size of cgc_outputs.dat since only qc written
        integer :: filesize
        
        inquire(unit=102, size=filesize)
        !print *, 'File size of cgc_outputs.dat: ', filesize
        !print *, 'Size of one qc array: ', qcSize
        numInputs=filesize/qcSize
        !print *, 'Number of inputs/outputs in files: ', numInputs

        return
    end subroutine


    subroutine openFiles
        open(100, file="cgc_inputs.dat", status='old', action='read', access='stream', form='unformatted')
        open(102, file="cgc_outputs.dat", status='old', action='read', access='stream', form='unformatted')
        return
    end subroutine

    subroutine closeFiles
        close(100)
        close(102)
        return
    end subroutine

    subroutine readRepeatKs

        ! integer :: unit, val, count, ios
        ! character(len=256) :: filename

        ! print *, 'Loading repeat ks for call ', callCount

        ! ! Build filename from global callCount
        ! write(filename, '(A,I0,A)') 'RepeatKFiles/repeatKs_call_', callCount, '.txt'

        ! ! Deallocate previous array if allocated
        ! if (allocated(callsRepeatKs)) then
        !     deallocate(callsRepeatKs)
        ! end if
        ! ! Allocate array with global npt
        ! allocate(callsRepeatKs(npt))
        ! count = 0

        ! ! Open file
        ! open(newunit=unit, file=filename, status='old', action='read')

        ! ! Read integers directly into data
        ! do
        !     read(unit, *, iostat=ios) val
        !     if (ios /= 0) exit
        !     count = count + 1
        !     callsRepeatKs(count) = val
        ! end do

        ! print *, 'Loaded'

        ! close(unit)

        integer :: unit, val, count, ios
        character(len=256) :: filename
        integer, allocatable :: temp(:)  ! temporary array, allocatable for move_alloc

        !print *, 'Loading repeat ks for call ', callCount

        ! Build filename from global callCount
        write(filename, '(A,I0,A)') 'RepeatKFiles/repeatKs_call_', callCount, '.txt'

        ! Open file
        open(newunit=unit, file=filename, status='old', action='read')

        ! Allocate temporary array with size npt
        allocate(temp(npt))

        ! Read integers directly into temp
        count = 0
        do
            read(unit, *, iostat=ios) val
            if (ios /= 0) exit
            count = count + 1
            if (count > npt) then
                print *, 'Warning: more repeat Ks in file than npt, truncating at npt'
                exit
            end if
            temp(count) = val
        end do

        !print *, 'read ', count, ' repeat ks from file.'

        close(unit)

        ! Allocate callsRepeatKs to exact size using move_alloc

        if (allocated(callsRepeatKs)) then
            deallocate(callsRepeatKs)
        end if
        allocate(callsRepeatKs(count))

        callsRepeatKs(:) = temp(1:count)

        deallocate(temp)

        !print *, 'Loaded ', count, ' repeat ks.'

    end subroutine

    subroutine getBalancedRepeatKs
         character(len=:), allocatable :: line
        integer :: unit, ios, i, nvals
        integer, allocatable :: tmp(:)

        ! Build filename from global callCount
        write(filename, '(A,I0,A)') 'balancedRepeatKs_call_', callCount, '.txt'

        ! ------------------------------------------------
        ! First pass: count how many non-empty lines
        ! ------------------------------------------------
        ngroups = 0
        open(newunit=unit, file=filename, status='old', action='read')

        allocate(character(200000) :: line)

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) > 0) ngroups = ngroups + 1
        end do
        close(unit)

        if (ngroups == 0) then
            print *, "WARNING: no groups found in ", filename
            return
        end if

        ! Allocate global groups array
        allocate(groups(ngroups))

        ! ------------------------------------------------
        ! Second pass: parse groups
        ! ------------------------------------------------
        open(newunit=unit, file=filename, status='old', action='read')

        i = 0
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle

            i = i + 1

            ! Count integers: number of spaces + 1
            nvals = count([(line(j:j) == ' ', j=1,len_trim(line))]) + 1

            allocate(tmp(nvals))

            ! Parse values
            read(line, *, iostat=ios) tmp
            if (ios /= 0) then
                print *, "ERROR parsing group ", i, ": ", trim(line)
                stop
            end if

            ! Save to global array
            allocate(groups(i)%values(nvals))
            groups(i)%values = tmp
        end do

        close(unit)

    end subroutine

    subroutine checkAllKIncluded
        !check if all ks from 1 to npt are included in callsRepeatKs and nonRepeatKs
        integer :: k, kk
        integer :: k_count(npt)

        print *, 'Checking all ks included for call ', callCount

        do k=1,npt
            k_count(k)=0
        enddo

        do kk=1,size(callsRepeatKs)
            k=callsRepeatKs(kk)
            !print *, 'k in callsRepeatKs: ', k
            k_count(k)=k_count(k)+1
        enddo

        do kk=1,size(nonRepeatKs)
            k=nonRepeatKs(kk)
            k_count(k)=k_count(k)+1
        enddo

        do k=1,npt
            if (k_count(k) /= 1) then
                print *, 'Error: k=', k, ' count=', k_count(k)
            endif
        enddo

    end subroutine checkAllKIncluded

    subroutine readInput
        ! Reads in the contour data from a file "cgc_inputs.dat"

        integer(kind=8) :: current_pos
        inquire(unit=100, pos=current_pos)
        !print *, 'callCount: ', callCount
        !print *, 'Reading input at position: ', current_pos

        read(100) x_arr(:, callCount)
        read(100) y_arr(:, callCount)
        read(100) z_arr(:, callCount)
        read(100) next_arr(:, callCount)
        read(100) npt_arr(callCount)
        return
    end subroutine

    subroutine readInputReversed
        !Reads in the contour data from a file "cgc_inputs.dat" in reverse order
        !to test work balance

        integer(kind=8):: filesize !kind=8 to hold large file sizes
        integer(kind=8):: pos
        
        inquire(unit=100, size=filesize)

        !print *, 'File size of cgc_inputs.dat: ', filesize
        !print *, 'Input block size: ', input_block_size
        !print *, 'Call count: ', callCount

        pos = 1_8 + filesize - int(input_block_size*callCount, kind=8)

        !print *, 'Reading input at position: ', pos

        read(100, pos=pos) x_arr(:, callCount)
        read(100) y_arr(:, callCount)
        read(100) z_arr(:, callCount)
        read(100) next_arr(:, callCount)
        read(100) npt_arr(callCount)

    end subroutine

    subroutine getNonRepeatKs
        ! !remove callsRepeatKs from array of all ks 1 to npt
        ! integer :: totalKs
        ! integer :: i, idx

        ! !print *, 'Calculating non-repeat ks for call ', callCount
        ! !print *, 'npt=', npt, ' size(callsRepeatKs)=', size(callsRepeatKs)

        ! !print*, 'Allocation of nonRepeatKs...'
        ! if (allocated(nonRepeatKs)) then
        !     !print *, 'Deallocating previous nonRepeatKs...'
        !     deallocate(nonRepeatKs)
        ! end if
        ! allocate(nonRepeatKs(npt - size(callsRepeatKs)))
        ! !print *, 'Allocated'

        ! !initialize nonRepeatKs to all ks
        ! !print *, 'Initializing non-repeat ks...'    
        ! do i=1,npt
        !     nonRepeatKs(i) = i
        ! enddo

        ! !print *, 'Flagging repeat ks...'
        ! !flag repeat ks for removal
        ! do i=1,size(callsRepeatKs)
        !     idx = callsRepeatKs(i)
        !     !print *, 'marking k ', idx, ' as repeat'
        !     nonRepeatKs(idx) = -1 !mark as removed
        ! enddo

        ! !print *, 'Compacting non-repeat ks...'
        ! !compact array to only non-repeat ks
        ! totalKs = 0
        ! do i=1,npt
        !     if (nonRepeatKs(i) /= -1) then
        !         totalKs = totalKs + 1
        !         !print *, 'keeping k ', nonRepeatKs(i), 'in position ', totalKs
        !         nonRepeatKs(totalKs) = nonRepeatKs(i)
        !     endif
        ! enddo

        ! Remove the ks listed in callsRepeatKs from 1..npt.
    ! callsRepeatKs is guaranteed to contain unique values.

    integer :: i, j
    logical, allocatable :: isRepeat(:)
    integer :: nNonRepeat

    ! Mask to mark repeats
    allocate(isRepeat(npt))
    isRepeat = .false.

    ! Mark repeat ks
    do i = 1, size(callsRepeatKs)
        isRepeat(callsRepeatKs(i)) = .true.
    end do

    ! Number of ks that are NOT in callsRepeatKs
    nNonRepeat = npt - size(callsRepeatKs)

    ! Allocate exact-size result array
    if (allocated(nonRepeatKs)) deallocate(nonRepeatKs)
    allocate(nonRepeatKs(nNonRepeat))

    ! Fill the array
    j = 0
    do i = 1, npt
        if (.not. isRepeat(i)) then
            j = j + 1
            nonRepeatKs(j) = i
        end if
    end do

    deallocate(isRepeat)

    end subroutine getNonRepeatKs

    subroutine con2grid_serial(qc)

         implicit double precision(a-h,o-z)
        implicit integer(i-n)

        !Passed arrays:
        double precision:: qc(ng,nt)
        !Local arrays:
        double precision:: qa(0:ngf+1,ntf)
        double precision:: qa_jp1(0:ngf+1,ntf)
        double precision:: qaend(ngf/2)
        integer:: ilm1(npt),ntc(npt)
        double precision:: cx(npt),cy(npt),cz(npt)
        double precision:: sq(npt)

        !timers (overall and one for each loop)
        double precision:: startTime, endTime, totalTime
        double precision:: l1Start, l1End, l1Time
        double precision:: l2Start, l2End, l2Time
        double precision:: l3Start, l3End, l3Time
        double precision:: l4Start, l4End, l4Time
        double precision:: l5Start, l5End, l5Time
        double precision:: l6Start, l6End, l6Time

        double precision:: combineStart, combineEnd, combineTime

        double precision:: preAvgStart, preAvgEnd, preAvgTime
        double precision:: avgStart, avgEnd, avgTime

        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

        startTime = omp_get_wtime()
        preAvgStart = startTime

        !Initialise crossing information:

	    
	    !print *, 'Loop 1...'
                !LOOP 1
                l1Start = omp_get_wtime()
                do k=1,npt
                    ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
                enddo
                l1End = omp_get_wtime()
                l1Time = l1End - l1Start

            !print *, 'Loop 2...'
                !LOOP 2
                l2Start = omp_get_wtime()
                do k=1,npt
                    ka=next(k)
                    cx(k)=z(k)*y(ka)-y(k)*z(ka)
                    cy(k)=x(k)*z(ka)-z(k)*x(ka)
                    cz(k)=x(k)*y(ka)-y(k)*x(ka)
                    ntc(k)=ilm1(ka)-ilm1(k)
                enddo    
                l2End = omp_get_wtime()
                l2Time = l2End - l2Start
                
		!print *, 'Loop 3...'
                !LOOP 3
                l3Start = omp_get_wtime()
                do k=1,npt
                    sig=sign(one,cz(k))
                    sq(k)=dq*sig
                    ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
                    if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
                        if (abs(cz(k)) .gt. zero) then
                            cx(k)=cx(k)/cz(k)
                            cy(k)=cy(k)/cz(k)
                        endif
                enddo
                l3End = omp_get_wtime()
                l3Time = l3End - l3Start

        !----------------------------------------------------------------------
        !Initialise PV jump array:
	!print *, 'Loop 4...'
        !LOOP 4
        l4Start = omp_get_wtime()        
        do i=1,ntf
            do j=0,ngf+1
                qa(j,i)=zero
                qa_jp1(j,i)=zero
            enddo
        enddo
        l4End = omp_get_wtime()
        l4Time = l4End - l4Start

        !Determine crossing indices:
	!LOOP 5
        l5Start = omp_get_wtime()
        !print *, 'Loop 5...'

        !do k=1,npt
        do k=1,npt
            if (ntc(k) .ne. 0) then
                jump=sign(1,ntc(k))
                ioff=ntf+ilm1(k)+(1+jump)/2
                ncr=0
                do while (ncr .ne. ntc(k))
                    i=1+mod(ioff+ncr,ntf)
                    ncr=ncr+jump

                    rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                    j=int(rlatc)+1
                    p=rlatc-dble(j-1)
                    qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                    qa(j+1,i)=qa(j+1,i)+    p*sq(k)
                enddo
            endif
        enddo

        l5End = omp_get_wtime()
        l5Time = l5End - l5Start

        !Get PV values, at half latitudes, by sweeping through latitudes:
        !LOOP 6
        l6Start = omp_get_wtime()
        do i=1,ntf
            do j=2,ngf
                qa(j,i)=qa(j,i)+qa(j-1,i)
            enddo
        enddo
        l6End = omp_get_wtime()
        l6Time = l6End - l6Start
        !Here, qa(j,i) stands for the PV at latitude j-1/2,
        !from j = 1, ..., ngf.

        preAvgEnd = omp_get_wtime()
        preAvgTime = preAvgEnd - preAvgStart
    
        !----------------------------------------------------------------------
        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

        avgStart = omp_get_wtime()

        !Average PV values on the fine grid to get corresponding 
        !values on the inversion grid (ng,nt):
        ngh=ngf
        nth=ntf
        
        do while (ngh .gt. ng) !need to precalc number of iters if want to parallelise with OMP
            !Pre-store PV adjacent to poles at complementary longitudes (+pi):
            nthh=nth/2
            nghp1=ngh+1
            do i=1,nthh
                ic=i+nthh
                qa(0,i)=qa(1,ic)
                qa(0,ic)=qa(1,i)
                qa(nghp1,i)=qa(ngh,ic)
                qa(nghp1,ic)=qa(ngh,i)
            enddo

            !Work from SP to NP to define PV at full latitudes from averages
            !at adjacent half latitudes:
            do i=1,nth
                do j=0,ngh
                qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
                enddo
            enddo

            !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

            !Next 1-2-1 average these values to define PV at half latitudes
            !on a grid twice as coarse:
            nghh=ngh/2
            do i=1,nth
                do j=1,nghh
                je=2*j
                qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
                enddo
            enddo

            !Now perform analogous longitudinal 1-2-1 average:
            do j=1,nghh
                qaend(j)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=1,nth-1
                ip1=i+1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
                enddo
            enddo

            do j=1,nghh
                qa(j,nth)=qaend(j)
            enddo
            !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

            !Average these on the twice coarser grid:
            nthh=nth/2

            do j=1,nghh
                qa(j,1)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=2,nthh
                io=2*i-1
                ie=io-1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,ie)+qa(j,io))
                enddo
            enddo

            ngh=nghh
            nth=nthh

        enddo

        !Finalise and take away f to define PV anomaly:
        do i=1,nt
            do j=1,ng
                qc(j,i)=qa(j,i)-fcor(j)
            enddo
        enddo

        avgEnd = omp_get_wtime()
        avgTime = avgEnd - avgStart

        endTime = omp_get_wtime()
        totalTime = endTime - startTime
        !print *, 'call', callCount, ' of con2grid took ', totalTime, ' seconds.'

        call accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, combineTime)
        
        return
        
    end subroutine

    subroutine con2grid_balancedRepeatKs(qc)
        ! Calculates the PV anomaly field (stored in qc) from the PV 
        ! contours (x,y,z).  Takes away Coriolis frequency (fcor).
	
	    use omp_lib

        implicit double precision(a-h,o-z)
        implicit integer(i-n)

        !Passed arrays:
        double precision:: qc(ng,nt)
        !Local arrays:
        double precision:: qa(0:ngf+1,ntf)
        double precision:: qa_jp1(0:ngf+1,ntf)
        double precision:: qaend(ngf/2)
        integer:: ilm1(npt),ntc(npt)
        double precision:: cx(npt),cy(npt),cz(npt)
        double precision:: sq(npt)

        !timers (overall and one for each loop)
        double precision:: startTime, endTime, totalTime
        double precision:: l1Start, l1End, l1Time
        double precision:: l2Start, l2End, l2Time
        double precision:: l3Start, l3End, l3Time
        double precision:: l4Start, l4End, l4Time
        double precision:: l5Start, l5End, l5Time
        double precision:: l6Start, l6End, l6Time

        double precision:: combineStart, combineEnd, combineTime

        double precision:: preAvgStart, preAvgEnd, preAvgTime
        double precision:: avgStart, avgEnd, avgTime

        integer:: groupCount, ki

        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

	    !print *, "qa size (bytes): ", size(qa) * storage_size(qa)/8

        startTime = omp_get_wtime()
        preAvgStart = startTime

        !Initialise crossing information:

        !!$OMP PARALLEL DEFAULT(NONE) SHARED(npt,dlfi,x,y,z,next,ilm1,cx,cy,cz,ntc,sq, zero) PRIVATE(k,ka,sig)
	    
	    !print *, 'Loop 1...'
            !!$OMP DO SCHEDULE(STATIC)
                !LOOP 1
                l1Start = omp_get_wtime()
                !!$OMP PARALLEL DO SCHEDULE(GUIDED)
                do k=1,npt
                    ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
                enddo
                !!$OMP END PARALLEL DO
                l1End = omp_get_wtime()
                l1Time = l1End - l1Start
            !!$OMP END DO

            !print *, 'Loop 2...'
	    !!$OMP DO SCHEDULE(STATIC)
                !LOOP 2
                l2Start = omp_get_wtime()
                do k=1,npt
                    ka=next(k)
                    cx(k)=z(k)*y(ka)-y(k)*z(ka)
                    cy(k)=x(k)*z(ka)-z(k)*x(ka)
                    cz(k)=x(k)*y(ka)-y(k)*x(ka)
                    ntc(k)=ilm1(ka)-ilm1(k)
                enddo    
                l2End = omp_get_wtime()
                l2Time = l2End - l2Start
                
		!print *, 'Loop 3...'
                !LOOP 3
                l3Start = omp_get_wtime()
                do k=1,npt
                    sig=sign(one,cz(k))
                    sq(k)=dq*sig
                    ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
                    if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
                        if (abs(cz(k)) .gt. zero) then
                            cx(k)=cx(k)/cz(k)
                            cy(k)=cy(k)/cz(k)
                        endif
                enddo
                l3End = omp_get_wtime()
                l3Time = l3End - l3Start
            !!$OMP END DO
        !!$OMP END PARALLEL

        !----------------------------------------------------------------------
        !Initialise PV jump array:
	!print *, 'Loop 4...'
        !LOOP 4
        l4Start = omp_get_wtime()        
        do i=1,ntf
            do j=0,ngf+1
                qa(j,i)=zero
                qa_jp1(j,i)=zero
            enddo
        enddo
        l4End = omp_get_wtime()
        l4Time = l4End - l4Start

        !Determine crossing indices:
	!LOOP 5
        l5Start = omp_get_wtime()
        !print *, 'Loop 5...'
    
    !$OMP PARALLEL PRIVATE(k,j,i,ioff,ncr,rlatc,p,jump, groupCount, ki)
        !split k=1,npt into repeat and non-repeat ks. 
        !hard coded load of pre balanced repeat ks
        !$omp do
        do groupCount = 1, ngroups
            do ki = 1, size(groups(groupCount)%values)
                k = groups(groupCount)%values(ki)
                if (ntc(k) .ne. 0) then
                    jump=sign(1,ntc(k))
                    ioff=ntf+ilm1(k)+(1+jump)/2
                    ncr=0
                    do while (ncr .ne. ntc(k))
                        i=1+mod(ioff+ncr,ntf)
                        ncr=ncr+jump

                        rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                        j=int(rlatc)+1
                        p=rlatc-dble(j-1)
                        qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                        qa_jp1(j+1,i)=qa_jp1(j+1,i)+    p*sq(k)
                    enddo
                endif
            enddo
        enddo
        !$omp end do
        
        !$OMP DO
        do kk=1,size(nonRepeatKs)
            k=nonRepeatKs(kk)
        !if (mod(k,1) .eq. 0) then
                ! !$OMP CRITICAL
                ! !print *, 'Thread ', omp_get_thread_num(), ' at k=', k
                ! !$OMP END CRITICAL
            !endif
            if (ntc(k) .ne. 0) then
                jump=sign(1,ntc(k))
                ioff=ntf+ilm1(k)+(1+jump)/2
                ncr=0
                do while (ncr .ne. ntc(k))
                    i=1+mod(ioff+ncr,ntf)
                    ncr=ncr+jump
                    
                    !check if i in thread's range (start to end)
                    !if (i < start .or. i > end) cycle

                    rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                    j=int(rlatc)+1
                    p=rlatc-dble(j-1)
                    !!$OMP CRITICAL
                    !print *, 'Thread ', omp_get_thread_num(), 'at k=', k, 'i, j: ', i, ',', j
                    !!$OMP ATOMIC
                    qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                    !!$OMP ATOMIC
                    qa_jp1(j+1,i)=qa_jp1(j+1,i)+    p*sq(k)
                    !!$OMP END CRITICAL

                    !print *, 'i,j', i, ',', j
                enddo
            endif
        enddo
        !$OMP END DO

        !combine qa and qa_jp1 into qa
        
        ! !$omp do collapse(2)
        ! do j = 0, ngf+1
        !     do i = 1, ntf
        !         qa(j,i) = qa(j,i) + qa_jp1(j,i)
        !     end do
        ! end do
        ! !$omp end do

        !!$OMP END PARALLEL DO
    !$OMP END PARALLEL

        ! combineStart = omp_get_wtime()
        ! !$omp parallel do collapse(2)
        ! do j = 0, ngf+1
        !     do i = 1, ntf
        !         qa(j,i) = qa(j,i) + qa_jp1(j,i)
        !     end do
        ! end do
        ! !$omp end parallel do
        ! combineEnd = omp_get_wtime()
        ! combineTime = combineEnd - combineStart

        ! combineStart = omp_get_wtime()
        !combine qa and qa_jp1 into qa serially
        qa = qa + qa_jp1
        combineEnd = omp_get_wtime()
        ! combineTime = combineEnd - combineStart

        l5End = omp_get_wtime()
        l5Time = l5End - l5Start

        !Get PV values, at half latitudes, by sweeping through latitudes:
        !LOOP 6
        l6Start = omp_get_wtime()
        do i=1,ntf
            do j=2,ngf
                qa(j,i)=qa(j,i)+qa(j-1,i)
            enddo
        enddo
        l6End = omp_get_wtime()
        l6Time = l6End - l6Start
        !Here, qa(j,i) stands for the PV at latitude j-1/2,
        !from j = 1, ..., ngf.

        preAvgEnd = omp_get_wtime()
        preAvgTime = preAvgEnd - preAvgStart
    
        !----------------------------------------------------------------------
        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

        avgStart = omp_get_wtime()

        !Average PV values on the fine grid to get corresponding 
        !values on the inversion grid (ng,nt):
        ngh=ngf
        nth=ntf
        
        do while (ngh .gt. ng) !need to precalc number of iters if want to parallelise with OMP
            !Pre-store PV adjacent to poles at complementary longitudes (+pi):
            nthh=nth/2
            nghp1=ngh+1
            do i=1,nthh
                ic=i+nthh
                qa(0,i)=qa(1,ic)
                qa(0,ic)=qa(1,i)
                qa(nghp1,i)=qa(ngh,ic)
                qa(nghp1,ic)=qa(ngh,i)
            enddo

            !Work from SP to NP to define PV at full latitudes from averages
            !at adjacent half latitudes:
            do i=1,nth
                do j=0,ngh
                qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
                enddo
            enddo

            !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

            !Next 1-2-1 average these values to define PV at half latitudes
            !on a grid twice as coarse:
            nghh=ngh/2
            do i=1,nth
                do j=1,nghh
                je=2*j
                qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
                enddo
            enddo

            !Now perform analogous longitudinal 1-2-1 average:
            do j=1,nghh
                qaend(j)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=1,nth-1
                ip1=i+1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
                enddo
            enddo

            do j=1,nghh
                qa(j,nth)=qaend(j)
            enddo
            !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

            !Average these on the twice coarser grid:
            nthh=nth/2

            do j=1,nghh
                qa(j,1)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=2,nthh
                io=2*i-1
                ie=io-1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,ie)+qa(j,io))
                enddo
            enddo

            ngh=nghh
            nth=nthh

        enddo

        !Finalise and take away f to define PV anomaly:
        do i=1,nt
            do j=1,ng
                qc(j,i)=qa(j,i)-fcor(j)
            enddo
        enddo

        avgEnd = omp_get_wtime()
        avgTime = avgEnd - avgStart

        endTime = omp_get_wtime()
        totalTime = endTime - startTime
        !print *, 'call', callCount, ' of con2grid took ', totalTime, ' seconds.'

        call accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, combineTime)
        
        return
        
    end subroutine 

    subroutine con2grid_masterRepeatKs(qc)
        ! Calculates the PV anomaly field (stored in qc) from the PV 
        ! contours (x,y,z).  Takes away Coriolis frequency (fcor).
	
	    use omp_lib

        implicit double precision(a-h,o-z)
        implicit integer(i-n)

        !Passed arrays:
        double precision:: qc(ng,nt)
        !Local arrays:
        double precision:: qa(0:ngf+1,ntf)
        double precision:: qa_jp1(0:ngf+1,ntf)
        double precision:: qaend(ngf/2)
        integer:: ilm1(npt),ntc(npt)
        double precision:: cx(npt),cy(npt),cz(npt)
        double precision:: sq(npt)

        !timers (overall and one for each loop)
        double precision:: startTime, endTime, totalTime
        double precision:: l1Start, l1End, l1Time
        double precision:: l2Start, l2End, l2Time
        double precision:: l3Start, l3End, l3Time
        double precision:: l4Start, l4End, l4Time
        double precision:: l5Start, l5End, l5Time
        double precision:: l6Start, l6End, l6Time

        double precision:: combineStart, combineEnd, combineTime

        double precision:: preAvgStart, preAvgEnd, preAvgTime
        double precision:: avgStart, avgEnd, avgTime

        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

	    !print *, "qa size (bytes): ", size(qa) * storage_size(qa)/8

        startTime = omp_get_wtime()
        preAvgStart = startTime

        !Initialise crossing information:

        !!$OMP PARALLEL DEFAULT(NONE) SHARED(npt,dlfi,x,y,z,next,ilm1,cx,cy,cz,ntc,sq, zero) PRIVATE(k,ka,sig)
	    
	    !print *, 'Loop 1...'
            !!$OMP DO SCHEDULE(STATIC)
                !LOOP 1
                l1Start = omp_get_wtime()
                !!$OMP PARALLEL DO SCHEDULE(GUIDED)
                do k=1,npt
                    ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
                enddo
                !!$OMP END PARALLEL DO
                l1End = omp_get_wtime()
                l1Time = l1End - l1Start
            !!$OMP END DO

            !print *, 'Loop 2...'
	    !!$OMP DO SCHEDULE(STATIC)
                !LOOP 2
                l2Start = omp_get_wtime()
                do k=1,npt
                    ka=next(k)
                    cx(k)=z(k)*y(ka)-y(k)*z(ka)
                    cy(k)=x(k)*z(ka)-z(k)*x(ka)
                    cz(k)=x(k)*y(ka)-y(k)*x(ka)
                    ntc(k)=ilm1(ka)-ilm1(k)
                enddo    
                l2End = omp_get_wtime()
                l2Time = l2End - l2Start
                
		!print *, 'Loop 3...'
                !LOOP 3
                l3Start = omp_get_wtime()
                do k=1,npt
                    sig=sign(one,cz(k))
                    sq(k)=dq*sig
                    ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
                    if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
                        if (abs(cz(k)) .gt. zero) then
                            cx(k)=cx(k)/cz(k)
                            cy(k)=cy(k)/cz(k)
                        endif
                enddo
                l3End = omp_get_wtime()
                l3Time = l3End - l3Start
            !!$OMP END DO
        !!$OMP END PARALLEL

        !----------------------------------------------------------------------
        !Initialise PV jump array:
	!print *, 'Loop 4...'
        !LOOP 4
        l4Start = omp_get_wtime()        
        do i=1,ntf
            do j=0,ngf+1
                qa(j,i)=zero
                qa_jp1(j,i)=zero
            enddo
        enddo
        l4End = omp_get_wtime()
        l4Time = l4End - l4Start

        !Determine crossing indices:
	!LOOP 5
        l5Start = omp_get_wtime()
        !print *, 'Loop 5...'

    ! !split i's across threads
    !     !i's in range 1 to ntf
    ! !$OMP PARALLEL PRIVATE(threadID, numThreads, chunk, start, end)
    !     threadID=omp_get_thread_num()
    !     numThreads=omp_get_num_threads()
    !     !allocate(thread_is(0:numThreads-1, 2))

    !     !get even and contiguous i's for each thread using threadID and mod
    !     threadID=omp_get_thread_num() !to get in range 1 to numThreads

	! !$OMP CRITICAL
    !         print *, 'Thread ID: ', threadID
    !         print *, 'Number of threads: ', numThreads
    !     !$OMP END CRITICAL

    !     chunk = ntf/numThreads

    !     start=threadID*chunk + 1
    !     if (threadID .ne. numThreads-1) then
    !         end=start + chunk - 1
    !     else
    !         end=ntf
    !     endif
    !     !thread_is(threadID, 1)=start
    !     !thread_is(threadID, 2)=end
    
    !$OMP PARALLEL PRIVATE(k,j,i,ioff,ncr,rlatc,p,jump)
	    !!$OMP PARALLEL DO SCHEDULE(GUIDED)!, REDUCTION(+:qa), PRIVATE(k,j,i,ioff,ncr,rlatc,p,jump)
        !do k=1,npt
        !split k=1,npt into repeat and non-repeat ks. have master do repeat ks serially
        !$OMP MASTER
        do kk=1,size(callsRepeatKs)
            k=callsRepeatKs(kk)
            if (ntc(k) .ne. 0) then
                jump=sign(1,ntc(k))
                ioff=ntf+ilm1(k)+(1+jump)/2
                ncr=0
                do while (ncr .ne. ntc(k))
                    i=1+mod(ioff+ncr,ntf)
                    ncr=ncr+jump

                    rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                    j=int(rlatc)+1
                    p=rlatc-dble(j-1)
                    qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                    qa_jp1(j+1,i)=qa_jp1(j+1,i)+    p*sq(k)
                enddo
            endif
        enddo
        !$OMP END MASTER
        
        !$OMP DO
        do kk=1,size(nonRepeatKs)
            k=nonRepeatKs(kk)
        !if (mod(k,1) .eq. 0) then
                ! !$OMP CRITICAL
                ! !print *, 'Thread ', omp_get_thread_num(), ' at k=', k
                ! !$OMP END CRITICAL
            !endif
            if (ntc(k) .ne. 0) then
                jump=sign(1,ntc(k))
                ioff=ntf+ilm1(k)+(1+jump)/2
                ncr=0
                do while (ncr .ne. ntc(k))
                    i=1+mod(ioff+ncr,ntf)
                    ncr=ncr+jump
                    
                    !check if i in thread's range (start to end)
                    !if (i < start .or. i > end) cycle

                    rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                    j=int(rlatc)+1
                    p=rlatc-dble(j-1)
                    !!$OMP CRITICAL
                    !print *, 'Thread ', omp_get_thread_num(), 'at k=', k, 'i, j: ', i, ',', j
                    !!$OMP ATOMIC
                    qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                    !!$OMP ATOMIC
                    qa_jp1(j+1,i)=qa_jp1(j+1,i)+    p*sq(k)
                    !!$OMP END CRITICAL

                    !print *, 'i,j', i, ',', j
                enddo
            endif
        enddo
        !$OMP END DO

        !combine qa and qa_jp1 into qa
        
        ! !$omp do collapse(2)
        ! do j = 0, ngf+1
        !     do i = 1, ntf
        !         qa(j,i) = qa(j,i) + qa_jp1(j,i)
        !     end do
        ! end do
        ! !$omp end do

        !!$OMP END PARALLEL DO
    !$OMP END PARALLEL

        ! combineStart = omp_get_wtime()
        ! !$omp parallel do collapse(2)
        ! do j = 0, ngf+1
        !     do i = 1, ntf
        !         qa(j,i) = qa(j,i) + qa_jp1(j,i)
        !     end do
        ! end do
        ! !$omp end parallel do
        ! combineEnd = omp_get_wtime()
        ! combineTime = combineEnd - combineStart

        ! combineStart = omp_get_wtime()
        !combine qa and qa_jp1 into qa serially
        qa = qa + qa_jp1
        combineEnd = omp_get_wtime()
        ! combineTime = combineEnd - combineStart

        l5End = omp_get_wtime()
        l5Time = l5End - l5Start

        !Get PV values, at half latitudes, by sweeping through latitudes:
        !LOOP 6
        l6Start = omp_get_wtime()
        do i=1,ntf
            do j=2,ngf
                qa(j,i)=qa(j,i)+qa(j-1,i)
            enddo
        enddo
        l6End = omp_get_wtime()
        l6Time = l6End - l6Start
        !Here, qa(j,i) stands for the PV at latitude j-1/2,
        !from j = 1, ..., ngf.

        preAvgEnd = omp_get_wtime()
        preAvgTime = preAvgEnd - preAvgStart
    
        !----------------------------------------------------------------------
        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

        avgStart = omp_get_wtime()

        !Average PV values on the fine grid to get corresponding 
        !values on the inversion grid (ng,nt):
        ngh=ngf
        nth=ntf
        
        do while (ngh .gt. ng) !need to precalc number of iters if want to parallelise with OMP
            !Pre-store PV adjacent to poles at complementary longitudes (+pi):
            nthh=nth/2
            nghp1=ngh+1
            do i=1,nthh
                ic=i+nthh
                qa(0,i)=qa(1,ic)
                qa(0,ic)=qa(1,i)
                qa(nghp1,i)=qa(ngh,ic)
                qa(nghp1,ic)=qa(ngh,i)
            enddo

            !Work from SP to NP to define PV at full latitudes from averages
            !at adjacent half latitudes:
            do i=1,nth
                do j=0,ngh
                qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
                enddo
            enddo

            !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

            !Next 1-2-1 average these values to define PV at half latitudes
            !on a grid twice as coarse:
            nghh=ngh/2
            do i=1,nth
                do j=1,nghh
                je=2*j
                qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
                enddo
            enddo

            !Now perform analogous longitudinal 1-2-1 average:
            do j=1,nghh
                qaend(j)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=1,nth-1
                ip1=i+1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
                enddo
            enddo

            do j=1,nghh
                qa(j,nth)=qaend(j)
            enddo
            !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

            !Average these on the twice coarser grid:
            nthh=nth/2

            do j=1,nghh
                qa(j,1)=f12*(qa(j,nth)+qa(j,1))
            enddo

            do i=2,nthh
                io=2*i-1
                ie=io-1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,ie)+qa(j,io))
                enddo
            enddo

            ngh=nghh
            nth=nthh

        enddo

        !Finalise and take away f to define PV anomaly:
        do i=1,nt
            do j=1,ng
                qc(j,i)=qa(j,i)-fcor(j)
            enddo
        enddo

        avgEnd = omp_get_wtime()
        avgTime = avgEnd - avgStart

        endTime = omp_get_wtime()
        totalTime = endTime - startTime
        !print *, 'call', callCount, ' of con2grid took ', totalTime, ' seconds.'

        call accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, combineTime)
        
        return
        
    end subroutine 

    subroutine accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, combineTime)

        !passed args
        double precision:: totalTime, preAvgTime, avgTime
        double precision:: l1Time, l2Time, l3Time, l4Time, l5Time
        double precision:: l6Time
        double precision:: combineTime

        !accumulate times into total timers
        con2gridToTTime = con2gridToTTime + totalTime
        preAvgTotTime = preAvgTotTime + preAvgTime
        avgTotTime = avgTotTime + avgTime
        l1TotTime = l1TotTime + l1Time
        l2TotTime = l2TotTime + l2Time
        l3TotTime = l3TotTime + l3Time
        l4TotTime = l4TotTime + l4Time
        l5TotTime = l5TotTime + l5Time
        combineTotTime = combineTotTime + combineTime
        l6TotTime = l6TotTime + l6Time
        return
    end subroutine

    subroutine readOutputs

        integer:: current_pos
        inquire(unit=102, pos=current_pos)
        !print *, 'callCount: ', callCount
        !print *, 'Reading output at position: ', current_pos

        read(102) qc_arr(:,:, callCount)
        return
    end subroutine

    subroutine readOutputsReversed
        integer:: filesize
        integer:: pos
        
        inquire(unit=102, size=filesize)
        pos = 1 + filesize - qcSize*callCount

        !print *, 'File size of cgc_outputs.dat: ', filesize
        !print *, 'QC block size: ', qcSize
        !print *, 'Call count: ', callCount

        read(102, pos=pos) qc_arr(:,:, callCount)

        return
    end subroutine

    subroutine compare_qcs
        ! ensure that qc computed matches qc from file
        implicit double precision(a-h,o-z)
        implicit integer(i-n)
        double precision:: qc_file(ng,nt)
        double precision:: max_diff
        integer:: i,j
        
        qc_file = qc_arr(:,:, callCount)
        
        max_diff = 0.0d0
        do j=1,nt
            do i=1,ng
                if (abs(qc(i,j) - qc_file(i,j)) > max_diff) then
                    max_diff = abs(qc(i,j) - qc_file(i,j))
                endif
            enddo
        enddo

        qcDiffs((iterCount-1)*numInputs + callCount) = max_diff

        return
    end subroutine

    subroutine finish
        call printTimes
        call deallocateVars
        return
    end subroutine

    subroutine printTimes
        print *, 'Total con2grid time: ', con2gridToTTime
        print *, 'preAvg time', preAvgTotTime
        print *, 'Avg time', avgTotTime
        print *, 'Loop 1 total time: ', l1TotTime
        print *, 'Loop 2 total time: ', l2TotTime
        print *, 'Loop 3 total time: ', l3TotTime
        print *, 'Loop 4 total time: ', l4TotTime
        print *, 'Loop 5 total time: ', l5TotTime
        !print *, 'Combine time: ', combineTotTime
        print *, 'Loop 6 total time: ', l6TotTime
        return
    end subroutine

    subroutine deallocateVars
        deallocate(x_arr)
        deallocate(y_arr)
        deallocate(z_arr)
        deallocate(next_arr)
        deallocate(npt_arr)
        deallocate(qc_arr)
        return
    end subroutine

end program cgcDev
