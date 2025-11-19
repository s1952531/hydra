program cgcDev
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
    integer:: next(0:npm),npt, callCount
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
    integer:: input_block_size= storage_size(x)/8 * size(x) + &
                                storage_size(y)/8 * size(y) + &
                                storage_size(z)/8 * size(z) + &
                                storage_size(next)/8 * size(next) + &
                                storage_size(npt)/8


    !Timers
    double precision:: con2gridToTTime=0.0d0, preAvgTotTime=0.0d0, avgTotTime=0.0d0
    double precision:: l1TotTime=0.0d0, l2TotTime=0.0d0, l3TotTime=0.0d0, l4TotTime=0.0d0, l5TotTime=0.0d0
    double precision:: l6TotTime=0.0d0, l7TotTime=0.0d0, l8TotTime=0.0d0, l9TotTime=0.0d0, l10TotTime=0.0d0
    double precision:: l11TotTime=0.0d0, l12TotTime=0.0d0, l13TotTime=0.0d0, l14TotTime=0.0d0, l15TotTime=0.0d0

    call init
    
    !main loop
    do iterCount=1,numIters
        do callCount = 1, numInputs
            x = x_arr(:, callCount)
            y = y_arr(:, callCount)
            z = z_arr(:, callCount)
            next = next_arr(:, callCount)
            npt = npt_arr(callCount)
            call con2grid(qc)
            call compare_qcs
        end do
    end do

    !print the max of the qcDiffs
    print *, 'Max difference in qc: ', maxval(qcDiffs)

    call finish

    contains

    subroutine init
        call openFiles
        call initVars
        
        do callCount = 1, numInputs
            !call readInput
            call readInputReversed
            !call readOutputs
            call readOutputsReversed
        end do
        call closeFiles
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
        print *, 'File size of cgc_outputs.dat: ', filesize
        print *, 'Size of one qc array: ', qcSize
        numInputs=filesize/qcSize
        print *, 'Number of inputs/outputs in files: ', numInputs

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


    subroutine readInput
        ! Reads in the contour data from a file "cgc_inputs.dat"

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

        integer:: filesize
        integer:: pos
        
        inquire(unit=100, size=filesize)

        pos = filesize - input_block_size*callCount

        read(100, pos=pos) x_arr(:, callCount)
        read(100) y_arr(:, callCount)
        read(100) z_arr(:, callCount)
        read(100) next_arr(:, callCount)
        read(100) npt_arr(callCount)

    end subroutine

    subroutine con2grid(qc)
        ! Calculates the PV anomaly field (stored in qc) from the PV 
        ! contours (x,y,z).  Takes away Coriolis frequency (fcor).

        implicit double precision(a-h,o-z)
        implicit integer(i-n)

        !Passed arrays:
        double precision:: qc(ng,nt)
        !Local arrays:
        double precision:: qa(0:ngf+1,ntf)
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

        double precision:: preAvgStart, preAvgEnd, preAvgTime
        double precision:: avgStart, avgEnd, avgTime

        !in while loop so accumulate
        double precision:: l7Start, l7End, l7Time=0.0d0
        double precision:: l8Start, l8End, l8Time=0.0d0
        double precision:: l9Start, l9End, l9Time=0.0d0
        double precision:: l10Start, l10End, l10Time=0.0d0
        double precision:: l11Start, l11End, l11Time=0.0d0
        double precision:: l12Start, l12End, l12Time=0.0d0
        double precision:: l13Start, l13End, l13Time=0.0d0
        double precision:: l14Start, l14End, l14Time=0.0d0
        !

        double precision:: l15Start, l15End, l15Time
        ! %      cumulative self     calls    
        ! 7.32    763.92    84.85    23810  __contours_MOD_con2grid
        ! 1.83   1067.81    21.25    23810  __contours_MOD_con2grid_avg

        startTime = omp_get_wtime()
        preAvgStart = startTime

        !Initialise crossing information:

        !!$OMP PARALLEL DEFAULT(NONE) SHARED(npt,dlfi,x,y,z,next,ilm1,cx,cy,cz,ntc,sq, zero) PRIVATE(k,ka,sig)
            !!$OMP DO SCHEDULE(STATIC)
                !LOOP 1
                l1Start = omp_get_wtime()
                do k=1,npt
                    ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
                enddo
                l1End = omp_get_wtime()
                l1Time = l1End - l1Start
            !!$OMP END DO

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
        !LOOP 4
        l4Start = omp_get_wtime()        
        do i=1,ntf
            do j=0,ngf+1
                qa(j,i)=zero
            enddo
        enddo
        l4End = omp_get_wtime()
        l4Time = l4End - l4Start

        !Determine crossing indices:
        !LOOP 5
        l5Start = omp_get_wtime()
        do k=1,npt
            if (ntc(k) .ne. 0) then
                jump=sign(1,ntc(k))
                ioff=ntf+ilm1(k)+(1+jump)/2
                ncr=0
                do while (ncr .ne. ntc(k))
                i=1+mod(ioff+ncr,ntf)
                rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
                j=int(rlatc)+1
                p=rlatc-dble(j-1)
                qa(j,i)=  qa(j,i)+(one-p)*sq(k)
                qa(j+1,i)=qa(j+1,i)+    p*sq(k)
                ncr=ncr+jump
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
            !LOOP 7
            l7Start = omp_get_wtime()
            do i=1,nthh
                ic=i+nthh
                qa(0,i)=qa(1,ic)
                qa(0,ic)=qa(1,i)
                qa(nghp1,i)=qa(ngh,ic)
                qa(nghp1,ic)=qa(ngh,i)
            enddo
            l7End = omp_get_wtime()
            l7Time = l7Time + l7End - l7Start

            !Work from SP to NP to define PV at full latitudes from averages
            !at adjacent half latitudes:
            !LOOP 8
            l8Start = omp_get_wtime()
            do i=1,nth
                do j=0,ngh
                qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
                enddo
            enddo
            l8End = omp_get_wtime()
            l8Time = l8Time + l8End - l8Start

            !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

            !Next 1-2-1 average these values to define PV at half latitudes
            !on a grid twice as coarse:
            nghh=ngh/2
            !LOOP 9
            l9Start = omp_get_wtime()
            do i=1,nth
                do j=1,nghh
                je=2*j
                qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
                enddo
            enddo
            l9End = omp_get_wtime()
            l9Time = l9Time + l9End - l9Start

            !Now perform analogous longitudinal 1-2-1 average:
            !LOOP 10
            l10Start = omp_get_wtime()
            do j=1,nghh
                qaend(j)=f12*(qa(j,nth)+qa(j,1))
            enddo
            l10End = omp_get_wtime()
            l10Time = l10Time + l10End - l10Start

            !LOOP 11
            l11Start = omp_get_wtime()
            do i=1,nth-1
                ip1=i+1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
                enddo
            enddo
            l11End = omp_get_wtime()
            l11Time = l11Time + l11End - l11Start

            !LOOP 12
            l12Start = omp_get_wtime()
            do j=1,nghh
                qa(j,nth)=qaend(j)
            enddo
            l12End = omp_get_wtime()
            l12Time = l12Time + l12End - l12Start
            !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

            !Average these on the twice coarser grid:
            nthh=nth/2

            !LOOP 13
            l13Start = omp_get_wtime()
            do j=1,nghh
                qa(j,1)=f12*(qa(j,nth)+qa(j,1))
            enddo
            l13End = omp_get_wtime()
            l13Time = l13Time + l13End - l13Start

            !LOOP 14
            l14Start = omp_get_wtime()
            do i=2,nthh
                io=2*i-1
                ie=io-1
                do j=1,nghh
                qa(j,i)=f12*(qa(j,ie)+qa(j,io))
                enddo
            enddo
            l14End = omp_get_wtime()
            l14Time = l14Time + l14End - l14Start

            ngh=nghh
            nth=nthh

        enddo

        !Finalise and take away f to define PV anomaly:
        !LOOP 15
        l15Start = omp_get_wtime()
        do i=1,nt
            do j=1,ng
                qc(j,i)=qa(j,i)-fcor(j)
            enddo
        enddo
        l15End = omp_get_wtime()
        l15Time = l15End - l15Start

        avgEnd = omp_get_wtime()
        avgTime = avgEnd - avgStart

        endTime = omp_get_wtime()
        totalTime = endTime - startTime
        print *, 'call', callCount, ' of con2grid took ', totalTime, ' seconds.'

        call accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, l7Time, l8Time, l9Time, l10Time, &
                              l11Time, l12Time, l13Time, l14Time, l15Time)
        
        return
        
    end subroutine 

    subroutine accumulateTimes(totalTime, preAvgTime, avgTime, &
                              l1Time, l2Time, l3Time, l4Time, l5Time, &
                              l6Time, l7Time, l8Time, l9Time, l10Time, &
                              l11Time, l12Time, l13Time, l14Time, l15Time)

        !passed args
        double precision:: totalTime, preAvgTime, avgTime
        double precision:: l1Time, l2Time, l3Time, l4Time, l5Time
        double precision:: l6Time, l7Time, l8Time, l9Time, l10Time
        double precision:: l11Time, l12Time, l13Time, l14Time, l15Time

        !accumulate times into total timers
        con2gridToTTime = con2gridToTTime + totalTime
        preAvgTotTime = preAvgTotTime + preAvgTime
        avgTotTime = avgTotTime + avgTime
        l1TotTime = l1TotTime + l1Time
        l2TotTime = l2TotTime + l2Time
        l3TotTime = l3TotTime + l3Time
        l4TotTime = l4TotTime + l4Time
        l5TotTime = l5TotTime + l5Time
        l6TotTime = l6TotTime + l6Time
        l7TotTime = l7TotTime + l7Time
        l8TotTime = l8TotTime + l8Time
        l9TotTime = l9TotTime + l9Time
        l10TotTime = l10TotTime + l10Time
        l11TotTime = l11TotTime + l11Time
        l12TotTime = l12TotTime + l12Time
        l13TotTime = l13TotTime + l13Time
        l14TotTime = l14TotTime + l14Time
        l15TotTime = l15TotTime + l15Time
        return
    end subroutine

    subroutine readOutputs
        read(102) qc_arr(:,:, callCount)
        return
    end subroutine

    subroutine readOutputsReversed
        integer:: filesize
        integer:: pos
        
        inquire(unit=102, size=filesize)
        pos = filesize - qcSize*callCount

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
        print *, 'Loop 6 total time: ', l6TotTime
        print *, 'Loop 7 total time: ', l7TotTime
        print *, 'Loop 8 total time: ', l8TotTime
        print *, 'Loop 9 total time: ', l9TotTime
        print *, 'Loop 10 total time: ', l10TotTime
        print *, 'Loop 11 total time: ', l11TotTime
        print *, 'Loop 12 total time: ', l12TotTime
        print *, 'Loop 13 total time: ', l13TotTime
        print *, 'Loop 14 total time: ', l14TotTime
        print *, 'Loop 15 total time: ', l15TotTime
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