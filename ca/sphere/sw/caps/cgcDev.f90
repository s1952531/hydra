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

    call init
    
    !main loop
    do iterCount=1,numIters
        do callCount = 1, numInputs
            x = x_arr(:, callCount)
            y = y_arr(:, callCount)
            z = z_arr(:, callCount)
            next = next_arr(:, callCount)
            npt = npt_arr(callCount)
	
	    write(10, *) 'Sample Call:', callCount
            call getContourIndiceRange
            call con2grid(qc)
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
            call readOutputs
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
        integer:: qcSize=size(qc) * storage_size(qc)/8 !div by 8 to get bytes
        integer :: filesize
        
        inquire(unit=102, size=filesize)
        !print *, 'File size of cgc_outputs.dat: ', filesize
        !print *, 'Size of one qc array: ', qcSize
        numInputs=filesize/qcSize
        !print *, 'Number of inputs/outputs in files: ', numInputs

        return
    end subroutine


    subroutine openFiles
        open(10, file="loop5_k_ij.txt", action='write', status='replace')
        open(100, file="cgc_inputs.dat", status='old', action='read', access='stream', form='unformatted')
        open(102, file="cgc_outputs.dat", status='old', action='read', access='stream', form='unformatted')
        return
    end subroutine

    subroutine closeFiles
        close(10)
        close(100)
        close(102)
        return
    end subroutine


    subroutine readInput
        ! Reads in the contour data from a file "cgc_inputs.dat"
        implicit double precision(a-h,o-z)
        implicit integer(i-n)

        read(100) x_arr(:, callCount)
        read(100) y_arr(:, callCount)
        read(100) z_arr(:, callCount)
        read(100) next_arr(:, callCount)
        read(100) npt_arr(callCount)
        return
    end subroutine

    subroutine getContourIndiceRange
        integer:: k

        do k=1,npt
            print *, 'next(', k, ') = ', next(k)
        enddo
        return
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

        !Initialise crossing information:

        !!$OMP PARALLEL DEFAULT(NONE) SHARED(npt,dlfi,x,y,z,next,ilm1,cx,cy,cz,ntc,sq, zero) PRIVATE(k,ka,sig)
            !!$OMP DO SCHEDULE(STATIC)
                do k=1,npt
                    ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
                enddo
            !!$OMP END DO

            !!$OMP DO SCHEDULE(STATIC)
                do k=1,npt
                    ka=next(k)
                    cx(k)=z(k)*y(ka)-y(k)*z(ka)
                    cy(k)=x(k)*z(ka)-z(k)*x(ka)
                    cz(k)=x(k)*y(ka)-y(k)*x(ka)
                    ntc(k)=ilm1(ka)-ilm1(k)
                enddo    
                
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
            !!$OMP END DO
        !!$OMP END PARALLEL

        !----------------------------------------------------------------------
        !Initialise PV jump array:
        do i=1,ntf
            do j=0,ngf+1
                qa(j,i)=zero
            enddo
        enddo

        !Determine crossing indices:
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

                write(10,*) 'k=',k,' i=',i,' j=',j
                enddo
            endif
        enddo

        !Get PV values, at half latitudes, by sweeping through latitudes:
        do i=1,ntf
            do j=2,ngf
                qa(j,i)=qa(j,i)+qa(j-1,i)
            enddo
        enddo
        !Here, qa(j,i) stands for the PV at latitude j-1/2,
        !from j = 1, ..., ngf.

        !----------------------------------------------------------------------
        !Average PV values on the fine grid to get corresponding 
        !values on the inversion grid (ng,nt):
        ngh=ngf
        nth=ntf

        do while (ngh .gt. ng)
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
        return
    end subroutine

    subroutine readOutputs
        read(102) qc_arr(:,:, callCount)
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
        deallocate(qcDiffs)
        deallocate(x_arr)
        deallocate(y_arr)
        deallocate(z_arr)
        deallocate(next_arr)
        deallocate(npt_arr)
        deallocate(qc_arr)
        return
    end subroutine

end program cgcDev
