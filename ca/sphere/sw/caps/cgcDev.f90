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
    integer:: i, j
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

    double precision:: qc(ng,nt), qcDiffs(23810) !there are 23810 calls for the default test case

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

    !!!!!!!!

    do callCount = 1, 23810
        call readInput
        call con2grid(qc)
        call compare_qcs
    end do

    !print the max of the qcDiffs
    print *, 'Max difference in qc: ', maxval(qcDiffs)

    contains

    subroutine readInput
        ! Reads in the contour data from a file "cgc_inputs.dat"
        implicit double precision(a-h,o-z)
        implicit integer(i-n)

        open(100, file="cgc_inputs_100.dat", status='old', action='read', access='stream', form='unformatted')
        read(100) x
        read(100) y
        read(100) z
        read(100) next
        read(100) npt
        close(100)

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
        do k=1,npt
        ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
        enddo

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

    subroutine compare_qcs
        ! ensure that qc computed matches qc from file
        implicit double precision(a-h,o-z)
        implicit integer(i-n)
        double precision:: qc_file(ng,nt)
        double precision:: max_diff
        integer:: i,j
        open(102, file="cgc_outputs_100.dat", status='old', action='read', access='stream', form='unformatted')
        read(102) qc_file
        close(102)

        max_diff = 0.0d0
        do j=1,nt
            do i=1,ng
                if (abs(qc(i,j) - qc_file(i,j)) > max_diff) then
                    max_diff = abs(qc(i,j) - qc_file(i,j))
                endif
            enddo
        enddo

        qcDiffs(callCount) = max_diff

        return
    end subroutine

end program cgcDev