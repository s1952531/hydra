!---------------------------------------------------------------------------
!      Computes a contour interval dq for an input field q from

!                    dq = (q_max - q_min)/n_q

!      dq is written to standard output, for use by the calling script.
!---------------------------------------------------------------------------

program contint
implicit none

integer,parameter:: ng=N_G, nt=2*ng
integer,parameter:: nq=N_Q
 !nq: number of contour intervals between q_min and q_max
integer,parameter:: ngridp=ng*nt,nbytes=4*(ngridp+1)

double precision:: qq(ng,nt),t,dq,qqmin,qqmax
integer:: i,j

 !Read initial gridded PV into qq:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(20,rec=1) t,qq
close(20)

qqmin=qq(1,1)
qqmax=qq(1,1)
do i=1,nt
  do j=1,ng
    qqmax=max(qqmax,qq(j,i))
    qqmin=min(qqmin,qq(j,i))
  enddo
enddo
dq=(qqmax-qqmin)/dble(nq)

write(*,'(1x,f16.12)') dq

end
