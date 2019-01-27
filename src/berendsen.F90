module nvt_setup

real(8) :: tomega

end module 


!-----------------------------------------------------------------------
subroutine berendsen(atype, v)
use atoms; use parameters
use nvt_setup
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
real(8) :: Ekinetic, ctmp
real(8) :: ctmp0, eGKE

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

Ekinetic= eGKE(atype, v)

ctmp0 = (treq*UTEMP0)/( GKE*UTEMP )
ctmp = sqrt (1 + (dt*UTIME)/tomega * (ctmp0 -1))

do i = 1, NATOMS
  v(i, 1:3) = v(i, 1:3)* ctmp
enddo 

!call LinearMomentum(atype, v)

return
end


!-----------------------------------------------------------------------
function eGKE(atype, v)
use atoms; use parameters 
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
real(8) :: Ekinetic
real(8) :: buf, Gbuf, eGKE

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

buf = Ekinetic
call MPI_ALLREDUCE (buf, Gbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
eGKE= Gbuf

return 
end function 

!----------------------------------------------------------------------c
      subroutine myrnd(rnd,dseed)
      real:: rnd
      integer::dseed
!----------------------------------------------------------------------c
!  Random-number generator.
!----------------------------------------------------------------------c
      real*8 d2p31m,d2p31
      data d2p31m/2147483647d0/
      data d2p31 /2147483648d0/
      call srand(dseed) 
      rnd = rand()
      return
      end

