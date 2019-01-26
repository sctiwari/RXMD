subroutine berendsen(atype, v)
use atoms; use parameters
!-----------------------------------------------------------------------
implicit none
real(8) :: atype(NBUFFER), v(NBUFFER,3)

integer :: i,ity
real(8) :: Ekinetic, ctmp
real(8) :: buf, Gbuf, ctmp0

do i=1, NATOMS
   ity=nint(atype(i))
   Ekinetic=0.5d0*mass(ity)*sum(v(i,1:3)*v(i,1:3))
enddo

buf = Ekinetic
call MPI_ALLREDUCE (buf, Gbuf, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
Ekinetic= Gbuf

ctmp0 = (treq*UTEMP0)/( GKE*UTEMP )
ctmp = sqrt (1 + (dt*UTIME)/10.0 * (ctmp0 -1))

do i = 1, NATOMS
  v(i, 1:3) = v(i, 1:3)* ctmp
enddo 

!call LinearMomentum(atype, v)

return
end
