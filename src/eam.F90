module eam_param
integer :: nrho, nfrho
real(8) :: cutforcesq_eam, rhomax
real(8)::  drho, dfrho, cut_eam, rdr, rdrho
real(8),allocatable :: frho(:), rhor(:), z2r(:)
real(8),allocatable :: frho_spline(:,:), rhor_spline(:,:), z2r_spline(:,:)
end module


subroutine eam_init()
use parameters; use atoms
use eam_param; use MemoryAllocator

implicit none
integer :: i, j, ii
real(8) ::  dif
real(8) :: sixth = 1.0/6.0

!======================================================================================
!======================================================================================
!======================================================================================
open (unit=111, file='Pb_real.set', form='formatted', status='old')
read (111, *); read (111, *) ! skip tow line
read (111, *); read (111, *)

read (111, *) nrho, drho, nfrho, dfrho, cut_eam 
read (111, *) 

call allocatord1d (frho, 1, nrho)
call allocatord1d(rhor, 1, nfrho)
call allocatord1d(z2r, 1, nfrho)

i=1
do while ( i < 2000)
read (111, *) frho(i), frho(i+1), frho(i+2), frho(i+3), frho(i+4) 
!print *, i+4, frho(i), frho(i+1), frho(i+2), frho(i+3), frho(i+4)
i=i+5
enddo 

i=1
do while ( i < 2000)
read (111, *) rhor(i), rhor(i+1), rhor(i+2), rhor(i+3), rhor(i+4)
!print *, i+4, rhor(i), rhor(i+1), rhor(i+2), rhor(i+3), rhor(i+4)
i=i+5
enddo

rhomax= maxval(rhor(:))

i=1
do while ( i < 2000)
read (111, *) z2r(i), z2r(i+1), z2r(i+2), z2r(i+3), z2r(i+4)
i=i+5
enddo

i=0 
cutforcesq_eam= cut_eam* cut_eam


call array2spline()
end subroutine eam_init

subroutine array2spline()
use eam_param; use MemoryAllocator

implicit none 
integer:: m

rdr= 1./dfrho
rdrho = 1.0/drho

call allocatord2d (frho_spline, 1, nrho, 0, 6)
call allocatord2d (rhor_spline, 1,  nfrho, 0, 6)
call AllocatorD2D (z2r_spline, 1, nfrho, 0, 6)

call interpolate (nrho, drho, frho, frho_spline)
call interpolate (nfrho, dfrho, rhor, rhor_spline)
call interpolate (nfrho, dfrho, z2r, z2r_spline)

!do m =1549, 1549
!print *, m, rhor(m), rhor_spline(m,:)
!enddo
end subroutine array2spline

subroutine interpolate (n, delta, f, spline)
implicit none
integer:: n, m
real(8) :: delta
real(8) :: f(1:n)
real(8) :: spline (1:n, 0:6)

do  m =1, n
   spline(m, 6) = f(m)
enddo 

spline(1,5) = spline(2,6) -spline(1,6)
spline(2,5) = 0.5* (spline(3,6) -spline(1,6))
spline(n-1,5) = 0.5* (spline(n,6)- spline(n-2,6))
spline(n,5) = spline(n,6) -spline(n-1, 6)

do m=3, n-2 
  spline (m,5)= ((spline(m-2, 6)-spline(m+2,6)) + &
                    8.0*(spline(m+1,6)-spline(m-1,6))) / 12.0;
enddo 

do m =1, n-1 
  spline (m,4) = 3.0 *(spline(m+1,6) - spline(m,6)) -2.0*spline(m,5) -spline(m+1,5) 
  spline (m,3) = spline(m,5) + spline(m+1,5) -2.0*(spline(m+1,6) -spline(m,6)) 
enddo 

 spline (n,4) =0.0
 spline (n,3) =0.0

do m =1, n 
   spline (m,2) = spline(m,5)/ delta
   spline (m,1) = 2.0* spline (m,4) /delta
   spline (m,0) = 3.0* spline (m,3) /delta
enddo 

!do m =1549, 1549 
!print *, m, f(m), spline(m,:)
!enddo 
end subroutine

subroutine eam_force (atype, pos, f, energy)
use parameters; use atoms
use eam_param; use MemoryAllocator

implicit none 
real(8),intent(in) :: atype(NBUFFER)
real(8),intent(in) :: pos(NBUFFER,3)
real(8),intent(inout) :: f(NBUFFER,3)


real(8) :: vdummy(1,1) !-- dummy v for COPYATOM. it works as long as the array dimension matches

integer :: i, j,j1, k, m, ity, jty, ii, jj
integer :: l2g, counter, tempvar, iid, jid
real(8) ::embed_eam, energy
real(8) :: dr(3), ff(3)
real(8) ::  rsq,r,p,rhoip,rhojp,z2,z2p,recip,phip,psip,phi;
real(8) :: coeff(0:6)
real(8) ::  qdummy(1)
real(8) :: QCopyDr(3)
integer ::i_error


integer :: itype(NBUFFER) !-- integer part of atype
integer :: gtype(NBUFFER) !-- global ID from atype


real(8) :: rho(NBUFFER), fp(NBUFFER)


!print *, copyptr(6)

do i=1, NBUFFER
   itype(i)=nint( atype(i))
enddo
do i=1, NBUFFER
   gtype(i)=l2g(atype(i))
enddo

counter=0
ff(:)= 0.0d0
rho(:) =0.0d0 
embed_eam= 0.0

!print *, rhor_spline(1549, 0:6)

    do i = 1, NATOMS 
        ity= itype(i)
        if (ity /=7) cycle 
        
        do j1= 1, nbplist(i, 0) 
            j= nbplist(i, j1)
            jty= itype(i)
            if (jty /=7) cycle
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            rsq = sum(dr(1:3)*dr(1:3))
            if (rsq < cutforcesq_eam) then 
                p= (sqrt(rsq) * rdr) +1 
                m= int(p) 
                m = min (m, nfrho-1) 
                p = p -m
                p = min (p, 1.0) 
                coeff(0:6) = rhor_spline(m, 0:6)
                !if (i==1) print *, m, coeff(3:6) 
                counter= counter+1 
                rho(i) = rho(i) + coeff(3)*p + coeff(4)*p + coeff(5)*p + coeff(6)
            endif
        enddo  
    enddo ! do i =1, NATOMS  
!do i =1, natoms 
!  print *, i, rho(i)
!enddo 

    do i=1, NATOMS
      ity= itype(i) 
      if (ity /=7 ) cycle 
      p= rho(i)*rdrho +1.0; 
      m = int (p) 
      m = max (1, min(m, nrho-1))
      p = p- m
      coeff(0:6) = frho_spline(m, 0:6)
      fp(i) = (coeff(0)*p + coeff(1))*p +coeff(2)
      phi= (( coeff(3)*p + coeff(4))* p + coeff(5))*p +coeff(6)
      if (rho(i) > rhomax) phi = phi+ fp(i) * (rho(i)-rhomax);
      embed_eam= embed_eam+ phi
   enddo 

call COPYATOMS_EAM(MODE_EAM, NMINCELL*lcsize(1:3), atype, pos, vdummy, f, qdummy, fp)

   do i=1, NATOMS 
        ity= itype(i)
         iid = gtype(i)
        if (ity /=7 ) cycle
        do j1= 1, nbplist(i, 0)
            j= nbplist(i, j1)
            jid = gtype(j)
            jty= itype(j)
            if (jty /=7) cycle
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            rsq = sum(dr(1:3)*dr(1:3))
            if ((rsq < cutforcesq_eam) .and. (iid/=jid)) then
              r= sqrt(rsq)
              p= r*rdr+1 
              m= int(p)
              m= min(m, nfrho-1)
              p = p -m  
              p = min (p, 1.0)
              coeff(0:6) = rhor_spline(m, 0:6)
              rhoip = (coeff(0) *p + coeff(1))*p + coeff(2)
              rhojp = rhoip 
              coeff(0:6) = z2r_spline(m, 0:6) 
              z2p = (coeff(0) *p + coeff(1))*p + coeff(2)
              z2 = (( coeff(3)*p + coeff(4))* p + coeff(5))*p +coeff(6)
              recip = 1/r 
              phi = z2 * recip 
              phip= z2p * recip - phi*recip 
              psip = fp(i)*rhojp + fp(j)*rhoip + phip
              embed_eam= embed_eam+ 0.5*phi
              ff(1:3) = dr(1:3)* psip*recip 
              f(i,1:3) = f(i,1:3) -ff(1:3)
!--- stress calculation
#ifdef STRESS
                ia=i; ja=j
                include 'stress'
#endif

             ! if (fp(j)==0) print *, i, j,iid,jid ,rsq, psip, rhoip, phip, fp(i), fp(j)
            endif ! cutoff distance 
        enddo ! j1= 1, nbplist
   enddo !i=1, NATOMS 

!do i =1, NATOMS
!  if(myid==0) print *, f(1, 1:3)
!enddo 

energy= energy + embed_eam
!do i =1, NATOMS 
!   print *, i, f(i, 1:3)
!enddo
end subroutine eam_force 

