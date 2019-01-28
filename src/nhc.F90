!-----------------------------------------------------------------------
module nhc_setup
!-----------------------------------------------------------------------
integer, parameter :: nresnmax=5
integer, parameter :: nyoshmax=5
integer, parameter :: nnosmax=5
integer :: nresn, nyosh, nnos
real(8) :: qomega, gnkt, gkt
real(8),allocatable :: qmass(:), wdti(:)
real(8),allocatable :: xlogs(:), vlogs(:), glogs(:)
real(8), allocatable :: wdti2(:), wdti4(:), wdti8(:)
end module 

!-----------------------------------------------------------------------
subroutine nhc_init()
use nhc_setup; use atoms; use parameters; use nvt_setup
!-----------------------------------------------------------------------
implicit none 

integer:: i, j, nn
integer:: iyosh, dseed
real:: wdti_1, vsqr, fscale

allocate (xlogs(0:nnosmax), vlogs(0:nnosmax), glogs(0:nnosmax))
allocate (qmass(0:nnosmax))

nresn= 1; nyosh= 3; nnos= 3

  qomega = 1/tomega
  gnkt = 3.0*GNATOMS * treq
  if (myid==0) print *, treq*UTEMP0
  gkt = treq
  qmass(0) = gnkt / qomega**2
  do i = 1, nnos-1 
    qmass(i) = qmass(0)/ (3.0*dble(GNATOMS))
  enddo

!-----set parameters for Yoshida-Suzuki Integration

! wdti_1 = dt/float(nresn)
! if(nyosh.eq.1) then
!   wdti(1) = wdti_1
! elseif(nyosh.eq.3) then
!   wdti(1) = (1d0/(2d0-2d0**(1d0/3d0)))*wdti_1
!   wdti(2) = 1d0-2d0*wdti(1)
!   wdti(3) = wdti(1)
! elseif(nyosh.eq.5) then
!   wdti(1) = (1d0/(4d0-4d0**(1d0/3d0)))*wdti_1
!   wdti(2) = wdti(1)
!   wdti(3) = 1d0-4d0*wdti(1)
!   wdti(4) = wdti(1)
!   wdti(5) = wdti(1)
! else
!   print*,'unsupported nyosh selected--now quitting'
!   stop
! endif

! do iyosh = 1, nyosh
!     wdti2(iyosh)=0.5d0*wdti(iyosh)
!     wdti4(iyosh)=0.5d0*wdti2(iyosh)
!     wdti8(iyosh)=0.5d0*wdti4(iyosh)
! enddo 

     do i=0,nnos-1
       xlogs(i)=0d0
       vlogs(i)=0d0
!       glogs(i)=0d0
     enddo
        vlogs(nnos)=0d0


 do i = 1, nnos-1 
    glogs(i)= (qmass(i-1) * vlogs(i-1)*vlogs(i-1) - treq)/ qmass(i)
 enddo 

!if (myid==0) print *, glogs(1:2)
!     dseed=23219
!     vsqr=0d0
!     do i=1,nnos
!       call myrnd(vlogs(i),dseed)
!       vsqr=vsqr+qmass(i)*vlogs(i)*vlogs(i)
!       !!if (myid==0) print *, "printing KE fictional mass", vsqr
!     enddo
!     fscale=dsqrt(3d0*nnos*gkt/vsqr)
!     do i=1,nnos
!       vlogs(i)=vlogs(i)*fscale
!     enddo

end subroutine nhc_init

!-----------------------------------------------------------------------
subroutine nhc(atype, v)
use nhc_setup; use atoms; use parameters
use nvt_setup
!-----------------------------------------------------------------------
implicit none 
integer:: i, j, ii
integer:: nnos1, iresn, iyosh, inos
real(8) :: scale, Ekinetic, eGKE, eKE2
real(8) :: atype(NBUFFER), v(NBUFFER,3)
real(8)  :: aanhc, ncfac, t_current, factor_eta
real(8) :: dt8, dt4, dt2

dt2= 0.5*dt*UTIME
dt4= 0.5*dt2
dt8= 0.5*dt4


  qomega = 1/tomega
  gnkt = 3.0*GNATOMS * treq
  gkt = treq
  qmass(0) = gnkt / qomega**2
  do i = 1, nnos-1
    qmass(i) = qmass(0)/ (3.0*dble(GNATOMS))
  enddo


! calculate kinetic energy
if (qmass(0) >0) then 
  Ekinetic= eGKE(atype, v)
  t_current= Ekinetic*UTEMP/GNATOMS
  eKE2= 2.d0* Ekinetic
  glogs(0) = (eKE2 - gnkt)/ qmass(0) 
  else 
    glogs(0) = 0.0 
endif 



if (myid==0) print '(4f18.6)', glogs(0), qmass(0), eKE2, t_current


ncfac= 1.0/nresn

 do iresn = 0, nresn-1

!   do iyosh= 1, nyosh
    do inos= nnos-1, 1, -1 
      aanhc= exp(-ncfac*dt8*vlogs(inos+1))
      vlogs(inos) = vlogs(inos)*aanhc
      vlogs(inos) = vlogs(inos) + glogs(inos) * ncfac*dt4
      vlogs(inos) = vlogs(inos) * aanhc
      !if (myid==0) print '(2I4, 4f12.8)' , inos, iresn, aanhc, vlogs(inos), dt8, glogs(inos)
    enddo 
      
      aanhc = exp(-ncfac*dt8*vlogs(1))
      vlogs(0) = vlogs(0)*aanhc
      vlogs(0) = vlogs(0) + glogs(0)*ncfac*dt4
      vlogs(0) = vlogs(0)* aanhc

      factor_eta = exp(-ncfac*dt2*vlogs(0))

     ! if (myid==0) print *, factor_eta
      do i = 1, NATOMS
        v(i, 1:3) = v(i, 1:3)* factor_eta
      enddo 
      eKE2 = eKE2*factor_eta* factor_eta
      t_current= t_current* factor_eta* factor_eta
      if (qmass(0) >0) then
        glogs(0) = (eKE2 - gnkt)/ qmass(0)
        if (myid==0) print '(6f18.6)', glogs(0), eKE2, gnkt, qmass(0), t_current, factor_eta
      else 
          glogs(0) = 0.0
      endif 

    do inos= 0, nnos-1 
      xlogs(inos) = xlogs(inos) + ncfac* dt2 * vlogs(inos)
    enddo 

    vlogs(0) = vlogs(0)* aanhc
    vlogs(0) = vlogs(0) + glogs(0)*ncfac*dt4
    vlogs(0) = vlogs(0)* aanhc

  !if (myid==0) print '(4f12.6)', vlogs(0), aanhc, ncfac, glogs(0)

    do inos =1 , nnos-1 
      aanhc = exp(-ncfac*dt8*vlogs(inos+1))
      vlogs(inos) = vlogs(inos)*aanhc
      glogs(inos)= (qmass(inos-1) * vlogs(inos-1)*vlogs(inos-1) - treq)/ qmass(inos)
      vlogs(inos) = vlogs(inos) + glogs(inos) * ncfac*dt4
      vlogs(inos) = vlogs(inos) * aanhc
    enddo 

  enddo 

return
end subroutine nhc
