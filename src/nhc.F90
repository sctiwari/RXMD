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
real(8) :: tomega=100.d0
end module 

!-----------------------------------------------------------------------
subroutine nhc_init()
use nhc_setup; use atoms; use parameters
!-----------------------------------------------------------------------
implicit none 

integer:: i, j, nn
integer:: iyosh, dseed
real:: wdti_1, vsqr, fscale

allocate (xlogs(nnosmax), vlogs(nnosmax), glogs(nnosmax))
allocate (qmass(nnosmax))
allocate (wdti(nyoshmax), wdti2(nyoshmax), wdti4(nyoshmax), wdti8(nyoshmax))

nresn= 3; nyosh= 3; nnos= 3

  qomega = UTIME/tomega
  gnkt = 2.0*GNATOMS * treq/UTEMP
  gkt = 2.0*treq/UTEMP 
  qmass(1) = gnkt / qomega**2
  do i = 2, nnos 
    qmass(i) = qmass(1)/ dble(GNATOMS)
  enddo

!-----set parameters for Yoshida-Suzuki Integration

  wdti_1 = dt/float(nresn)
  if(nyosh.eq.1) then
    wdti(1) = wdti_1
  elseif(nyosh.eq.3) then
    wdti(1) = (1d0/(2d0-2d0**(1d0/3d0)))*wdti_1
    wdti(2) = 1d0-2d0*wdti(1)
    wdti(3) = wdti(1)
  elseif(nyosh.eq.5) then
    wdti(1) = (1d0/(4d0-4d0**(1d0/3d0)))*wdti_1
    wdti(2) = wdti(1)
    wdti(3) = 1d0-4d0*wdti(1)
    wdti(4) = wdti(1)
    wdti(5) = wdti(1)
  else
    print*,'unsupported nyosh selected--now quitting'
    stop
  endif

  do iyosh = 1, nyosh
      wdti2(iyosh)=0.5d0*wdti(iyosh)
      wdti4(iyosh)=0.5d0*wdti2(iyosh)
      wdti8(iyosh)=0.5d0*wdti4(iyosh)
  enddo 

      do i=1,nnos
        xlogs(i)=0d0
        vlogs(i)=0d0
        glogs(i)=0d0
      enddo
      dseed=23219d0
      vsqr=0d0
      do i=1,nnos
        call myrnd(vlogs(i),dseed)
        vsqr=vsqr+qmass(i)*vlogs(i)*vlogs(i)
        !!if (myid==0) print *, "printing KE fictional mass", vsqr
      enddo
      fscale=dsqrt(3d0*nnos*gkt/vsqr)
      do i=1,nnos
        vlogs(i)=vlogs(i)*fscale
      enddo

end subroutine nhc_init

!-----------------------------------------------------------------------
subroutine nhc(atype, v)
use nhc_setup; use atoms; use parameters
!-----------------------------------------------------------------------
implicit none 
integer:: i, j, ii
integer:: nnos1, iresn, iyosh, inos
real(8) :: scale, Ekinetic, eGKE, eKE2
real(8) :: atype(NBUFFER), v(NBUFFER,3)
real(8)  :: aanhc

  nnos1= nnos+1 
  scale = 1d0

! calculate kinetic energy
  Ekinetic= eGKE(atype, v)
  eKE2= 2.d0* Ekinetic
  glogs(1) = (eKE2 - gnkt)/ qmass(1) 
  
  do i = 1, nnos-1
    glogs(i+1) = (qmass(i) * vlogs(i) * vlogs(i) - gkt)/ qmass(i+1)
  enddo 

!if (myid==0) print *, gnkt, eKE2, qmass(1:3), glogs(1:3), vlogs(1:3)

  do iresn = 1, nresn

    do iyosh= 1, nyosh

!---------First half update of the thermostat velocities
        vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)
          do inos=1,nnos-1
            aanhc=dexp(-wdti8(iyosh)*vlogs(nnos1-inos))
            vlogs(nnos-inos)=( vlogs(nnos-inos)*aanhc + wdti4(iyosh)*glogs(nnos-inos) )*aanhc
          enddo


!---------Update the atom-velocity scaling factor
          aanhc=dexp(-wdti2(iyosh)*vlogs(1))
          scale=scale*aanhc

!---------Update the thermostat positions
          do inos=1,nnos
            xlogs(inos)=xlogs(inos)+vlogs(inos)*wdti2(iyosh)
          enddo

!---------Second half update of the thermostat velocities
          glogs(1)=(scale*scale*eKE2 - gnkt)/qmass(1)
          do inos=1,nnos-1
            aanhc=dexp(-wdti8(iyosh)*vlogs(inos+1))
            vlogs(inos)= ( vlogs(inos)*aanhc  + wdti4(iyosh)*glogs(inos) )*aanhc
            glogs(inos+1)=(qmass(inos)*vlogs(inos)*vlogs(inos)-gkt) /qmass(inos+1)
          enddo
          vlogs(nnos)=vlogs(nnos)+glogs(nnos)*wdti4(iyosh)

    enddo 
  enddo 
      v(1:natoms, 1:3) = v(1:natoms, 1:3)* scale

end subroutine nhc
