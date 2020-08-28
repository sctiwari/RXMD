module buck_param

real(8), allocatable :: TBL_buck_V(:,:,:), TBL_buck_F(:,:,:)
real(8):: UDR_buck, UDRi_buck
real(8), allocatable :: sigma(:,:), abuck(:,:)
real(8), parameter:: rcut_buck= 7.0
real(8), parameter:: rcut2_buck= 49.0
integer :: NTABLE_BUCK
end module buck_param


subroutine buck_init()
use parameters; use atoms
use buck_param

implicit none
integer :: i, j, ii
real(8) ::  dif
real(8) :: sixth = 1.0/6.0

NTABLE_BUCK = 5000

allocate(sigma(nso, nso), abuck(nso, nso))
allocate(TBL_buck_V(nso, nso, NTABLE_BUCK), TBL_buck_F(nso,nso,NTABLE_BUCK))

TBL_buck_V(:,:,:) = 0.0
TBL_buck_F(:,:,:) = 0.0

sigma(1,7) = 5.0; abuck(1,7) = 100.0
sigma(2,7) = 5.0; abuck(2,7) = 100.0
sigma(3,7) = 5.0; abuck(3,7) = 100.0
sigma(4,7) = 5.0; abuck(4,7) = 100.0
sigma(5,7) = 5.0; abuck(5,7) = 100.0
sigma(6,7) = 5.0; abuck(6,7) = 100.0


sigma(7,1) = 5.0; abuck(7,1) = 100.0
sigma(7,2) = 5.0; abuck(7,2) = 100.0
sigma(7,3) = 5.0; abuck(7,3) = 100.0
sigma(7,4) = 5.0; abuck(7,4) = 100.0
sigma(7,5) = 5.0; abuck(7,5) = 100.0
sigma(7,6) = 5.0; abuck(7,6) = 100.0

end subroutine buck_init


subroutine buck_potential_table()
use parameters
use buck_param
use atoms

implicit none 
real(8) :: dr2, dr1, temp1, temp2, temp3
real(8) :: ri, ri2, ri6, ri12, voffset, foffset
integer :: ity, jty, i
real(8) :: sig, buck, buck1, rexp, crexp


UDR_buck= rcut_buck*rcut_buck / dble(NTABLE_BUCK)
UDRi_buck= 1.d0/UDR_buck

open(10,file='buck_table.txt')


do ity = 1, 6
    jty= 7
    do i= 1, NTABLE_BUCK
      dr2 = UDR_buck*i
      dr1 = sqrt(dr2)
      sig= sigma(ity, jty) 
      buck= abuck(ity, jty) 
      buck1= buck/sig
      rexp= exp(-dr1/sig) 
      crexp= exp(-rcut_buck/sig) 
      TBL_buck_V(ity, jty, i) = buck*rexp - buck*crexp + (dr1- rcut_buck )* buck1*crexp
      TBL_buck_F(ity, jty, i) = buck1*rexp - buck1*crexp
      TBL_buck_V(jty, ity, i) = buck*rexp - buck*crexp + (dr1- rcut_buck )* buck1*crexp
      TBL_buck_F(jty, ity, i) = buck1*rexp - buck1*crexp 
      !print *, myid, ity, jty, i
    enddo     
enddo   

do i=1, ntable 
    if (myid==0) write (10, '(I10, 6F17.10)' ), i, TBL_buck_V(1,7,i), TBL_buck_F(1,7,i), &
             &   TBL_buck_V(2,7,i), TBL_buck_F(2,7,i), TBL_buck_V(3,7,i), TBL_buck_F(3,7,i)
enddo 
close(10) 
end subroutine buck_potential_table


subroutine buck_force (atype, pos, f, energy)
use parameters; use atoms
use buck_param; use MemoryAllocator

implicit none
real(8),intent(in) :: atype(NBUFFER)
real(8),intent(in) :: pos(NBUFFER,3)
real(8),intent(inout) :: f(NBUFFER,3)


real(8) :: vdummy(1,1) !-- dummy v for COPYATOM. it works as long as the array dimension matches

integer :: i, j,j1, k, m, ity, jty, ii, jj
integer :: l2g, iid, jid, itb, drtb
real(8) :: energy, rsq
real(8) :: dr(0:3), vr, fr, ff(1:3)

integer :: itype(NBUFFER) !-- integer part of atype
integer :: gtype(NBUFFER) !-- global ID from atype


do i=1, NBUFFER
   itype(i)=nint( atype(i))
enddo
do i=1, NBUFFER
   gtype(i)=l2g(atype(i))
enddo


   do i=1, NATOMS
        ity= itype(i)
         iid = gtype(i)
        if (ity ==7 ) cycle
        do j1= 1, nbplist(i, 0)
            j= nbplist(i, j1)
            jid = gtype(j)
            jty= itype(j)
            if (jty /=7) cycle
            dr(1:3) = pos(i,1:3) - pos(j,1:3)
            rsq = sum(dr(1:3)*dr(1:3))
            if ((rsq < rcut2_buck) ) then

               itb = int(rsq*UDRi_buck)
               drtb = rsq - itb*UDR_buck
               drtb = drtb*UDRi_buck

                  vr = TBL_buck_V(ity,jty,itb) - drtb*(TBL_buck_V(ity,jty,itb) - TBL_buck_V(ity,jty,itb+1))
                  fr = TBL_buck_F(ity,jty,itb) - drtb*(TBL_buck_F(ity,jty,itb) - TBL_buck_F(ity,jty,itb+1))

                  energy = energy + vr
                  !write(*,*) atype(i), drtb, itb, atype(j), vr, fr, TE_LJ

                  ff(1:3) = fr*dr(1:3)
                  f( i, 1:3) = f(i, 1:3) + ff(1:3)
                  f( j, 1:3) = f(j, 1:3) - ff(1:3)
                  !print *, i, j, itype(i), itype(j), f( i, 1:3), f( j, 1:3)
#ifdef stress
                  ia=i; ja=j
                  include 'stress'
#endif

            endif ! cutoff distance 
        enddo ! j1= 1, nbplist
   enddo !i=1, NATOMS 


end subroutine buck_force
