program cond
implicit none
include 'mpif.h'

!===========================================================================
!   This code calculates imaginary part of conductivity of dwave 
! using finite bandwidth.  Run dwaveborn0X.f90 first to obtain self energies.
! See PRB 72 035125 (2005) and my attachment for the equation.
! (check .dat file and make sure *maxw* in line 11 is the same)
!===========================================================================


integer :: rank, nump, ierr   ! these for mpi
integer, parameter :: king = 0  ! rank of the main processor
integer, parameter :: maxk = 32 ! # of k points (-maxk:maxk) 


!-------------------------------------------------------------------
integer, parameter :: maxw = 100 ! # of omega for plots ********
integer, parameter :: maxw2 = 250 ! see dat file (-maxw:maxw) for integral 

!  Note that   dnu=maxomg/maxw & dw=dw2=maxomg2/maxw2
!  dnu = dw must be true

double precision, parameter :: maxomg2 = 1.0d0 ! max w from .dat (integral)
double precision, parameter :: maxomg = 0.4d0 !  max w for plot ******
!------------------------------------------------------------------------


double precision, parameter :: sig = 1.0e-10  ! small # to avoid NaN  
double precision, parameter :: sig2 = 1.0e-10
double precision, parameter :: temp = 10.0d0  ! temperature (see .dat)***
double precision, parameter :: kb = 8.617332385e-5
double precision, parameter :: charge = 1.602176e-19
double precision, parameter :: hbar = 6.582e-16
integer :: i, j, l, m, w, v, numk , ww

double precision, parameter :: pi = 3.141592654d0
double precision, parameter :: t0 = 0.0989d0
double precision, parameter :: t1 =-0.5908d0
double precision, parameter :: t2 = 0.0962d0
double precision, parameter :: t3 =-0.1306d0
double precision, parameter :: t4 =-0.0507d0
double precision, parameter :: t5 = 0.0939d0

double precision :: kx, ky, dk, deps2, nu, omg, dw, nf1, nf2, beta
double precision :: const, area, omg2, ddeps, spec, dv, staff, omg3, nf3
double complex :: g11, g12, g22, detg, gg11, gg12, gg22, detgg
double complex ::  sigma, func1, func2, green, func0
double precision :: lam, tau, pk, pre, num1, num2, sub, func
double precision, dimension(-maxk:maxk,-maxk:maxk) :: eps
double precision,dimension(-maxw2:maxw2):: omega, phi1, phi2, kai1, kai2, Z1,Z2
double complex, dimension(-maxw2:maxw2) :: phi, kai, Z 
double complex, dimension(-maxw:maxw) :: sgm  
double precision, dimension(-maxw:maxw) ::  subterm1, subterm2, subterm3
double precision, dimension(-maxw:maxw) :: term1, term2,  term3
double precision, dimension(-maxw:maxw) ::  sgm1, rterm1, rterm2

!------------------------------------------------------------------------

call mpi_init(ierr)
call mpi_comm_rank (mpi_comm_world, rank, ierr)
call mpi_comm_size (mpi_comm_world, nump, ierr)

!constructing eps(i,j) band dispersion
do i = -maxk, maxk
do j = -maxk, maxk
        kx = dfloat(i)/dfloat(maxk)*pi
        ky = dfloat(j)/dfloat(maxk)*pi
   
   eps(i,j) = t0 + 0.5d0*t1*(cos(kx) + cos(ky)) + t2*cos(kx)*cos(ky) &
    &   + 0.5d0*t3*(cos(2.0d0*kx) + cos(2.0d0*ky)) &
    &   + 0.5d0*t4*(cos(2.0d0*kx)*cos(ky) + cos(kx)*cos(2.0d0*ky)) &
    &   + t5*cos(2.0d0*kx)*cos(2.0d0*ky)

end do 
end do  

! reading self energy from .dat file**************
open(11,file='phi-test.dat') !*************
open(12,file='kai-test.dat')  !*************
open(13,file='zbar-test.dat')  !*****************
do w = -maxw2, maxw2
      read(11,*) omega(w), phi1(w), phi2(w)
      read(12,*) omega(w), kai1(w), kai2(w)
      read(13,*) omega(w), Z1(w), Z2(w)
        phi(w) = dcmplx(phi1(w),phi2(w))
        kai(w) = dcmplx(kai1(w),kai2(w))
        Z(w) = dcmplx(Z1(w),Z2(w))
end do
close(11)
close(12)
close(13)



!------------------ Real part of condcutivity -------------------------

  numk = (2*maxk)*(2*maxk)  ! total # of k points
  dw = maxomg2/dfloat(maxw2)  ! domega in integral
  beta = 1.0d0/kb/temp
  pre = 1.0d0/dfloat(numk)
  dv = maxomg/dfloat(maxw)

     if (rank==0) then
        print*, ' dw = dv needs to be true '
        print*, 'dw = ', dw, ' and dv = ', dv
     endif


  ! constant proportionarity (dont need this actually)
  const = charge*charge/(2.0d0*hbar*hbar*pi) 



do v = 1, maxw  ! frequency nu for ploting **********
    
     nu = dfloat(v)*dv

   subterm1(v) = 0.0d0
   subterm2(v) = 0.0d0
   subterm3(v) = 0.0d0

do i = 0, maxk   ! -maxk+1,maxk
      kx = dfloat(i)/dfloat(maxk)*pi
do j = 0, maxk  ! -maxk+1,maxk
      ky = dfloat(j)/dfloat(maxk)*pi

         ! group velocity square
     deps2 = ( -0.5d0*t1*sin(kx) - t2*cos(ky)*sin(kx) - t3*sin(2.0d0*kx) &
        + t4*0.5d0*( -2.0d0*cos(ky)*sin(2.0d0*kx) - cos(2.0d0*ky)*sin(kx) ) &
        -2.0d0*t5*cos(2.0d0*ky)*sin(2.0d0*kx) )**2

         ! second derivative of eps(i,j)
     ddeps = -t1*0.5d0*cos(kx) -t2*cos(ky)*cos(kx) - t3*2.0d0*cos(2.0d0*kx) &
        - 0.5d0*t4*(4.0d0*cos(ky)*cos(2.0d0*kx) + cos(2.0d0*ky)*cos(kx)) &
        -4.0d0*t5*cos(2.0d0*ky)*cos(2.0d0*kx)

      pk = 0.5d0*( cos(kx) - cos(ky) )

do w = -maxw2+rank, maxw2, nump ! integrate over w (omega)

    if (  abs(w+v) <= maxw2 ) then

     omg = dw*dfloat(w)

     g11 = Z(w+v)+eps(i,j)+kai(w+v)
     g12 = phi(w+v)*pk
     g22 = Z(w+v)-eps(i,j)-kai(w+v)
     gg11 = Z(w)+eps(i,j)+kai(w)
     gg12 = phi(w)*pk
     gg22 = Z(w)-eps(i,j)-kai(w)
     detg = Z(w+v)**2 - (eps(i,j)+kai(w+v))**2 - (phi(w+v)*pk)**2   !---
     detgg = Z(w)**2 - (eps(i,j)+kai(w))**2 - (phi(w)*pk)**2   !----


 !*************************************************
     func0 =  (g11*gg11+2.0d0*g12*gg12+g22*gg22)/detg/detgg
     func1 = (g11*dconjg(gg11)+2.0d0*g12*dconjg(gg12)+g22*dconjg(gg22)) &
         /detg/dconjg(detgg)
     func2 = (dconjg(g11)*gg11+2.0d0*dconjg(g12)*gg12+dconjg(g22)*gg22)  &
       /dconjg(detg)/detgg

   !  func0 =  (g11*gg11+g12*gg12)/detg/detgg
   !  func1 = (g11*dconjg(gg11)+g12*dconjg(gg12)) &
   !      /detg/dconjg(detgg)
   !  func2 = (dconjg(g11)*gg11+dconjg(g12)*gg12)  &
   !    /dconjg(detg)/detgg

   ! func0 =  (g22*gg22+g12*gg12)/detg/detgg
   !  func1 = (g22*dconjg(gg22)+g12*dconjg(gg12)) &
   !      /detg/dconjg(detgg)
   !  func2 = (dconjg(g22)*gg22+dconjg(g12)*gg12)  &
   !    /dconjg(detg)/detgg


     sub = real(-func0 + func1 + func2 - dconjg(func0) )
 !****************************************************


     omg3 = nu+omg  
     nf1 = 1.0d0/( dexp(omg*beta) + 1.0d0 )
     nf3 = 1.0d0/( dexp(omg3*beta) + 1.0d0 )

     num2 = sub*(nf1-nf3)

     staff = abs(omg+nu+maxomg2+sig)/abs(omg+nu-maxomg2+sig)

     subterm2(v) = subterm2(v) + pre/pi*deps2*dw*num2/nu*dlog(staff)


    end if  !*******************


do ww = -maxw2, maxw2  ! integrate over w2
     omg2 = dw*dfloat(ww)  

     g11 = Z(ww)+eps(i,j)+kai(ww)
     g12 = phi(ww)*pk
     g22 = Z(ww)-eps(i,j)-kai(ww)
     gg11 = Z(w)+eps(i,j)+kai(w)
     gg12 = phi(w)*pk
     gg22 = Z(w)-eps(i,j)-kai(w)
     detg = Z(ww)**2 - (eps(i,j)+kai(ww))**2 - (phi(ww)*pk)**2   !---
     detgg = Z(w)**2 - (eps(i,j)+kai(w))**2 - (phi(w)*pk)**2   !----
    
 !*****************************************
     func0 =  (g11*gg11+2.0d0*g12*gg12+g22*gg22)/detg/detgg
     func1 = (g11*dconjg(gg11)+2.0d0*g12*dconjg(gg12)+g22*dconjg(gg22)) &
         /detg/dconjg(detgg)
     func2 = (dconjg(g11)*gg11+2.0d0*dconjg(g12)*gg12+dconjg(g22)*gg22)  &
       /dconjg(detg)/detgg

   !  func0 =  (g11*gg11+g12*gg12)/detg/detgg
   !  func1 = (g11*dconjg(gg11)+g12*dconjg(gg12)) &
   !      /detg/dconjg(detgg)
   !  func2 = (dconjg(g11)*gg11+dconjg(g12)*gg12)  &
   !    /dconjg(detg)/detgg


    !func0 =  (g22*gg22+g12*gg12)/detg/detgg
    ! func1 = (g22*dconjg(gg22)+g12*dconjg(gg12)) &
    !     /detg/dconjg(detgg)
    ! func2 = (dconjg(g22)*gg22+dconjg(g12)*gg12)  &
    !   /dconjg(detg)/detgg


     func = real( -func0 + func1 + func2 - dconjg(func0) )
 !*****************************************

     nf2 = 1.0d0/( dexp(omg2*beta) + 1.0d0 ) 
     num1 = func*(nf1-nf2)
 
     subterm1(v) =  subterm1(v) + pre/pi*deps2*dw*dw*(num1-num2)/nu &
        /(nu+omg-omg2+sig)


end do  ! end ww

 
    ! green = 2.0d0*Z(w) / &  !************
      green = (Z(w)   +eps(i,j)+kai(w))/ &  !********** 
    &  ( Z(w)**2 - (eps(i,j)+kai(w))**2 - (phi(w)*pk)**2 )

     spec = dimag(green)

     subterm3(v) = subterm3(v) + 4.0d0*pre*dw*nf1*spec*ddeps/nu

   !   end if  !************

end do  ! end w
end do  ! end j
end do  ! end i

call mpi_barrier(mpi_comm_world,ierr)
call mpi_allreduce(subterm1(v), term1(v), 1, mpi_double_precision, mpi_sum, &
        &  mpi_comm_world, ierr )
call mpi_allreduce(subterm2(v), term2(v), 1, mpi_double_precision, mpi_sum, &
        &  mpi_comm_world, ierr )
call mpi_allreduce(subterm3(v), term3(v), 1, mpi_double_precision, mpi_sum, &
        &  mpi_comm_world, ierr )
call mpi_barrier(mpi_comm_world,ierr)


if (rank==0) then
   print*, v, term1(v)+term2(v), term3(v)
end if

     rterm1(v) = term1(v) + term2(v)
     rterm2(v) = 2.0d0*term3(v)   ! two from spin 

     rterm1(-v) = -rterm1(v)
     rterm2(-v) = -rterm2(v)

end do  ! end v

!--------------------- End of real part conductivity ----------------


if (rank==0) then
open(21,file='im-cond.OUT')  !**********
do v = -maxw,maxw  
        nu = dfloat(v)*dv
     write(21,*) nu, 4.0d0*rterm1(v) - 4.0d0*rterm2(v)  
        ! 4.0 from even kx ky
end do
close(21)
end if


call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)
stop
end program cond






