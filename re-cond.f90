program cond
implicit none
include 'mpif.h'

!===========================================================================
!   This code calculates conductivity of dwave superconductor using
! finite bandwidth. Imaginary part is obtained from Kramer-Kronig integral. 
!  Run dwaveborn0X.f90 first to obtain self energies.
! See PRB 72 035125 (2005) and my attachment for the equation.
! (check .dat file and make sure *maxw* in line 11 is the same)
!===========================================================================


integer :: rank, nump, ierr   ! these for mpi
integer, parameter :: king = 0  ! rank of the main processor
integer, parameter :: maxk = 64  ! # of k points (-maxk:maxk) ******
integer, parameter :: maxw = 1000 ! see dat file (-maxw:maxw)--------

double precision, parameter :: maxomg = 1.0d0 ! omega max (see dat file)***
double precision, parameter :: sig = 1.0e-10  ! small # to avoid NaN  
double precision, parameter :: sig2 = 1.0e-10
double precision, parameter :: temp = 10.0d0  ! temperature (see .dat)***
double precision, parameter :: kb = 8.617332385e-5
double precision, parameter :: charge = 1.602176e-19
double precision, parameter :: hbar = 6.582e-16
integer :: i, j, l, m, w, v, numk , ww, maxk1, maxk2

double precision, parameter :: pi = 3.141592654d0
double precision, parameter :: t0 = 0.0989d0
double precision, parameter :: t1 =-0.5908d0
double precision, parameter :: t2 = 0.0962d0
double precision, parameter :: t3 =-0.1306d0
double precision, parameter :: t4 =-0.0507d0
double precision, parameter :: t5 = 0.0939d0

double precision :: kx, ky, dk, deps2, nu, omg, dw, ferm1, ferm2, beta
double precision :: const, area
double complex :: g11, g12, g22, detg, gg11, gg12, gg22, detgg, func, func0
double complex ::  sigma, func1, func2
double precision :: lam, tau, pk, pre
double precision, dimension(-maxk:maxk,-maxk:maxk) :: eps
double precision,dimension(-maxw:maxw):: omega, phi1, phi2, kai1, kai2, Z1,Z2
double complex, dimension(-maxw:maxw) :: phi, kai, Z 
double complex, dimension(-maxw:maxw) :: sgm  ! real conductivity sigma
double complex, dimension(-maxw:maxw) :: sub_sgm  
double precision, dimension(-maxw:maxw) :: sgm2, sgm1

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
do w = -maxw, maxw
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
  dw = maxomg/dfloat(maxw)  ! domega in integral
  beta = 1.0d0/kb/temp
  pre = 1.0d0/dfloat(numk)

  ! constant proportionarity (dont need this actually)
  const = charge*charge/(2.0d0*hbar*hbar*pi) 


do v = 1,maxw  ! frequency nu for ploting

     sub_sgm(v) = dcmplx( 0.0d0, 0.0d0 )
     nu = dfloat(v)*dw

      if (rank==king) then
        print*, v
      end if

do i = -maxk+1, maxk ! sum over kx
      kx = dfloat(i)/dfloat(maxk)*pi 
do j = -maxk+1, maxk ! sum over ky
     ky = dfloat(j)/dfloat(maxk)*pi
do w = -maxw+rank, maxw, nump   ! integrate over w (omega)
 
   ww = v+w ! -v+w  is also ok since sigma1 is even function

   if ( abs(ww) <= maxw ) then

     dk = pi/dfloat(maxk)  ! dk
    deps2 = ( -0.5d0*t1*sin(kx) - t2*cos(ky)*sin(kx) - t3*sin(2.0d0*kx) &
        + t4*0.5d0*( -2.0d0*cos(ky)*sin(2.0d0*kx) - cos(2.0d0*ky)*sin(kx) ) &
        -2.0d0*t5*cos(2.0d0*ky)*sin(2.0d0*kx) )**2
 
     pk = 0.5d0*( cos(kx) - cos(ky) )
     
     omg = dfloat(w)*dw 
     g11 = Z(ww)+eps(i,j)+kai(ww)
     g12 = phi(ww)*pk
     g22 = Z(ww)-eps(i,j)-kai(ww)
     gg11 = Z(w)+eps(i,j)+kai(w)
     gg12 = phi(w)*pk
     gg22 = Z(w)-eps(i,j)-kai(w)
     detg = Z(ww)**2 - (eps(i,j)+kai(ww))**2 - (phi(ww)*pk)**2 + sig!--
     detgg = Z(w)**2 - (eps(i,j)+kai(w))**2 - (phi(w)*pk)**2 + sig !----
    
  !**********************************
     func0 =  (g11*gg11+2.0d0*g12*gg12+g22*gg22)/detg/detgg
     func1 = (g11*dconjg(gg11)+2.0d0*g12*dconjg(gg12)+g22*dconjg(gg22)) &
         /detg/dconjg(detgg)
     func2 = (dconjg(g11)*gg11+2.0d0*dconjg(g12)*gg12+dconjg(g22)*gg22)  &
       /dconjg(detg)/detgg


    ! func0 =  (g11*gg11+g12*gg12 )/detg/detgg
    ! func1 = (g11*dconjg(gg11)+ g12*dconjg(gg12) ) &
    !   /detg/dconjg(detgg)
    ! func2 = (dconjg(g11)*gg11+ dconjg(g12)*gg12 )  &
     !  /dconjg(detg)/detgg

    ! func0 =  (g22*gg22+g12*gg12 )/detg/detgg
    ! func1 = (g22*dconjg(gg22)+ g12*dconjg(gg12) ) &
    !   /detg/dconjg(detgg)
    ! func2 = (dconjg(g22)*gg22+ dconjg(g12)*gg12 )  &
    !   /dconjg(detg)/detgg


  !************************************

     func = -func0 + func1 + func2 - dconjg(func0)

     ferm1 = 1.0d0/( dexp(omg*beta) + 1.0d0 )
     ferm2 = 1.0d0/( dexp((dfloat(ww)*dw)*beta) + 1.0d0 ) 

    ! add sig(small parameter) to avoid singularity
     sub_sgm(v) = sub_sgm(v) + pre*dw*deps2*func*(ferm1-ferm2)/(nu + sig)

   end if

end do  ! end w
end do  ! end j
end do  ! end i

call mpi_barrier(mpi_comm_world,ierr)
call mpi_allreduce(sub_sgm(v), sgm(v), 1, mpi_double_complex, mpi_sum, &
        &  mpi_comm_world, ierr )
call mpi_barrier ( mpi_comm_world, ierr )

sgm1(v) = real(sgm(v))
sgm1(-v) = sgm1(v)  ! since sigma1 is even

end do  ! end v

!--------------------- End of real part conductivity ----------------


!******************************************************************
!        This needs to be adjusted until it agrees with the 
!       result obtained from im-cond.f90
 
           sgm1(0) = 2200.0d0  ! 1100 for 1mev grid

!       It can be approximated using least square to fit Lorentzian
!       since  sigma(0) ~ delta(0)
!******************************************************************


area = 0.0d0
do v = 0,maxw-1  ! integrating sigma (needed for lambda, tao)
     area = area + (sgm1(v+1)+sgm1(v))*dw/2.0d0
end do




!---------- Starting imaginary part of conductivity ---------------

if (rank==king) then
do v = -maxw,maxw
        sgm2(v) = 0.0d0
        nu = dfloat(v)*dw
do w = -maxw,maxw
        omg = dfloat(w)*dw
        sgm2(v) = sgm2(v) - 1.0d0/pi*dw*(sgm1(w)-sgm1(v))/(omg-nu+sig2)
end do
   sgm2(v) = sgm2(v) -sgm1(v)/pi*dlog( abs((maxomg-nu+sig2)/(maxomg+nu+sig2) ))
end do
end if


if (rank==king) then
open(21,file='con-original.OUT')
open(22,file='opt-original.OUT')
do v = -maxw,maxw
        nu = dfloat(v)*dw
        sigma = dcmplx( sgm1(v), sgm2(v) )
        
        lam = 2.0d0*area/pi*dimag( 1.0d0/sigma ) + nu
        tau = 2.0d0*area/pi*dreal( 1.0d0/sigma )

        write(21,*) nu,  dreal(sigma), dimag(sigma)
        write(22,*) nu, lam, -tau
end do
close(21)
close(22)
end if


call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)
stop
end program cond






