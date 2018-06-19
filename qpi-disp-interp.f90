program main
implicit none
include 'mpif.h'

! Note:  this is the qpi using interpolated data from matlab code


!=========================================================================
! The following code will calculate QPI density of states
! for mono-layer FeSe.
! Obtain self energy from matlab code
! see Wang etc, PRB, 67, 020511, 2003
!=========================================================================

!  NOTE:  values are interporated usng matlab
integer, parameter :: numw = 160   ! # of omega (-numw:numw)  **********
integer, parameter :: numk = 510  !  # of k  (-numk-1:numk-1) *********
!integer, parameter :: freq = 130  ! type energy of your interest [meV] *****
!  above  assum grid size 1meV

!check matlab kx(2)-kx(1)
!double precision, parameter :: dk = 0.012271846303085d0  


double precision, parameter :: imp = 0.01d0  ! impurity potential [eV]
double precision, parameter :: t0 = -1.4878d0  ! chemical potential ******

double precision, parameter :: t1 = 0.38d0  ! default positive *****


integer :: rank, nump, ierr
double precision ::  omg, kx, ky,  cee, spec, totalk, dk
double complex :: pdet, pgee0, pgee3, pgee1, mdet, mgee0, mgee3, mgee1


double precision, allocatable, dimension(:,:,:) ::  phi1, phi2, kai1, kai2
double precision, allocatable, dimension(:,:,:) ::  zbar1, zbar2
double complex, allocatable, dimension(:,:,:) :: phi, kai, zbar
double precision, allocatable, dimension(:,:) :: eps
double complex, allocatable, dimension(:,:,:,:) :: pgren0, mgren0
double complex, dimension(2,2) :: ptmat, mtmat      ! T-matrix 
double complex, dimension(2,2) :: pgren1, mgren1, pgren2, mgren2
double precision, allocatable, dimension(:,:) :: dos, subdos

double precision, parameter :: pi = 3.141592654d0

integer :: i, j, l, m, w , k1, k2, kk1, kk2, q1, q2, p1, p2, kk3, kk4, p3, p4


!--------------------------------------------------------------------------

call mpi_init( ierr)
call mpi_comm_rank(mpi_comm_world,rank,ierr)
call mpi_comm_size(mpi_comm_world,nump,ierr)

cee = 1.0d0/imp
totalk = dfloat( (2*numk+1)*(2*numk+1) )
dk = 0.012271846303085d0/2.0d0  ! for 128x4 grids


allocate( phi1(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( phi2(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( kai1(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( kai2(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( zbar1(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( zbar2(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( phi(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( kai(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )
allocate( zbar(-numk-1:numk-1,-numk-1:numk-1,-numw:numw)    )

if (rank==0)then
print*, 'begin reading file'
endif


open(11,file='self-interp-lam81x4.dat')  !*************
do i = -numk-1,numk-1
do j = -numk-1,numk-1
do w = -numw,numw
       read(11,*) kx, ky, omg, &
        zbar1(i,j,w), zbar2(i,j,w), &
        kai1(i,j,w), kai2(i,j,w), phi1(i,j,w), phi2(i,j,w)

     phi(i,j,w) = dcmplx( phi1(i,j,w), phi2(i,j,w) )
     kai(i,j,w) = dcmplx( kai1(i,j,w), kai2(i,j,w) )
     zbar(i,j,w) = dcmplx( zbar1(i,j,w), zbar2(i,j,w) ) 
end do
end do
end do 
close(11)
deallocate ( phi1, phi2, kai1, kai2, zbar1, zbar2 )

if (rank==0)then
print*, 'end reading file'
endif


!  constructing band dispersion
allocate( eps(-numk-1:numk-1,-numk-1:numk-1) )
do i = -numk-1,numk-1
do j = -numk-1,numk-1
      kx = dfloat(i)*dk
      ky = dfloat(j)*dk

        eps(i,j) = -2.0d0*t1*( cos(kx)+cos(ky) ) - t0

end do
end do


allocate( pgren0(-numk-1:numk-1,-numk-1:numk-1,2,2) )
allocate( mgren0(-numk-1:numk-1,-numk-1:numk-1,2,2) )

allocate( dos(-numk-1:numk-1,-numw:numw)   )
allocate( subdos(-numk-1:numk-1,-numw:numw)  )



!if (rank==0) then
!open(31,file='spec.dat')
!do k1 = -numk,numk
!     kx = dfloat(k1)*dk
!do w = -numw,numw

 !    pdet = zbar(k1,k1,w)**2 - (eps(k1,k1)+kai(k1,k1,w))**2 &
  !           - phi(k1,k1,w)**2

   ! spec = -1.0d0/pi*dimag( (zbar(k1,k1,w)+eps(k1,k1)+kai(k1,k1,w))/pdet  ) 
           
   !  write(31,*)   kx,  dfloat(w),  spec

!end do
!end do
!close(31)
!end if

!if (rank==0) then
!open(32,file='spec2.dat')
!do w=-numw,numw

 !       pdet = zbar(0,0,w)**2 - (eps(0,0)+kai(0,0,w))**2 &
 !            - phi(0,0,w)**2

  !  spec = -1.0d0/pi*dimag( (zbar(0,0,w)+eps(0,0)+kai(0,0,w))/pdet  )

  !   write(32,*)     dfloat(w),  spec
!end do
!close(32)
!end if



call mpi_barrier(mpi_comm_world,ierr)

!-----------------------------------------------------------------

do w = -numw+rank,numw,nump

print*, w

! construcitng greens0 function matrix
! p stands for positive freqency, m for minus
do i = -numk-1,numk-1
do j = -numk-1,numk-1

     pdet = zbar(i,j,w)**2 - (eps(i,j)+kai(i,j,w))**2 &
             - phi(i,j,w)**2
     
     pgren0(i,j,1,1) = ( zbar(i,j,w)+eps(i,j)+kai(i,j,w) )/pdet
     pgren0(i,j,1,2) = phi(i,j,w)/pdet
     pgren0(i,j,2,1) = pgren0(i,j,1,2)
     pgren0(i,j,2,2) = ( zbar(i,j,w)-eps(i,j)-kai(i,j,w) )/pdet

     mdet = zbar(i,j,-w)**2 - (eps(i,j)+kai(i,j,-w))**2 &
             - phi(i,j,-w)**2

     mgren0(i,j,1,1) = ( zbar(i,j,-w)+eps(i,j)+kai(i,j,-w) )/mdet
     mgren0(i,j,1,2) = phi(i,j,-w)/mdet
     mgren0(i,j,2,1) = mgren0(i,j,1,2)
     mgren0(i,j,2,2) = ( zbar(i,j,-w)-eps(i,j)-kai(i,j,-w) )/mdet

end do 
end do



! summing components of greens function for T-matrix

pgee0 = dcmplx(0.0d0,0.0d0)
pgee3 = dcmplx(0.0d0,0.0d0)
pgee1 = dcmplx(0.0d0,0.0d0)
mgee0 = dcmplx(0.0d0,0.0d0)
mgee3 = dcmplx(0.0d0,0.0d0)
mgee1 = dcmplx(0.0d0,0.0d0)
do i = -numk-1,numk-1
do j = -numk-1,numk-1

    pdet = zbar(i,j,w)**2 - (eps(i,j)+kai(i,j,w))**2  &
            - phi(i,j,w)**2

    pgee0 = pgee0 + 1.0d0/totalk*zbar(i,j,w)/pdet
    pgee3 = pgee3 + 1.0d0/totalk*(eps(i,j)+kai(i,j,w))/pdet
    pgee1 = pgee1 + 1.0d0/totalk*phi(i,j,w)/pdet

     mdet = zbar(i,j,-w)**2 - (eps(i,j)+kai(i,j,-w))**2  &
              - phi(i,j,-w)**2

    mgee0 = mgee0 + 1.0d0/totalk*zbar(i,j,-w)/mdet
    mgee3 = mgee3 + 1.0d0/totalk*(eps(i,j)+kai(i,j,-w))/mdet
    mgee1 = mgee1 + 1.0d0/totalk*phi(i,j,-w)/mdet

end do
end do



! T-matrix,  see Hussey 
ptmat(1,1) = (pgee0+cee-pgee3)/( (cee-pgee3)**2 - pgee0*pgee0 + pgee1*pgee1 )
ptmat(1,2) = pgee1/(  (cee-pgee3)**2 - pgee0*pgee0 + pgee1*pgee1 )
ptmat(2,1) = ptmat(1,2)
ptmat(2,2) = (pgee0-cee+pgee3)/( (cee-pgee3)**2 - pgee0*pgee0 + pgee1*pgee1 )

mtmat(1,1) = (mgee0+cee-mgee3)/( (cee-mgee3)**2 - mgee0*mgee0 + mgee1*mgee1 )
mtmat(1,2) = mgee1/(  (cee-mgee3)**2 - mgee0*mgee0 + mgee1*mgee1 )
mtmat(2,1) = mtmat(1,2)
mtmat(2,2) = (mgee0-cee+mgee3)/( (cee-mgee3)**2 - mgee0*mgee0 + mgee1*mgee1 )

!ptmat = dcmplx(0.0d0,0.01d0)  !*****************
!mtmat = dcmplx(0.0d0,0.01d0)


!---- begin calculating qpi-dos ---------
do q1 = -numk-1, numk-1    ! momentum transfer q = k+k'

     q2 = q1
     subdos(q1,w) = 0.0d0

do k1 = -numk-1,numk-1      ! momentum k prime
do k2 = -numk-1,numk-1
     
     kk1 = q1+k1
     kk2 = q2+k2
     kk3 = -q1+k1
     kk4 = -q2+k2     

     if ( abs(kk1) < numk-1 .and. abs(kk2) < numk-1 .and. &
        abs(kk3) < numk-1 .and.  abs(kk4) < numk-1   ) then ! if q inside BZ

     do i = 1,2   ! matrix indices
     do j = 1,2
       pgren1(i,j) = dcmplx(0.0d0,0.0d0)
       mgren1(i,j) = dcmplx(0.0d0,0.0d0)
       pgren2(i,j) =  dcmplx(0.0d0,0.0d0)
       mgren2(i,j) =  dcmplx(0.0d0,0.0d0)
        
     do l = 1,2    ! matrix multiplication
     do m = 1,2

       pgren1(i,j) = pgren1(i,j) + &
         pgren0(k1,k2,i,l)*ptmat(l,m)*pgren0(kk1,kk2,m,j)
      
       mgren1(i,j) = mgren1(i,j) + &
         mgren0(k1,k2,i,l)*mtmat(l,m)*mgren0(kk1,kk2,m,j)
       
       pgren2(i,j) = pgren2(i,j) + &
         pgren0(k1,k2,i,l)*ptmat(l,m)*pgren0(kk3,kk4,m,j)
       
       mgren2(i,j) = mgren2(i,j) + &
         mgren0(k1,k2,i,l)*mtmat(l,m)*mgren0(kk3,kk4,m,j)

     end do   ! end m
     end do   ! end l
     end do   ! end j
     end do   ! end i

  
     subdos(q1,w) = subdos(q1,w) + dimag( pgren1(1,1) + mgren1(2,2) &
        - dconjg(pgren2(1,1)) - dconjg(mgren2(2,2))      )  
   
     end if



     if ( kk1 < -numk-1 ) then  ! if q is outside BZ
        p1 = kk1+2*numk
     elseif ( kk1 >= numk-1 ) then
        p1 = kk1-2*numk
     else
        p1 = kk1
     end if

     if (kk2 < -numk-1 ) then
        p2 = kk2+2*numk
     elseif ( kk2 >= numk-1 ) then
        p2 = kk2 - 2*numk
     else
        p2 = kk2
     endif

     if ( kk3 < -numk-1 ) then  ! if q is outside BZ
        p3 = kk3+2*numk
     elseif ( kk3 >= numk-1 ) then
        p3 = kk3-2*numk
     else
        p3 = kk3
     end if

     if (kk4 < -numk-1 ) then
        p4 = kk4+2*numk
     elseif ( kk4 >= numk-1 ) then
        p4 = kk4 - 2*numk
     else
        p4 = kk4
     endif



!if ( abs(p2)>numk-1  .or.  abs(p1)>numk-1 &
!     .or.  abs(p3)>numk-1  .or.  abs(p4)>numk-1 )then !check if p1,p2 are in BZ
!print*,  p1, p2, p3, p4
!endif

     do i = 1,2   ! matrix indices
     do j = 1,2
       pgren1(i,j) = dcmplx(0.0d0,0.0d0)
       mgren1(i,j) = dcmplx(0.0d0,0.0d0)
       pgren2(i,j) =  dcmplx(0.0d0,0.0d0)
       mgren2(i,j) =  dcmplx(0.0d0,0.0d0)

     do l = 1,2    ! matrix multiplication
     do m = 1,2

       pgren1(i,j) = pgren1(i,j) + &
         pgren0(k1,k2,i,l)*ptmat(l,m)*pgren0(p1,p2,m,j)
       mgren1(i,j) = mgren1(i,j) + &
         mgren0(k1,k2,i,l)*mtmat(l,m)*mgren0(p1,p2,m,j)
       pgren2(i,j) = pgren2(i,j) + &
         pgren0(k1,k2,i,l)*ptmat(l,m)*pgren0(p3,p4,m,j)
       mgren2(i,j) = mgren2(i,j) + &
         mgren0(k1,k2,i,l)*mtmat(l,m)*mgren0(p3,p4,m,j)

     end do   ! end m
     end do   ! end l
     end do   ! end j
     end do   ! end i


     subdos(q1,w) = subdos(q1,w) + dimag( pgren1(1,1) + mgren1(2,2) &
        - dconjg(pgren2(1,1)) - dconjg(mgren2(2,2))      )


end do   ! end k2
end do   ! end k1

end do   ! end q1
end do   ! end w (omega)


deallocate( pgren0, mgren0   )
deallocate(  phi, kai, zbar, eps )

    
call mpi_barrier(mpi_comm_world,ierr)
call mpi_allreduce(subdos,dos,size(subdos), mpi_double_precision, mpi_sum,&
        mpi_comm_world, ierr )
!call mpi_allreduce(subdos2,dos2,size(subdos2), mpi_double_precision, mpi_sum,&
 !       mpi_comm_world, ierr )
call mpi_barrier(mpi_comm_world,ierr)


if (rank==0)then
print*, 'producing output file'
endif


if (rank==0) then
open(21,file='qdisp-consist-lam81x4.dat')
do q1 = -numk-1,numk-1
     kx = dfloat(q1)*dk
do w = -numw,numw 

     write(21,*)   kx,  dfloat(w),  abs(dos(q1,w))/totalk

end do
end do
close(21)
end if

!if(rank==0)then
!open(22,file='qdisp2.dat')
!do w=-numw,numw
!   write(22,*)   dfloat(w), abs( dos(0,w) )/totalk
!end do
!close(22)
!end if  


deallocate( dos, subdos    )



call mpi_barrier(mpi_comm_world,ierr)
call mpi_finalize(ierr)
stop

end program main











