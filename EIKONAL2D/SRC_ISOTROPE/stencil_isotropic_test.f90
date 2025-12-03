!***************************************************************************
!
! test for isotropic case
!
! modification for norder=1 only
!
!***************************************************************************
program stencil_isotropic_test

implicit none

real(kind=8) :: xa,za,xb,zb,xc,zc
real(kind=8) :: va,va1
real(kind=8) :: vb,vb1
real(kind=8) :: vc
real(kind=8) :: Time_A,Time_B,Time_C
real(kind=8) :: x0,z0
integer(kind=4) :: norder

write(*,*) ' isotropic case '
100 continue
write(*,*) ' enter norder '
read(*,*) norder

if (norder /=1) stop

write(*,*) ' enter xa,za '
read(*,*) xa,za
write(*,*) ' enter xb,zb '
read(*,*) xb,zb
write(*,*) ' enter xc,zc '
read(*,*) xc,zc

write(*,*) ' enter the source position '
read(*,*) x0,z0

write(*,*) ' va,va1 '
read(*,*) va,va1
write(*,*) ' vb,vb1 '
read(*,*) vb,vb1
write(*,*) ' vc '
read(*,*) vc

Time_A=dsqrt((xa-x0)*(xa-x0)+(za-z0)*(za-z0))/va
Time_B=dsqrt((xb-x0)*(xb-x0)+(zb-z0)*(zb-z0))/vb


call triangle_solver(xa,za,xb,zb,xc,zc, &
                   va,va1,vb,vb1,vc,  &
                Time_A,Time_B,Time_C)

write(*,*) ' Time_A,Time_B,Time_C '
write(*,*) Time_A,Time_B,Time_C

goto 100
stop
end program stencil_isotropic_test
