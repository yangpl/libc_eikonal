!############################################################################
!# The main subroutine in fast sweeping loop. for each grid point we need to 
!# acquire eight travel times for eight triangles. 
!# Each triangle plays the role of finite-difference stencil.
!############################################################################
! defining variables:
! n1,n2,h1,h2 : grid size in direction 1,2 and the grid step size.
! TT : traveltime map.
! neig_ind : the matrix which contains the position of each neighbor.
! trg : the index of triangles we use in finite difference.
!
! Compute the least time at the point (i1,i2) return in the variable output
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE local_solver(neig_ind,i1,i2,n1,n2,xh1_s,xh2_s,TT_s,vel_p_s,output)
  IMPLICIT NONE
  REAL(kind=4),INTENT(in)                          :: xh1_s,xh2_s
  INTEGER(kind=4),INTENT(IN)                       :: i1,i2,n1,n2
  REAL(kind=4),DIMENSION(0:n1+1,0:n2+1),INTENT(in) :: TT_s
  real(kind=4),INTENT(out)                         :: output
  real(kind=4),DIMENSION(0:n1+1,0:n2+1),INTENT(in) :: vel_p_s

  INTEGER(kind=4)                             :: trg
  INTEGER(kind=4),DIMENSION(8,8),INTENT(IN)   :: neig_ind
  REAL(kind=8)                                :: Far_time
  real(kind=8)                                :: d_ac,d_bc
  real(kind=8)                                :: vpa,vpb,vpc
  real(kind=8)                                :: Time_A,Time_B,Time_C
  real(kind=8)                                :: xa,za,xb,zb,xc,zc
  real(kind=8)                                :: xnx,xnz
  real(kind=8),dimension(8)                   :: temps

  Far_time=1.0d21    ! high values for testing

  !============================ current point i1,i2

  output=1.0d23

  !write(*,*) ' entering in local solver ',i1,i2,n1,n2

  do trg=1,8 ! loop for covering all the elementary triangular stencils. 

     TIME_C=1.0d21

     !write(*,*) ' trg ',trg,i1+neig_ind(trg,1),i2+neig_ind(trg,2),i1+neig_ind(trg,3),i2+neig_ind(trg,4)

     Time_A=dble(TT_s(i1+neig_ind(trg,1),i2+neig_ind(trg,2)))
     Time_B=dble(TT_s(i1+neig_ind(trg,3),i2+neig_ind(trg,4)))

     !write(*,*) 'Time A & B ',Time_A,Time_B

     ! a decision here: do we need to compute the new value ... 
     if(Time_A < Far_time .and. Time_B < Far_time) then ! we have previously updated values

        vpa=vel_p_s(i1+neig_ind(trg,1),i2+neig_ind(trg,2))   ! get medium property at point A
        xa=dble(xh1_s)*dfloat(neig_ind(trg,1))                     ! and the geometry
        za=dble(xh2_s)*dfloat(neig_ind(trg,2))

        vpb=vel_p_s(i1+neig_ind(trg,3),i2+neig_ind(trg,4))   ! get medium property at point B
        xb=dble(xh1_s)*dfloat(neig_ind(trg,3))                     ! and the geometry
        zb=dble(xh2_s)*dfloat(neig_ind(trg,4))

        vpc=vel_p_s(i1,i2)                                   ! get medium property at point C
        xc=0.0d0                      ! and the geometry (center of the quadrangle stencil)
        zc=0.0d0
        !================================ call for the stencil  first order for the moment
        call triangle_solver(xa,za,xb,zb,xc,xc, &
             vpa, &
             vpb, &
             vpc, &
             Time_A,Time_B,Time_C)

     elseif(Time_A > Far_time .and. Time_B < Far_time) then ! only from B (time in B set)

        !============= from B
        vpb=vel_p_s(i1+neig_ind(trg,3),i2+neig_ind(trg,4))  ! define the parameters at point B
        xb=dble(xh1_s)*dfloat(neig_ind(trg,3))                     ! define the geometry
        zb=dble(xh2_s)*dfloat(neig_ind(trg,4))
        xc=0.0d0
        zc=0.0d0
        xnx=(xc-xb)
        xnz=(zc-zb)
        d_bc=dsqrt(xnx*xnx+xnz*xnz)
        Time_C=Time_B+d_bc/vpb   ! time as if coming from B with phase velocity vqp

     elseif(Time_A < Far_time .and. Time_B > Far_time) then !  only from A (time in A set)

        !============= from A
        vpa=vel_p_s(i1+neig_ind(trg,1),i2+neig_ind(trg,2))   ! define the parameters at point A
        xa=dble(xh1_s)*dfloat(neig_ind(trg,1))               ! define the geometry
        za=dble(xh2_s)*dfloat(neig_ind(trg,2))
        xc=0.0d0
        zc=0.0d0
        xnx=(xc-xa)
        xnz=(zc-za)
        d_ac=dsqrt(xnx*xnx+xnz*xnz)
        Time_C=Time_A+d_ac/vpa   ! time as if coming from A with phase velocity vqp

     endif ! end of checking uninitialized values

     temps(trg)=Time_C    ! we store these estimated times

     !if(iabs(i1-i1_src) < 2 .and. iabs(i2-i2_src) < 2) then
     !            write(*,*) Time_A,Time_B,neig_ind(trg,1),neig_ind(trg,2),neig_ind(trg,3),neig_ind(trg,4)
     !            write(*,*) i1_src,i2_src,i1,i2,' temps(i) ',trg,temps(trg)
     !endif

  enddo ! loop over triangles  trg

  output=DMIN1(temps(1),temps(2),temps(3),temps(4),temps(5),temps(6),temps(7),temps(8))

  return
END SUBROUTINE local_solver
