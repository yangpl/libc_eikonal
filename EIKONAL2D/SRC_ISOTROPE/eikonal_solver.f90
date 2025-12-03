SUBROUTINE eikonal_solver(n1,n2,Max_iter,x1_src_s,x2_src_s,xh1_s,xh2_s,vel_p_s,TT_s)
  !=========================================================================
  !n1,n2 : model size in grid point number in direction 1 and 2.
  !x2_src_s,x1_src_s : position of the source in meter in directions 2,1. """NO NEGATIVE NUMBER"""
  !xh1_s,xh2_s : grid step in direction 1 and in direction 2.
  !
  !vel_ps_s,del_p_s,eps_ps_s,tet_p_s : medium properties
  !
  !start_time,end_time : start and end time for code.
  !trg : index of triangles in finite difference scheme.
  !Taprx  : solution of local solver for each grid point.
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ! NOTES : Isotropic case ...
  !========================================================================
  IMPLICIT NONE   

  real(kind=4),DIMENSION(0:n1+1,0:n2+1),INTENT(in)     :: vel_p_s
  real(kind=4),DIMENSION(0:n1+1,0:n2+1),intent(inout)  :: TT_s   !  extended boundaries but solution between 1:n1;1:n2


  real(kind=8)                                     :: T_bound1,T_bound2,T_bound3,T_bound4
  real(kind=4),DIMENSION(0:n1+1,0:n2+1)            :: TT_old_s
  real(kind=4)                                     :: Taprx_s

  !=================================================================== sweep iteration
  integer(kind=4)                            :: Max_iter

  !=================================================================== source

  real(kind=4),intent(in)                    :: x1_src_s,x2_src_s
  real(kind=8)                               :: x1_src,x2_src
  !=================================================================== mesh
  integer(kind=4),intent(in)                 :: n1,n2
  real(kind=4),intent(in)                    :: xh1_s,xh2_s
  real(kind=8)                               :: x1_node,x2_node,xn1,xn2
  real(kind=8)                               :: vp
  real(kind=8)                               :: d
  integer(kind=4),DIMENSION(8,8)             :: neig_ind

  real(kind=8)                               :: Max_fsm


  real(kind=8),PARAMETER                     :: PI=3.141592653589793D0
  real(kind=8),PARAMETER                     :: Far_time=1.0d21        ! far value
  real(kind=8),PARAMETER                     :: err=1.0d-12            ! smallest real number
  real(kind=8),PARAMETER                     :: fsm_stp=1.0d-14         ! convergence of sweep 

  real(kind=8)                               :: start_time,end_time

  integer(kind=4)                            :: it_fsm,i1,i2,i1_src,i2_src
  integer(kind=4)                            :: sflag


  integer(kind=4) :: i11,i22,i1step,i2step,i1shift,i2shift

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> mapping for neighboring nodes of the stencil (first order)
  CALL sub_neighbors(neig_ind)                  ! eight triangles and shifts of two vectors (four values)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>Converting inputs to the Double precision.

  x2_src = DBLE(x2_src_s)
  x1_src = DBLE(x1_src_s)

  !  x1_src=FLOOR(x1_src/dble(xh1_s))*xh1_s                 ! we start the origin at zero
  !  x2_src=FLOOR(x2_src/dble(xh2_s))*xh2_s

  write(*,*) ' source ',x1_src,x2_src


  DO i2=0,n2+1    
     DO i1=0,n1+1                             
        TT_s(i1,i2)=Far_time                              ! set far value everywhere in the time table
     END DO
  END DO

  CALL CPU_TIME(start_time) 

  !~~~~~>>>>>>>here, we define the velocity of wave if source is not located on grid point.
  ! In the cases that the source is not located on grid, the boundary condition ( zero traveltime
  ! at the source position) will be replaced with 4 analytical value for surrounded grid points.

  IF (DABS(MOD(x2_src,dble(xh1_s))) < err .AND. DABS(MOD(x1_src,dble(xh2_s))) < err) THEN 
     sflag=1                                            ! we are on a grid node with the precision err
     i1_src=FLOOR(x1_src/dble(xh1_s))+1                 ! we start the origin at zero
     i2_src=FLOOR(x2_src/dble(xh2_s))+1

     TT_s(i1_src,i2_src)=0.0d0                    ! set the value of time at zero at the source node

    !============================================ first point (i1_src,i2_src)
i2step=0                                         ! step 0,2,4,6 ...
i1step=0
do i22=0,i2step,1
do i11=0,i1step,1

     i2shift=-i2step/2+i22-1                     !   si step=0 une valeur de -1
     i1shift=-i1step/2+i11-1                     !   si step=2 trois valeurs (-2,-1,0)
                                                 !   si step=4 cinq valeurs  (-3,-2,-1,0,2)
     x1_node=dble(xh1_s)*(i1_src-i1shift)
     x2_node=dble(xh2_s)*(i2_src-i2shift)
     xn1=x1_node-x1_src
     xn2=x2_node-x2_src
     d=DSQRT(xn1*xn1+xn2*xn2)                           ! distance source-n1-n1  first point of the rectangle
     vp=DBLE(vel_p_s(i1_src,i2_src))                    ! to be checked for indices  we start the origin at (1,1)
     TT_s(i1_src+1-i1shift,i2_src+1-i2shift)=d/vp
!write(*,*)  'i,j source, i,j sol, t',i1_src,i2_src,i1_src+1-i1shift,i2_src+1-i2shift,TT_s(i1_src+1-i1shift,i2_src+1-i2shift),vp
enddo
enddo

     WRITE(*,*) "Source is located on a node : (i1_src,i2_src) with spread", i1_src,i2_src,i1step,i2step

  ELSE 
     sflag=2

     i1_src=FLOOR(x1_src/dble(xh1_s))+1                  ! we start the origin at zero
     i2_src=FLOOR(x2_src/dble(xh2_s))+1

     !============================================ first point (i1_src,i2_src)
     x1_node=dble(xh1_s)*(i1_src-1)
     x2_node=dble(xh2_s)*(i2_src-1)
     xn1=x1_node-x1_src
     xn2=x2_node-x2_src
     d=DSQRT(xn1*xn1+xn2*xn2)    ! distance source-n1-n1  first point of the rectangle
     vp=DBLE(vel_p_s(i1_src,i2_src))                    ! to be checked for indices  we start the origin at (1,1)
     T_bound1=d/vp
     TT_s(i1_src,i2_src)=T_bound1

     !============================================ second point (i1_src,i2_src+1)
     x1_node=dble(xh1_s)*(i1_src-1)
     x2_node=dble(xh2_s)*(i2_src)
     xn1=x1_node-x1_src
     xn2=x2_node-x2_src
     d=DSQRT(xn1*xn1+xn2*xn2)    ! distance source-second point of the rectangle
     vp=DBLE(vel_p_s(i1_src,i2_src))                    ! to be checked for indices  we start the origin at (1,1)
     T_bound2=d/vp
     TT_s(i1_src,i2_src+1)=T_bound2

     !============================================ third point (i1_src+1,i2_src+1)
     x1_node=dble(xh1_s)*(i1_src)
     x2_node=dble(xh2_s)*(i2_src)
     xn1=x1_node-x1_src
     xn2=x2_node-x2_src
     d=DSQRT(xn1*xn1+xn2*xn2)    ! distance source-third point of the rectangle
     vp=DBLE(vel_p_s(i1_src,i2_src))                    ! to be checked for indices  we start the origin at (1,1)
     T_bound3=d/vp
     TT_s(i1_src+1,i2_src+1)=T_bound3

     !============================================ fourth point (i1_src+1,i2_src)
     x1_node=dble(xh1_s)*(i1_src)
     x2_node=dble(xh2_s)*(i2_src-1)
     xn1=x1_node-x1_src
     xn2=x2_node-x2_src
     d=DSQRT(xn1*xn1+xn2*xn2)    ! distance source-fourth point of the rectangle
     vp=DBLE(vel_p_s(i1_src,i2_src))                    ! to be checked for indices  we start the origin at (1,1)
     T_bound4=d/vp
     TT_s(i1_src+1,i2_src)=T_bound4


     WRITE(*,*) '^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^'
     WRITE(*,*) "Source is not located on a grid point",x1_src,x2_src,i1_src,i2_src
     WRITE(*,*) "The velocities around the source are as follows:"
     WRITE(*,*) vp                         ! isotrope here ...
     WRITE(*,*) "Times set at corners",T_bound1,T_bound2,T_bound3,T_bound4
     WRITE(*,*) '______________________________________________________________________________________'
  ENDIF
  !______

  !~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>the loop of "Fast Sweeping Method (FSM)" 

  Max_fsm=Far_time
  it_fsm=0


  DO WHILE (Max_fsm >= fsm_stp .and. it_fsm < Max_iter)   ! always two criteria for exit - convergence of maximum nber of loops
     !  DO WHILE (it_fsm < Max_iter)   ! always two criteria for exit - convergence of maximum nber of loops
     it_fsm=it_fsm + 1
     write(*,*) 'iter ',it_fsm,Max_iter,Max_fsm
     TT_old_s=TT_s

     Max_fsm=0.0D0

     !~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> Sweep 1
     DO i2=1,n2
        DO i1=1,n1
           CALL local_solver(neig_ind,i1,i2,n1,n2,xh1_s,xh2_s,TT_s,vel_p_s,Taprx_s)
           TT_s(i1,i2)=DMIN1(TT_s(i1,i2),Taprx_s)
        ENDDO
     END DO
    !~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> Sweep 2 
     DO i2=1,n2
        DO i1=n1,1,-1
           CALL local_solver(neig_ind,i1,i2,n1,n2,xh1_s,xh2_s,TT_s,vel_p_s,Taprx_s)
           TT_s(i1,i2)=DMIN1(TT_s(i1,i2),Taprx_s)
        END DO
     END DO
     !~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> Sweep 3
     DO i2=n2,1,-1
        DO i1=n1,1,-1
           CALL local_solver(neig_ind,i1,i2,n1,n2,xh1_s,xh2_s,TT_s,vel_p_s,Taprx_s)
           TT_s(i1,i2)=DMIN1(TT_s(i1,i2),Taprx_s)
        END DO
     END DO
     !~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> Sweep 4
     DO i2=n2,1,-1
        DO i1=1,n1
           CALL local_solver(neig_ind,i1,i2,n1,n2,xh1_s,xh2_s,TT_s,vel_p_s,Taprx_s)
           TT_s(i1,i2)=DMIN1(TT_s(i1,i2),Taprx_s)
 !======================== #####################################################    convergence ... 
              Max_fsm = DMAX1(Max_fsm, dabs(dble(TT_old_s(i1,i2)-TT_s(i1,i2))))
              !======================== ############################" to be included in the last sweep for testing convergence
        END DO
      ENDDO
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>>>>> End of FSM loop
  END DO
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  CALL CPU_TIME(end_time) 

  WRITE(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
  WRITE(*,*) "Convergence value max(TT_old,TT)",Max_fsm
  WRITE(*,*) "Number of iterations            ",it_fsm 
  WRITE(*,*) "++++++++++++++++++++++++++++++++++++++"
  WRITE(*,*) "The CPU run time in sec",end_time-start_time
  WRITE(*,*) "Stopping criterion for FSM", fsm_stp
  WRITE(*,*) "Maximum of iterations for FSM",Max_iter 
  WRITE(*,*) "Model size in grid points n1 n2",n1,n2
  WRITE(*,*) "++++++++++++++++++++++++++++++++++++++"
  return
END SUBROUTINE eikonal_solver
