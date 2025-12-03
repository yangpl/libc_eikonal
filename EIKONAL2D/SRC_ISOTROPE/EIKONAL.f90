!##########################################################################
!# This code calculates the first arrival travel times for Isotropic      #
!#   models by using a fast sweeping technique                            #
!#   one use non-linear upwind difference and Gauss-Seidel iterations     #
!#   with alternating sweeping ordering                                   #
!#                                                                        #
!#                                                                        #
!#                                                                        #
!##########################################################################
! n1,n2 : size of the model by grid point number
! i1,i2 : indexes in direction 1 and 2.
! Max_iter : Maximum iteration for Fixed point iteration as an threshold.
! v0s,theta_maps,epss,dels : Subsurface models <velocity,tilted angles,epsilon and delta>
! TT_map : computed traveltime map.
! TT.bin : Output file of traveltime map.
! x1_src_s,x2_src_s : position of the source in meter, direction 1, 2.
! xh1_s,xh2_s  : grid step size. 
! NOTE : The subscript "s" for some of above variables means they are in single precision.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> 
Program main

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>> 
  IMPLICIT NONE 
  INTEGER(kind=4)                                :: n1,n2,Max_iter,i1,i2
  REAL(kind=4)                                   :: x2_src_s,x1_src_s,xh1_s,xh2_s
  REAL(kind=4), DIMENSION(:,:),ALLOCATABLE       :: v0s,TT_map

  integer(kind=4)                                :: unitf=11

  call read_inputs(n1,n2,xh1_s,xh2_s,Max_iter,x1_src_s,x2_src_s,unitf)

  allocate(TT_map(0:n1+1,0:n2+1))
  allocate(v0s(0:n1+1,0:n2+1))

  call read_medium(n1,n2,v0s,unitf)

  call eikonal_solver(n1,n2,Max_iter,x1_src_s,x2_src_s,xh1_s,xh2_s,v0s,TT_map)

  OPEN(1,file='TT.bin',access='direct',recl=n1*n2*4)
!  WRITE(1,rec=1) ((TT_map(i1,i2),i2=n2,1,-1),i1=1,n1)        ! writing the final map of Travel-Time in "TT_test.bin" file
  WRITE(1,rec=1) ((TT_map(i1,i2),i1=1,n1),i2=1,n2)        ! writing the final map of Travel-Time in "TT_test.bin" file

  deallocate(v0s,TT_map)

  END PROGRAM main

  
  
  

 
