!#############################################################################
!# this subroutine define the some offsets (i_shf) which help us to to define the 
!# neighbors of each grid point, which have roles in elementary stencils.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ we have eight elementary triangles
!
!
!
!Finite difference stencil: 
!
!
!                               X                N               Y
!                                 x        1     !     8       x
!                                   x            !           x
!                                     x          !         x             ^
!                                       x        !       x               |
!                                         x      !     x                 |
!                                    2      x    !   x       7           |======>            
!                                             x  ! x      
!                               W----------------C-----------------E
!                                             x  ! x
!                                           x    !   x       6
!                                   3     x      !     x
!                                       x        !       x
!                                     x          !         x
!                                   x      4     !    5      x
!                                 x              !             x
!                               T                S               Z
!
!
!
!
!
!                    
!C  ---->> C(i1_C,i2_C)
!A  ---->> A(i1_C+i1_shftA,i2_c+i2_shftA)
!B  ---->> B(i1_C+i1_shftB,i2_c+i2_shftB)
!
!#############################################################################
SUBROUTINE sub_neighbors(neig_ind)
  IMPLICIT NONE
  INTEGER(kind=4),DIMENSION(8,8),intent(out)   ::neig_ind

  INTEGER(kind=4)                         ::trg
!  INTEGER(kind=4)                         ::i1_shfA,i2_shfA,i1_shfB,i2_shfB

  do trg=1,8
  !_______________________Triangle 1  NCX   A=N and B=X
  IF(trg==1) THEN
!     i1_shfA=0
!     i2_shfA=1
!     i1_shfB=-1
!     i2_shfB=1
!neig_ind(1,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(1,1:4)= (/  0, 1,-1, 1 /)
neig_ind(1,5:8)= (/  0, 2,-2, 2 /)
     !_______________________Triangle XCW   A=X and B=W
  ELSEIF(trg==2) THEN 
!     i1_shfA=-1
!     i2_shfA=1
!     i1_shfB=-1
!     i2_shfB=0
!neig_ind(2,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(2,1:4)= (/ -1, 1,-1, 0 /)
neig_ind(2,5:8)= (/ -2, 2,-2, 0 /)
     !_______________________Triangle WCT   A=W and B=T
  ELSEIF(trg==3) THEN  
!     i1_shfA=-1
!     i2_shfA=0
!     i1_shfB=-1
!     i2_shfB=-1
!neig_ind(3,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(3,1:4)= (/ -1, 0,-1,-1 /)
neig_ind(3,5:8)= (/ -2, 0,-2,-2 /)
     !_______________________Triangle TCS   A=T and B=S
  ELSEIF(trg==4) THEN
!     i1_shfA=-1
!     i2_shfA=-1
!     i1_shfB=0
!     i2_shfB=-1 
!neig_ind(4,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(4,1:4)= (/ -1,-1, 0,-1 /)
neig_ind(4,5:8)= (/ -2,-2, 0,-2 /)
     !_______________________Triangle SCZ   A=S and B=Z
  ELSEIF(trg==5) THEN
!     i1_shfA=0
!     i2_shfA=-1
!     i1_shfB=1
!     i2_shfB=-1 
!neig_ind(5,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(5,1:4)= (/  0,-1, 1,-1 /)
neig_ind(5,5:8)= (/  0,-2, 2,-2 /)
     !_______________________Triangle ZCE   A=Z and B=E
  ELSEIF(trg==6) THEN
!     i1_shfA=1
!     i2_shfA=-1
!     i1_shfB=1
!     i2_shfB=0 
!neig_ind(6,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(6,1:4)= (/  1,-1, 1, 0 /)
neig_ind(6,5:8)= (/  2,-2, 2, 0 /)
     !_______________________Triangle ECY   A=E and B=Y
  ELSEIF(trg==7) THEN
!     i1_shfA=1
!     i2_shfA=0
!     i1_shfB=1
!     i2_shfB=1 
!neig_ind(7,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(7,1:4)= (/  1, 0, 1, 1 /)
neig_ind(7,5:8)= (/  2, 0, 2, 2 /)
     !_______________________Triangle YCN   A=Y and B=N
  ELSEIF(trg==8) THEN
!     i1_shfA=1
!     i2_shfA=1
!     i1_shfB=0
!     i2_shfB=1 
!neig_ind(8,1:4)= (/ i1_shfA,i2_shfA,i1_shfB,i2_shfB /)
neig_ind(8,1:4)= (/  1, 1, 0, 1 /)
neig_ind(8,5:8)= (/  2, 2, 0, 2 /)
     !______________________

  ENDIF



  end do

END SUBROUTINE sub_neighbors
