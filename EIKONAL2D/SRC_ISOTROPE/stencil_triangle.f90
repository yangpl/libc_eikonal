!############################################################################
!# this subroutine for each grid point C provides the expected time time_C
!# for the triangle B-C-A
!
! Isotropic medium : 8-nodes FD stencil 
!
! SRC dedicated to training ...  July 2017 Jean Virieux
!
!############################################################################
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! norder = 1     elementary stencil ... for cartesian grid ... 8 angular sections
! three points   point A, point B where known travel times and positions
!                point C where the travel time should be computed
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SUBROUTINE triangle_solver(xa,za,xb,zb,xc,zc, &
     vpa, &
     vpb, &
     vpc, &
     Time_A,Time_B,Time_C)

  use Polynomial2RootSolvers, ONLY : quadraticRoots

  IMPLICIT NONE

  INTEGER, PARAMETER :: wp = KIND(0.0D0)
  !INTEGER, PARAMETER :: wp = KIND(0.0)

  real(kind=8),INTENT(in)          :: xa,za,xb,zb,xc,zc
  real(kind=8),INTENT(in)          :: vpa
  real(kind=8),INTENT(in)          :: vpb
  real(kind=8),INTENT(in)          :: vpc
  real(kind=8),INTENT(in)          :: Time_A,Time_B
  real(kind=8),INTENT(out)         :: Time_C

  real(kind=8)                     :: mat_a,mat_b,mat_c,mat_d,det_inv,tampon
  real(kind=8)                     :: mat_m,mat_n,mat_p,mat_q
  real(kind=8)                     :: tc,d_ac,d_bc,tau ! ,diff_tau,tamat_1,tamat_2
  real(kind=8)                     :: up,up2,dispersion
  real(kind=8)                     :: px_c,pz_c,pinv
  real(kind=8)                     :: xnx,xnx_c,xnz,xnz_c

  real(kind=8),parameter           :: eps=1.0d-12

!==================== use the wp for consistency with quadratic code of TOMS 954
  real(kind=wp)                    :: q1,q0
  real(kind=wp)                    :: racine(1:2,1:2)   ! racine(n,1) real part and racine(n,2) imag. part
  integer(kind=4)                  :: nReal

  real(kind=8)                     :: Time_C_OK,Time_AC,Time_BC
  integer(kind=4)                  :: n_d,iopt,i,norder

  norder=1                         ! set it to first-order approximation of derivatives
                                   ! not accurate enough for applications ...
  !====================================================================================
  ! A-C-B      
  ! A (xa,za) B(xb,zb)  C(xc,zc)
  ! one of elementary stencil
  !========================== linear expansion of (px,pz) wrt Tc
  !====================================================================================
  mat_a=xc-xa; mat_b=zc-za
  mat_c=xc-xb; mat_d=zc-zb
  !                                            matrice composing the slowness contribution
  !   (                            )
  !   ( (xc-xa)      (zc-za)       )
  !   ( ((xc-xb)     (zc-zb)       )
  !   (                            )
  !
  !======================= get the inverse of the matrix 2x2
  det_inv=1.0d0/(mat_a*mat_d-mat_b*mat_c)     !
  tampon=mat_a
  mat_a=mat_d*det_inv                         !   1     (d   -b)   (a'  b')
  mat_d=tampon*det_inv                        ! -----   (      ) = (       )
  mat_c=-mat_c*det_inv                        !(ad-bc)  (-c   a)   (c'  d')
  mat_b=-mat_b*det_inv                        !

  ! get expression of Tc     px=M Tc + N
  !                          pz=P Tc + Q
  !                      to be found
  !
  !   a'(Tc-Ta) + b'(Tc-Tb) = M Tc + N
  !   c'(Tc-Ta) + a'(Tc-Tb) = P Tc + Q
  !
  !   M = a'+b'     N = -a'*Ta-b'*Tb   (1st order)
  !   P = c'+d'     Q = -c'*Ta-d'*Tb
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ first-order stencil
  if(norder == 1) then
     mat_m=mat_a+mat_b
     mat_n=-mat_a*Time_A-mat_b*Time_B
     mat_p=mat_d+mat_c
     mat_q=-mat_c*Time_A-mat_d*Time_B
  else
     write(*,*) ' error in the ordre 1 or 2 while it is: ',norder
     stop 
  endif

  !========================= equation of second degree
  !
  !(M*M+P*P) Tc*Tc + 2 (MN+PQ) Tc +N*N+Q*Q-u*u=0  eikonal
  !
  ! transformed into
  !
  !  Tc*Tc + q1 Tc + q0 = 0
  !
  !=========================

  up=1.0d00/vpc           ! phase slowness
  up2=up*up              ! square of the phase slowness

  racine(:,:)=0.0D0

  if(abs(mat_m*mat_m + mat_p*mat_p) > eps) then
     !============ quadratic
     q1=2.0*(mat_m*mat_n + mat_p*mat_q)/(mat_m*mat_m + mat_p*mat_p)
     q0=(mat_n*mat_n + mat_q*mat_q-up2)/(mat_m*mat_m + mat_p*mat_p)
     call quadraticRoots (q1, q0, nReal, racine) 
     n_d=2
  else
     !============ linear because (mat_m*mat_m + mat_p*mat_p) = 0
     racine(1,1)=-0.5*(mat_n*mat_n + mat_q*mat_q-up2)/(mat_m*mat_n + mat_p*mat_q)
     racine(1,2)=0.0
     nReal=1
     n_d=1
     write(*,*) ' degenerate case'
  endif


  Time_C_OK=1.d0/eps   ! set to zero
  iopt=0
  do i=1,n_d           ! nReal = 0 when complex roots
     !========================== consider only real roots but not those equal to zero
     if(racine(i,2) < eps  .and. racine(i,1) > eps ) then

        !        write(*,*) ' ===================================== '
        Time_C=racine(i,1)   ! maybe the solution
        !========================== check the causality
        !========================== compute the ray velocity (group velocity) tangent
        !========================== to the ray (not normal to waveform which defines the phase velocity) 
        px_c=mat_m*Time_C+mat_n
        pz_c=mat_p*Time_C+mat_q

        pinv=1.0d0/dsqrt(px_c*px_c+pz_c*pz_c)

        xnx_c=px_c*pinv
        xnz_c=pz_c*pinv

        !        write(*,*) ' px,pz,p,1/p,xnx_c,xnz_c '
        !        write(*,*) px_c,pz_c,p,1.0/p,xnx_c,xnz_c

        dispersion=up*pinv - 1.0d0   ! dispersion p=slowness ... attention pinv est l'inverse de lenteur en fait

        if(dabs(dispersion) < 0.000001) then     ! OK we are on P phase slowness
           !           write(*,*) ' entering in the causality cone dispersion',dispersion
           !========================== check the causality
           !========================== compute the ray velocity 
           !========================== to the ray (not normal to waveform // to slowness p)
           !=================================
           !   the point defined by x=xc+tc*pxc should belong to the line (AB)
           !                        z=zc+tc*pzc
           !   equation of line (AB)   (x-xa)*(zb-za)-(z-za)*(xb-xa)=0
           !   we deduce     tc [(zb-za)*pxc-(xb-xa)*pzc]=(zc-za)*(xb-xa)-(xc-xa)*(zb-za)
           !=================================
           tc=((xb-xa)*zc-(zb-za)*xc+zb*xa-xb*za)/((zb-za)*px_c-(xb-xa)*pz_c)   ! parameter on the line (c,p_c)
           !                                if(abs(tc) > 1.d21) write(*,*) ' tc ',tc
           if(tc < 0.0d0) then   ! on pointe dans la bonne direction  (time in C 
              !                       should increase and, therefore, we have to look backward from C toward A and B 
              if(dabs(xa-xb) < eps) then
                 tau=(zc-za+tc*pz_c)/(zb-za)    ! (zb-za) should be different from zero
              else
                 tau=(xc-xa+tc*px_c)/(xb-xa)
              endif
              !======================================== inside the causality cone (between A and B)
              if(tau <= 1.0d00 .and. tau >=0.d00) then  
                 Time_C_OK=dmin1(Time_C,Time_C_OK)
                 iopt=1
                 !                write(*,*) ' ONE SOLUTION ===> ',Time_C," tc",tc
              endif
           endif
        endif ! end of the slowness surface P
     endif ! end over real solutions
  enddo ! end of the loop over the four solutions
  !=================================
  ! get the final value which should exist in any case
  !=================================
  if(iopt > 0) then
     Time_C=1.0d12
     do i=1,iopt
        Time_C=dmin1(Time_C,Time_C_OK)    ! take the causal smallest value
     enddo
     return
  else
     !==================================== should get the viscous solution (two options)
     !============= from A   @@@@@@@@@@@@@@@@@@@@@@@@@
     xnx=(xc-xa)
     xnz=(zc-za)
     d_ac=dsqrt(xnx*xnx+xnz*xnz)
     Time_AC=Time_A+d_ac*up   ! time as if coming from A with phase velocity vqp
     !============= from B   @@@@@@@@@@@@@@@@@@@@@@@@@
     xnx=(xc-xb)
     xnz=(zc-zb)
     d_bc=dsqrt(xnx*xnx+xnz*xnz)
     Time_BC=Time_B+d_bc*up   ! time as if coming from B with phase velocity vqp
     !============= set the minimum time so we always get one solution
     Time_C=dmin1(Time_AC,Time_BC)

  endif  ! end with an existing solution (causal or viscous)
  return
end subroutine triangle_solver




