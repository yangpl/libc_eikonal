!############################################################################
!Module for reading the inputs models. The models should be in binary format.
!the inputs should be write in following form in a file named by "inputs" :
!================================================ inputs 
!n1 n2                  ! Size of the model
!h1 h2                  ! Grid spscing
!max_iter               ! Number of iteration
!x1_src, x2_src         ! Position of the source (z,x)
!v_ver_file             ! name of the vertical(zero) velocity model
!angles_file            ! name of the angle map
!eps_map_file           ! value of Epsilone (Thomsen Parameter) 
!del_map_file           ! value of Delta (Thomsen Parameter)
!=================================================
!############################################################################
! v0s           : "name of the vertical velocity file" (m/sec)
! theta_maps          : "name of the rotation map (theta) map "
! epss         : "Value of the Epsilon"
! dels         : "Value of the Delta"
! x1_src,x2_src   : "source position in meter in direction 1 , 2 " 
! n1 ,n2          :  model dimension in grid point
! Max_iter        :The maximum number of iteration.
! h               :Grid step size (descretization step)
 
  subroutine read_inputs(n1,n2,xh1_s,xh2_s,max_iter,x1_src,x2_src,unitf)
    implicit none

    integer(kind=4),intent(out) :: n1,n2, max_iter
    real(kind=4),intent(out) :: x1_src,x2_src,xh1_s,xh2_s
    character(len=250)  :: config_file
    integer(kind=4),intent(inout)  :: unitf
    integer(kind=4) :: ierr
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    ! Read input file
    config_file = 'inputs'
    open(unit = unitf, file = config_file, status = 'old', action='read', iostat=ierr)
    if (ierr /= 0) then
      write(*,*) " error : could not open file : ", trim(adjustl(config_file))
      stop
    endif

    read (unitf,*)   n1, n2          ! "number of grid points on z and x: n1 n2"
    read (unitf,*)   xh1_s,xh2_s     ! "grid spacing "
    read (unitf,*)   max_iter        ! "number of maximum iterations "
    read (unitf,*)   x1_src, x2_src          ! "source position in meter in direction 1 , 2 "

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>

    write(*,*) "_____________________________________________"
    write(*,*)  "check point1: inputs reading done"
    write(*,*) " file is still open"
    write(*,*) "============================================="

    return
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>
  end subroutine read_inputs


  subroutine read_medium(n1,n2,v0s,unitf)


    real(kind=4),dimension(0:n1+1,0:n2+1)     :: v0s
    character(len=250)  :: v_ver_file
    integer(kind=4)                           :: unitf, ierr,i1,i2
    real(kind=4) :: vhomo
    integer(kind=4) :: ioption



    real(kind=4),dimension(:,:),allocatable            :: valeur

!=================================== continue the reading of file 'inputs'

    read (unitf,*)   v_ver_file      ! "name of the vertical velocity file"
    read (unitf,*)   ioption,vhomo
    close(unitf)

    write(*,*) "_____________________________________________"
    write(*,*)  "reading medium properties"
    write(*,*)  " file is now closed"
    write(*,*) "============================================="


    allocate(valeur(1:n1,1:n2))

    if(ioption == 1) then

       !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>
       ! Read subsurface model parameters files
       open(unit = unitf,file=v_ver_file,access='direct',recl=n1*n2*4, iostat=ierr)
       if (ierr /= 0) then
          write(*,*) " error : could not open file : ", trim(adjustl(v_ver_file))
          stop
       endif
       read(unitf,rec=1) valeur(:,:)          ! the vertical velocity for each point.
    else
       valeur(:,:)=vhomo
    endif

    do i2=1,n2
       do i1=1,n1
          v0s(i1,i2)=valeur(i1,i2)
       enddo
    enddo

    do i1=1,n1
       v0s(i1,0)=v0s(i1,1)
    enddo
    do i2=1,n2
       v0s(0,i2)=v0s(1,i2)
    enddo
    v0s(0,0)=v0s(1,1)
    v0s(n1+1,n2+1)=v0s(n1,n2)
    v0s(n1+1,0)=v0s(n1,1)
    v0s(0,n2+1)=v0s(1,n2)
    close(unitf)

    deallocate(valeur)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>

    write(*,*) "_____________________________________________"
    write(*,*)  "end of reading medium properties"
    write(*,*) "============================================="

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>
  end subroutine read_medium
