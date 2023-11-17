!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_transfer 
!***********************************************************************
!                                                                      !
!  Transfer data between analysis and filter grid                      !
!                                                                      !
! Modules: kinds, mg_parameter, mg_intstate, mg_bocos, mg_interpolate, !
!          mg_timers, mg_mppstuff                                      !
!                                                     M. Rancic (2021) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter
use mg_intstate, only: VALL,WORKA,WORK
use mg_intstate, only: cvf1,cvf2,cvf3,cvf4,lref
use mg_timers
use mg_mppstuff, only:  mype,ierror,mpi_comm_world
use mg_mppstuff, only: nx,my,mpi_comm_comp,finishMPI
use mg_interpolate, only: lwq_vertical_adjoint_spec
use mg_interpolate, only: lwq_vertical_direct_spec
use mg_interpolate, only: l_vertical_adjoint_spec
use mg_interpolate, only: l_vertical_direct_spec
use mg_interpolate, only: l_vertical_adjoint_spec2
use mg_interpolate, only: l_vertical_direct_spec2
!TEST
!use mg_output, only: output_spec1_2dd
!TEST
implicit none

integer(i_kind):: n,m,l,k,i,j

public anal_to_filt_all
public filt_to_anal_all

public anal_to_filt_all2
public filt_to_anal_all2

public anal_to_filt
public filt_to_anal

public stack_to_composite
public composite_to_stack

public S2C_ens
public C2S_ens


!TEST
!use mg_output, only: output_spec1_2dd
!TEST


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt_all
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
real(r_kind),allocatable,dimension(:,:,:,:):: A3D
real(r_kind),allocatable,dimension(:,:,:):: A2D

real(r_kind),allocatable,dimension(:,:,:,:):: F3D

!----------------------------------------------------------------------

allocate(A3D(km3_all,1:nm,1:mm,lm_a))
allocate(F3D(km3_all,1:nm,1:mm,lm))

                                                 call btim(  s2c_tim)

     call S2C_ens(WORKA,A3D,1,nm,1,mm,lm_a,km_a,km_a_all)

                                                 call etim(  s2c_tim)

                                                 call btim(  vadj_tim)
  if(lm_a>lm) then
    if(l_lin_vertical) then
       call l_vertical_adjoint_spec(km3_all,lm_a,lm,1,nm,1,mm,A3D,F3D)
    else
       call lwq_vertical_adjoint_spec(km3_all,lm_a,lm,1,nm,1,mm,               &
                                      cvf1,cvf2,cvf3,cvf4,lref,A3D,F3D)
    endif
  else 

    do L=1,lm
      F3D(:,:,:,L)=A3D(:,:,:,L)
    enddo

  endif
                                                 call etim(  vadj_tim)

                                                 call btim(  c2s_tim)

      call C2S_ens(F3D,WORK,1,nm,1,mm,lm,km,km_all)

                                                 call etim(  c2s_tim)

                                                 call btim(  a2f_tim)
     call anal_to_filt
                                                 call etim(  a2f_tim)


deallocate(A3D,F3D)
!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal_all
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
real(r_kind),allocatable,dimension(:,:,:,:):: A3D
real(r_kind),allocatable,dimension(:,:,:,:):: F3D

!----------------------------------------------------------------------

allocate(F3D(km3_all,1:nm,1:mm,lm))

allocate(A3D(km3_all,1:nm,1:mm,lm_a))

    call filt_to_anal
   
    call S2C_ens(WORK,F3D,1,nm,1,mm,lm,km,km_all)
  
 if(lm_a>lm) then
   if(l_lin_vertical) then
     call l_vertical_direct_spec(km3_all,lm,lm_a,1,nm,1,mm,F3D,A3D)
   else 
     call lwq_vertical_direct_spec(km3_all,lm,lm_a,1,nm,1,mm,              &
                                   cvf1,cvf2,cvf3,cvf4,lref,F3D,A3D)
   endif
 else

   do L=1,lm
     A3D(:,:,:,L)=F3D(:,:,:,L)
   enddo
 
 endif  

    call C2S_ens(A3D,WORKA,1,nm,1,mm,lm_a,km_a,km_a_all)

deallocate(A3D,F3D)
!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt_all2
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************


                                                 call btim(  vadj_tim)
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'before WORK=WORKA'
!     endif
!     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  if(lm_a>lm) then
     call l_vertical_adjoint_spec2(km3*n_ens,lm_a,lm,1,nm,1,mm,WORKA,WORK)
  else 
     WORK = WORKA
  endif
                                                 call etim(  vadj_tim)


                                                 call btim(  a2f_tim)
     call anal_to_filt
                                                 call etim(  a2f_tim)


!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal_all2
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************


    call filt_to_anal
   
  if(lm_a>lm) then
     call l_vertical_direct_spec2(km3*n_ens,lm,lm_a,1,nm,1,mm,WORK,WORKA)
  else 
     WORKA = WORK
  endif

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite                   &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(in):: ARR_ALL
real(r_kind),dimension(km3,1-hx:im+hx,1-hy:jm+hy,lm),intent(out):: A3D
real(r_kind),dimension(km2,1-hx:im+hx,1-hy:jm+hy)   ,intent(out):: A2D
!----------------------------------------------------------------------
    do L=1,lm
      do j=1-hy,jm+hy
      do i=1-hx,im+hx
        do k=1,km3
          A3D(k,i,j,L)=ARR_ALL( (k-1)*lm+L,i,j )
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          A2D(k,:,:)=ARR_ALL(km3*lm+k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack                   &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km2,1-hx:im+hx,1-hy:jm+hy),   intent(in):: A2D
real(r_kind),dimension(km3,1-hx:im+hx,1-hy:jm+hy,lm),intent(in):: A3D
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(out):: ARR_ALL
integer(i_kind):: i,j,L
!----------------------------------------------------------------------
    do L=1,lm
      do j=1-hy,jm+hy
      do i=1-hx,im+hx
        do k=1,km3
          ARR_ALL( (k-1)*lm+L,i,j )=A3D(k,i,j,L)
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          ARR_ALL(km3*lm+k,:,:)=A2D(k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine S2C_ens                              &
!***********************************************************************
!                                                                      !
! General transfer data from stack to composite variables for ensemble !
!                                                                      !
!***********************************************************************
(ARR_ALL,A3D,imn,imx,jmn,jmx,lmx,kmx,kmx_all)
!----------------------------------------------------------------------
use mg_parameter, only: km3,n_ens,km3_all
implicit none
integer, intent(in):: imn,imx,jmn,jmx,lmx,kmx,kmx_all
real(r_kind),dimension(kmx_all,imn:imx,jmn:jmx)    ,intent(in):: ARR_ALL
real(r_kind),dimension(km3_all,imn:imx,jmn:jmx,lmx),intent(out):: A3D
integer(i_kind):: k,L
integer(i_kind):: n,n_inc
!----------------------------------------------------------------------
  do n=1,n_ens
    n_inc = kmx*(n-1)

    do L=1,lmx
      do j=jmn,jmx
      do i=imn,imx 
        do k=1,km3
          A3D(km3*(n-1)+k,i,j,L)=ARR_ALL(n_inc+(k-1)*lmx+L,i,j)
        enddo
      enddo
      enddo
    enddo

   enddo
 
!----------------------------------------------------------------------
                        endsubroutine S2C_ens                         

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine C2S_ens                              &
!***********************************************************************
!                                                                      !
! General transfer data from composite to stack variables for ensemble !
!                                                                      !
!***********************************************************************
(A3D,ARR_ALL,imn,imx,jmn,jmx,lmx,kmx,kmx_all)
!----------------------------------------------------------------------
use mg_parameter, only: km2,km3,n_ens,km2_all,km3_all
implicit none
integer, intent(in):: imn,imx,jmn,jmx,lmx,kmx,kmx_all
real(r_kind),dimension(km3_all,imn:imx,jmn:jmx,lmx),intent(in):: A3D
real(r_kind),dimension(kmx_all,imn:imx,jmn:jmx)    ,intent(out):: ARR_ALL
integer(i_kind):: k,L
integer(i_kind):: n,n_inc
!----------------------------------------------------------------------

  do n=1,n_ens
    n_inc = kmx*(n-1)

     do L=1,lmx
       do j=jmn,jmx
       do i=imn,imx 
         do k=1,km3
           ARR_ALL(n_inc+(k-1)*lmx+L,i,j )= A3D(km3*(n-1)+k,i,j,L)
         enddo
       enddo
       enddo
     enddo

  enddo
!----------------------------------------------------------------------
                        endsubroutine C2S_ens                              

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_adjoint_offset
use mg_interpolate, only: quad_adjoint_offset
use mg_interpolate, only: lin_adjoint_offset
use mg_bocos, only:  bocoT_2d
implicit none

integer(i_kind):: ibm,jbm


!----------------------------------------------------------------------
       VALL=0.

!     ibm = ib
!     jbm = jb

                                                call btim(  aintp_tim)
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'before ..._adjoint_offset'  
!     endif
!     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       if(l_lin_horizontal) then
         ibm=1
         jbm=1
         call lin_adjoint_offset(WORK,VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),km_all,ibm,jbm)
       else &
       if(l_quad_horizontal) then
         ibm=2
         jbm=2
         call quad_adjoint_offset(WORK,VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),km_all,ibm,jbm)
       else
         ibm=3
         jbm=3
         call lsqr_adjoint_offset(WORK,VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),km_all,ibm,jbm)
       endif


                                                 call etim(  aintp_tim)


!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    

                                                 call btim(  bocoT_tim)

       if(.not.l_lin_horizontal) then
         call bocoT_2d(VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),km_all,im,jm,ibm,jbm)
       endif

                                                 call etim(  bocoT_tim)
 

!----------------------------------------------------------------------
                        endsubroutine anal_to_filt

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_direct_offset
use mg_interpolate, only: quad_direct_offset
use mg_interpolate, only: lin_direct_offset
use mg_bocos, only:  boco_2d
implicit none

integer(i_kind):: ibm,jbm

!----------------------------------------------------------------------
!     ibm = ib
!     jbm = jb
    
       if(l_lin_horizontal) then
         ibm=1
         jbm=1
       else &
       if(l_quad_horizontal) then
         ibm=2
         jbm=2
       else 
         ibm=3
         jbm=3
       endif

!***
!***  Supply boundary conditions for VALL
!***
                                                 call btim(  boco_tim)

       if(.not.l_lin_horizontal) then
         call boco_2d(VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),km_all,im,jm,ibm,jbm)
       endif

                                                 call etim(  boco_tim)


                                                 call btim(   intp_tim)
       if(l_lin_horizontal) then
         call lin_direct_offset(VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),WORK,km_all,ibm,jbm)
       else &
       if(l_quad_horizontal) then
         call quad_direct_offset(VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),WORK,km_all,ibm,jbm)
       else
         call lsqr_direct_offset(VALL(1:km_all,1-ibm:im+ibm,1-jbm:jm+jbm),WORK,km_all,ibm,jbm)
       endif

                                                 call etim(   intp_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_transfer
