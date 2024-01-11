!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_filtering
!***********************************************************************
!                                                                      !
! Contains all multigrid filtering prodecures                          ! 
!                                                                      ! 
!                                                     M. Rancic (2020) !
!  Updates:                                           M. Rancic (2023) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter, only: im,jm,hx,hy,hz,lm,gm,Fimax,Fjmax
use mg_parameter, only: mgbf_line,lquart,lhelm
use mg_parameter, only: km2,km3,km
use mg_parameter, only: l_vertical_filter
!use mpimod, only: mype,ierror
use mg_mppstuff, only: mype,ierror
use mg_mppstuff, only: l_hgen,my_hgen,finishMPI,barrierMPI
use mg_generations, only: upsending_all,downsending_all,weighting_all
use mg_generations, only: upsending_ens,downsending_ens,weighting_ens
use mg_generations, only: upsending2_ens,downsending2_ens
use mg_transfer, only: stack_to_composite,composite_to_stack
use mg_transfer, only: S2C_ens,C2S_ens
use mg_bocos, only: boco_2d,bocoT_2d
use mg_bocos, only: boco_3d, bocoT_3d
use mg_bocos, only: bocox,bocoy
use mg_bocos, only: bocoTx,bocoTy
use jp_pbfil, only: rbeta,rbetaT
use jp_pbfil3, only: dibetat,dibeta
use mg_output

use mg_timers

public filtering_procedure 

private filtering_rad2_bckg
private filtering_lin2_bckg
private filtering_fast_bckg

private filtering_rad3
private filtering_lin3

private filtering_rad2_ens
private filtering_lin2_ens
private filtering_fast_ens

! For new development
private filtering_rad2_z_loc_g3
private filtering_rad2_z_loc_g4

private sup_vrbeta1
private sup_vrbeta1T
private sup_vrbeta3
private sup_vrbeta3T

private sup_vrbeta1_ens
private sup_vrbeta1T_ens

private sup_vrbeta1_bck
private sup_vrbeta1T_bck

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_procedure(mg_filt) 
!***********************************************************************
!                                                                      !
! Driver for Multigrid filtering procedures with Helmholtz operator    !
!                                                                      !
!   1, 2, 3: Radial filter                                             !
!         1: 2d radial filter for all variables                        !
! ->      2: 2d radial filter with 1d in vertical for 3d variables     !
!         3: 3d radial filter for 3d variables                         !
!                                                                      !
!   4, 5, 6: Line filter                                               !
!         4: 2d line filter for all variables                          !
!         5: 2d line filter with 1d in vertical for 3d variables       !
!         6: 3d line filter for 3d variables                           !
!                                                                      !
!                                                                      !
!***********************************************************************
implicit none 

integer(i_kind),intent(in):: mg_filt
!-----------------------------------------------------------------------
!  if(mgbf_line) then
!    if(mg_filt<4) then
!       print*,'("Line filters have options 4-6")'
!       stop
!    endif
!  else
!    if(mg_filt>3.and.mg_filt<7) then
!       print*,'("Radial filters have options 1-3,7,8")'
!       stop
!    endif
!  endif 

      select case(mg_filt)
        case(1)
          call filtering_rad2_bckg
        case(2)
          call filtering_lin2_bckg   
        case(3)
          call filtering_fast_bckg
        case(4)
          call filtering_rad3
        case(5)
          call filtering_lin3
        case(6)
          call filtering_rad2_ens
        case(7)
          call filtering_lin2_ens
        case(8)
          call filtering_fast_ens
       end select

!-----------------------------------------------------------------------
                        endsubroutine filtering_procedure    

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad3
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure                                        !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d radial filter                                               !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp2,pasp3,ss2,ss3
use mg_intstate, only: VALL,HALL
implicit none


real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D


integer(i_kind) L,i,j

!----------------------------------------------------------------------
allocate(VM3D(km3,1-hx:im+hx,1-hy:jm+hy,lm))                 ; VM3D=0.
allocate(VM2D(km2,1-hx:im+hx,1-hy:jm+hy   ))                 ; VM2D=0.
allocate(HM3D(km3,1-hx:im+hx,1-hy:jm+hy,lm))                 ; HM3D=0.
allocate(HM2D(km2,1-hx:im+hx,1-hy:jm+hy   ))                 ; HM2D=0.

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)


!***
!*** Apply adjoint of Beta filter at all generations 
!***



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call btim(   hfiltT_tim)
!
! Adjoint filtering
!
      call stack_to_composite(VALL,VM2D,VM3D)
        call rbetaT(km2,hx,1,im,hy,1,jm,pasp2,ss2,VM2D)
        call sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)

    if(l_hgen) then
      call stack_to_composite(HALL,HM2D,HM3D)
        call rbetaT(km2,hx,1,im,hy,1,jm,pasp2,ss2,HM2D)
        call sup_vrbeta3T(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL)
    endif 
                                                 call etim(   hfiltT_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call btim(    bocoT_tim)

        call bocoT_2d(VALL,km,im,jm,hx,hy)
        call bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(    bocoT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                 call btim( weight_tim)

      call weighting_all(VALL,HALL,lhelm)

                                                 call etim( weight_tim)


!***
!*** Apply Beta filter at all generations 
!***

                                                 call btim( boco_tim)

      call boco_2d(VALL,km,im,jm,hx,hy)
      call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(    boco_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call btim(   hfilt_tim)
!
! Filtering
!
      call stack_to_composite(VALL,VM2D,VM3D)
        call rbeta(km2,hx,1,im,hy,1,jm,pasp2,ss2,VM2D(:,:,:))
        call sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)
  if(l_hgen)  then
      call stack_to_composite(HALL,HM2D,HM3D)
        call rbeta(km2,hx,1,im,hy,1,jm,pasp2,ss2,HM2D(:,:,:))
        call sup_vrbeta3(km3,hx,hy,hz,im,jm,lm,pasp3,ss3,HM3D)
      call composite_to_stack(HM2D,HM3D,HALL)
  endif
                                                 call etim(   hfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!***
!***  Downsend, interpolate and add 
!***  Then zero high generations 
!***

                                                 call btim(   dnsend_tim)

       call downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)
deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)


!-----------------------------------------------------------------------
                        endsubroutine filtering_rad3   

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_lin3
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 6:                                     !
!                                                                      !
!     - Multiple of 2D  line filter                                    !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 3d line filter                      
!                                                                      !
!***********************************************************************
!TEST
use, intrinsic :: ieee_arithmetic
!TEST
use mg_parameter, only: nfil
use mg_intstate, only: dixs,diys,dizs,hss2,vpasp3
use mg_intstate, only: qcols,dixs3,diys3,dizs3
use mg_intstate, only: VALL,HALL
use jp_pkind2, only: fpi
implicit none

integer(i_kind) k,i,j,L
integer(i_kind) icol,iout,jout,lout
logical:: ff

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

real(r_kind), allocatable, dimension(:,:,:,:):: W
real(r_kind), allocatable, dimension(:,:,:,:):: H

integer(fpi), allocatable, dimension(:,:,:):: JCOL


allocate(VM3D(km3,1-hx:im+hx,1-hy:jm+hy,lm))                 ; VM3D=0.
allocate(VM2D(km2,1-hx:im+hx,1-hy:jm+hy   ))                 ; VM2D=0.
allocate(HM3D(km3,1-hx:im+hx,1-hy:jm+hy,lm))                 ; HM3D=0.
allocate(HM2D(km2,1-hx:im+hx,1-hy:jm+hy   ))                 ; HM2D=0.

allocate(W(km3,1-hx:im+hx,1-hy:jm+hy,1-hz:lm+hz))            ; W=0.
allocate(H(km3,1-hx:im+hx,1-hy:jm+hy,1-hz:lm+hz))            ; H=0.

allocate(JCOL(1:im,1:jm,1:Lm))                               ; JCOL=0

!-----------------------------------------------------------------------


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend
!***
     
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***

!
! From single stack to composite variables
!

       call stack_to_composite(VALL,VM2D,VM3D)
     if(l_hgen)  then
       call stack_to_composite(HALL,HM2D,HM3D)
     endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply adjoint filter to 2D variables first
!

       do icol=3,1,-1
                                                 call btim(   hfiltT_tim)
         call dibetat(km2,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
                                                 call etim(   hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_2d(VM2D,km2,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)
       enddo
    
     do icol=3,1,-1
                                                 call btim(   hfiltT_tim)
       if(l_hgen)  then
         call dibetat(km2,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
                                                 call etim(   hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)
     enddo

!
! Create and apply adjoint filter to extended 3D variables
!

         W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)

       do icol=7,1,-1
         do L=1,hz
           W(:,:,:,1-L )= W(:,:,:,1+L )
           W(:,:,:,LM+L)= W(:,:,:,LM-L)
         end do
                                                 call btim(   hfiltT_tim)
         call dibetat(km3,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
                                                 call etim(   hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_3d(W,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax)
                                                 call etim(    bocoT_tim)
       enddo

     if(l_hgen)  then
          H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
     endif

     do icol=7,1,-1
                                                 call btim(   hfiltT_tim)
       if(l_hgen)  then
         do L=1,hz
           H(:,:,:,1-L )= H(:,:,:,1+L )
           H(:,:,:,LM+L)= H(:,:,:,LM-L)
         end do

         call dibetat(km3,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                     ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
       endif
                                                 call etim(   hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_3d(H,km3,im,jm,Lm,hx,hy,hz,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)
     enddo


!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

       VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
       call composite_to_stack(VM2D,VM3D,VALL)

     if(l_hgen)  then
       HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
       call composite_to_stack(HM2D,HM3D,HALL)
     endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***
!
! From single stacked to composite variables
!

       call stack_to_composite(VALL,VM2D,VM3D)
     if(l_hgen)  then
       call stack_to_composite(HALL,HM2D,HM3D)
     endif



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
!  Apply filter to 2D variables first
!
       do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(VM2D,km2,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim(   hfilt_tim)
         call dibeta(km2,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VM2D, ff, iout,jout)
                                                 call etim(   hfilt_tim)
       enddo

     do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(HM2D,km2,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
                                                 call btim(   hfilt_tim)
       if(l_hgen)  then
         call dibeta(km2,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HM2D, ff, iout,jout)
       endif
                                                 call etim(   hfilt_tim)
     enddo

!
! Create and apply filter to extended 3D variables
!

           W(:,:,:,1:lm)=VM3D(:,:,:,1:lm)
        do L=1,hz
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j, 1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do

       do icol=1,7
                                                 call btim( boco_tim)
         call boco_3d(W,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax)
                                                 call etim( boco_tim)
                                                 call btim(   hfilt_tim)
         call dibeta(km3,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, W, ff, iout,jout,lout)
                                                 call etim(   hfilt_tim)
        enddo 

     if(l_hgen)  then
           H(:,:,:,1:lm)=HM3D(:,:,:,1:lm)
        do L=1,hz
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
              H(:,i,j,1-L )=H(:,i,j, 1+L)
              H(:,i,j,LM+L)=H(:,i,j,LM-L)
          end do
          end do
        end do
     endif
       do icol=1,7
                                                 call btim( boco_tim)
         call boco_3d(H,km3,im,jm,lm,hx,hy,hz,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
                                                 call btim(   hfilt_tim)
         if(l_hgen)  then
         call dibeta(km3,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, 1-hz,1,lm,lm+hz,icol, nfil  &
                    ,qcols,dixs3,diys3,dizs3,JCOL,vpasp3, H, ff, iout,jout,lout)
         endif
                                                 call etim(   hfilt_tim)
       enddo

!
! Go back from extended 3D variables and combine them with 2D variables in one stacked variable
!

       VM3D(:,:,:,1:lm)= W(:,:,:,1:lm)
       call composite_to_stack(VM2D,VM3D,VALL)

     if(l_hgen)  then
       HM3D(:,:,:,1:lm)=H(:,:,:,1:lm)
       call composite_to_stack(HM2D,HM3D,HALL)
     endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)


!-----------------------------------------------------------------------

deallocate(VM3D)
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

deallocate(W)
deallocate(H)

deallocate(JCOL)

!-----------------------------------------------------------------------
                        endsubroutine filtering_lin3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad2_bckg
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Apply vertical filter before and after horizontal              ! 
!     - 2d radial filter                                               !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
use mg_intstate, only: VALL,HALL
implicit none

integer(i_kind) L,i,j
!-----------------------------------------------------------------------


!***
!*** Adjoint of beta filter in vertical direction
!***
                                                 call btim(    vfiltT_tim)
if(l_vertical_filter) then

      call sup_vrbeta1T_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1, VALL)

endif
                                                 call etim(    vfiltT_tim)

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)


!***
!*** Apply adjoint of Beta filter at all generations 
!***

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call btim(    hfiltT_tim)

      call rbetaT(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbetaT(km,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
  endif
                                                 call etim(    hfiltT_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call btim(    bocoT_tim)

        call bocoT_2d(VALL,km,im,jm,hx,hy)
        call bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(    bocoT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call weighting_all(VALL,HALL,lhelm)

                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( boco_tim)

      call boco_2d(VALL,km,im,jm,hx,hy)
      call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(    boco_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
                                                 call btim(    hfilt_tim)

      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
  endif

                                                 call etim(    hfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)
!***
!*** Apply beta filter in vertical direction
!***

                                                 call btim(    vfilt_tim)
if(l_vertical_filter) then

      call sup_vrbeta1_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VALL)

endif
                                                 call etim(    vfilt_tim)

!-----------------------------------------------------------------------
                        endsubroutine filtering_rad2_bckg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_lin2_bckg
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 8:                                     !
!                                                                      !
!     - Apply vertical filter before and after horizontal              !
!     - 2d line filter
!                                                                      !
!***********************************************************************
use mg_parameter, only: nfil
use mg_intstate, only: dixs,diys,hss2
use mg_intstate, only: VALL,HALL
use mg_intstate, only: pasp1,ss1
implicit none

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff

!-----------------------------------------------------------------------
!***
!*** Adjoint of beta filter in vertical direction
!***
                                                 call btim(    vfiltT_tim)
if(l_vertical_filter) then

      call sup_vrbeta1T_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1, VALL)

endif
                                                 call etim(    vfiltT_tim)

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend (Step 1)
!***
     
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)

!***
!*** Apply adjoint of Beta filter at all generations 
!***


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontal
!
       do icol=3,1,-1
                                                 call btim(    hfiltT_tim)
         call dibetat(km,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_2d(VALL,km,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)
       enddo

     do icol=3,1,-1
                                                 call btim(    hfiltT_tim)
       if(l_hgen)  then
         call dibetat(km,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
         call bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)
     enddo

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call weighting_all(VALL,HALL,lhelm)

                                                call etim( weight_tim)

!***
!*** Apply Beta filter at all generations
!***



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontal
!
       do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(VALL,km,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
         call dibeta(km,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
                                                 call etim( hfilt_tim)
       enddo

     do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
       if(l_hgen)  then
         call dibeta(km,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
       endif
                                                 call etim( hfilt_tim)
     enddo

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL,lquart)

                                                 call etim(   dnsend_tim)
!***
!*** Apply beta filter in vertical direction
!***

                                                 call btim(    vfilt_tim)
if(l_vertical_filter) then

      call sup_vrbeta1_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VALL)

endif
                                                 call etim(    vfilt_tim)

!-----------------------------------------------------------------------
                        endsubroutine filtering_lin2_bckg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_fast_bckg
!***********************************************************************
!                                                                      !
! Fast istotropoic multigrid filtering procedure:                      !
!                                                                      !
!     - Apply adjoint of vertical filter before and directec vertical  !
!       filter after horizontal                                        !
!     - 1d+1d horizontal filter                                        !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,paspx,paspy,ss1,ssx,ssy
use mg_intstate, only: VALL,HALL
implicit none

integer(i_kind) L,i,j
!-----------------------------------------------------------------------
!***
!*** Adjoint of beta filter in vertical direction
!***
                                                 call btim(    vfiltT_tim)
if(l_vertical_filter) then

      call sup_vrbeta1T_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1, VALL)

endif
                                                 call etim(    vfiltT_tim)

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
                                                 call btim( upsend_tim)
       call upsending_all(VALL,HALL,lquart)
                                                 call etim( upsend_tim)

!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Horizontally
!

                                                 call btim(    hfiltT_tim)
    do j=1,jm
      call rbetaT(km,hx,1,im,paspx,ssx,VALL(:,:,j))
    enddo
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
       call bocoTx(VALL,km,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)

                                                 call btim(    hfiltT_tim)
    do i=1,im
      call rbetaT(km,hy,1,jm,paspy,ssy,VALL(:,i,:))
    enddo
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
      call bocoTy(VALL,km,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)

                                                 call btim(    hfiltT_tim)
  if(l_hgen)  then
    do j=1,jm
      call rbetaT(km,hx,1,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
      call bocoTx(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)

                                                 call btim(    hfiltT_tim)
  if(l_hgen)  then
    do i=1,im
      call rbetaT(km,hy,1,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
      call bocoTy(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)
      call weighting_all(VALL,HALL,lhelm)
                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
! Horizonatally

                                                 call btim( boco_tim)
      call bocoy(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
  if(l_hgen)  then
    do i=1,im
      call rbeta(km,hy,1,jm,paspy,ssy,HALL(:,i,:))
    enddo
  endif
                                                 call etim( hfilt_tim)

                                                 call btim( boco_tim)
      call bocox(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
  if(l_hgen)  then
    do j=1,jm
      call rbeta(km,hx,1,im,paspx,ssx,HALL(:,:,j))
    enddo
  endif
                                                 call etim( hfilt_tim)

                                                 call btim( boco_tim)
      call bocoy(VALL,km,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
    do i=1,im
      call rbeta(km,hy,1,jm,paspy,ssy,VALL(:,i,:))
    enddo
                                                 call etim( hfilt_tim)

                                                 call btim( boco_tim)
     call bocox(VALL,km,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
    do j=1,jm
      call rbeta(km,hx,1,im,paspx,ssx,VALL(:,:,j))
    enddo
                                                 call etim( hfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       call downsending_all(HALL,VALL,lquart)
                                                 call etim(   dnsend_tim)

!***
!*** Apply beta filter in vertical direction
!***
                                                 call btim(    vfilt_tim)

if(l_vertical_filter) then

       call sup_vrbeta1_bck(km,km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VALL)

endif

                                                 call etim(    vfilt_tim)

!-----------------------------------------------------------------------
                        endsubroutine filtering_fast_bckg

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad2_ens
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 11 for ensemble                        !
!                                                                      !
!     - Apply vertical filter before and after horizontal              !
!     - 2d radial filter                                               !
!     - Version for localization of ensemble                           !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
use mg_intstate, only: VALL,HALL
use mg_parameter, only: km3_all,km_all,l_filt_g1
implicit none

integer(i_kind) L,i,j
!-----------------------------------------------------------------------


!***
!*** Adjoint of beta filter in vertical direction
!***

                                                 call btim(    vfiltT_tim)
if(l_vertical_filter) then

        call sup_vrbeta1T_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif

                                                 call etim(    vfiltT_tim)

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       if(lquart) then
         call upsending2_ens(VALL,HALL,km_all)
       else
         call upsending_ens(VALL,HALL,km_all)
       endif
                                                 call etim( upsend_tim)


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    hfiltT_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

  if(l_filt_g1) then
      call rbetaT(km_all,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
  endif
  if(l_hgen)  then
      call rbetaT(km_all,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    hfiltT_tim)


                                                 call btim(     bocoT_tim)
  if(l_filt_g1) then
        call bocoT_2d(VALL,km_all,im,jm,hx,hy)
  endif
        call bocoT_2d(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(     bocoT_tim)


!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call weighting_ens(VALL,HALL,km_all)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations 
!***

                                                 call btim(     boco_tim)

  if(l_filt_g1) then
      call boco_2d(VALL,km_all,im,jm,hx,hy)
  endif
      call boco_2d(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)

                                                 call etim(     boco_tim)


                                                 call btim( hfilt_tim)

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

  if(l_filt_g1) then
      call rbeta(km_all,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
  endif
  if(l_hgen)  then
      call rbeta(km_all,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    hfilt_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       if(lquart) then
         call downsending2_ens(HALL,VALL,km_all)
       else 
         call downsending_ens(HALL,VALL,km_all)
       endif
                                                 call etim(   dnsend_tim)

!***
!*** Apply beta filter in vertical direction
!***
                                                 call btim(    vfilt_tim)
if(l_vertical_filter) then

        call sup_vrbeta1_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif
                                                 call etim(    vfilt_tim)


!-----------------------------------------------------------------------
                        endsubroutine filtering_rad2_ens

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_lin2_ens
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure                                        !
!                                                                      !
!     - Vertical filter before and after horizontal                    !
!     - Line filters in horizontal                                     !
!     - Version for localization of ensemble                           !
!                                                                      !
!***********************************************************************
use mg_parameter, only: nfil
use mg_parameter, only: km3_all,km_all,l_filt_g1
use mg_intstate, only: dixs,diys,hss2
use mg_intstate, only: VALL,HALL
use mg_intstate, only: pasp1,ss1
implicit none

integer(i_kind) L,i,j
integer(i_kind) icol,iout,jout
logical:: ff

!----------------------------------------------------------------------
!***
!*** Adjoint of beta filter in vertical direction
!***
                                                 call btim(    vfiltT_tim)

if(l_vertical_filter) then

       call sup_vrbeta1T_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif

                                                 call etim(    vfiltT_tim)

!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
      if(lquart) then
         call upsending2_ens(VALL,HALL,km_all)
       else
         call upsending_ens(VALL,HALL,km_all)
       endif
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Adjoint filtering
!

    if(l_filt_g1) then
       do icol=3,1,-1
                                                 call btim( hfiltT_tim)
         call dibetat(km_all,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
                                                 call etim( hfiltT_tim)
                                                 call btim( bocoT_tim)
         call bocoT_2d(VALL,km_all,im,jm,hx,hy)
                                                 call etim(  bocoT_tim)
       enddo
    endif

     do icol=3,1,-1
       if(l_hgen)  then
                                                 call btim( hfiltT_tim)
         call dibetat(km_all,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                      dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
                                                 call etim( hfiltT_tim)
       endif
                                                 call btim(    bocoT_tim)
         call bocoT_2d(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)
      call weighting_ens(VALL,HALL,km_all)
                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations
!***



!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!
    if(l_filt_g1) then
       do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(VALL,km_all,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
         call dibeta(km_all,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), VALL, ff, iout,jout)
                                                 call etim( hfilt_tim)
       enddo
    endif

     do icol=1,3
                                                 call btim( boco_tim)
         call boco_2d(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
       if(l_hgen)  then
                                                 call btim( hfilt_tim)
         call dibeta(km_all,1-hx,1,im,im+hx, 1-hy,1,jm,jm+hy, nfil,  &
                     dixs(:,:,icol),diys(:,:,icol),hss2(:,:,icol), HALL, ff, iout,jout)
                                                 call etim( hfilt_tim)
       endif
     enddo


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!***
!***  Downsend, interpolate and add, then zero high generations 
!***

                                                 call btim(   dnsend_tim)
       if(lquart) then
         call downsending2_ens(HALL,VALL,km_all)
       else
         call downsending_ens(HALL,VALL,km_all)
       endif

                                                 call etim(   dnsend_tim)

!***
!*** Apply beta filter in vertical direction
!***
                                                 call btim(    vfilt_tim)
if(l_vertical_filter) then

       call sup_vrbeta1T_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif

                                                 call etim(    vfilt_tim)



!-----------------------------------------------------------------------
                        endsubroutine filtering_lin2_ens

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_fast_ens
!***********************************************************************
!                                                                      !
! Fast multigrid filtering procedure:                                  !
!                                                                      !
!     - Apply vertical filter before and after horizontal              !
!     - 1d+1d horizontal filter + 1d vertical filter                   !
!     - Version fo localizaiton of ensemble                            !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,paspx,paspy,ss1,ssx,ssy
use mg_intstate, only: VALL,HALL
use mg_parameter, only: km3_all,km_all,l_filt_g1
use mg_parameter, only: km3,n_ens
implicit none

integer(i_kind) L,i,j
!-----------------------------------------------------------------------
!***
!*** Adjoint of beta filter in vertical direction
!***
                                                 call btim(    vfiltT_tim)
if(l_vertical_filter) then

        call sup_vrbeta1T_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif

                                                 call etim(    vfiltT_tim)


!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
       if(lquart) then
         call upsending2_ens(VALL,HALL,km_all)
       else
         call upsending_ens(VALL,HALL,km_all)
       endif
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------


!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

  if(l_filt_g1) then
                                                 call btim(    hfiltT_tim)
    do j=1,jm
      call rbetaT(km_all,hx,1,im,paspx,ssx,VALL(:,:,j))
    enddo
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
       call bocoTx(VALL,km_all,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)

                                                 call btim(    hfiltT_tim)
    do i=1,im
      call rbetaT(km_all,hy,1,jm,paspy,ssy,VALL(:,i,:))
    enddo
                                                 call etim(    hfiltT_tim)
                                                 call btim(    bocoT_tim)
      call bocoTy(VALL,km_all,im,jm,hx,hy)
                                                 call etim(    bocoT_tim)
  endif


  if(l_hgen)  then
                                                 call btim(    hfiltT_tim)
    do j=1,jm
      call rbetaT(km_all,hx,1,im,paspx,ssx,HALL(:,:,j))
    enddo
                                                 call etim(    hfiltT_tim)
  endif
                                                 call btim(    bocoT_tim)
      call bocoTx(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)

  if(l_hgen)  then
                                                 call btim(    hfiltT_tim)
    do i=1,im
      call rbetaT(km_all,hy,1,jm,paspy,ssy,HALL(:,i,:))
    enddo
                                                 call etim(    hfiltT_tim)
  endif
                                                 call btim(    bocoT_tim)
      call bocoTy(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim(    bocoT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)
      call weighting_ens(VALL,HALL,km_all)
                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations 
!***


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

  if(l_filt_g1) then

                                                 call btim( boco_tim)
      call bocoy(VALL,km_all,im,jm,hx,hy)
                                                 call etim( boco_tim)

                                                 call btim( hfilt_tim)
    do i=1,im
      call rbeta(km_all,hy,1,jm,paspy,ssy,VALL(:,i,:))
    enddo
                                                 call etim( hfilt_tim)

                                                 call btim( boco_tim)
     call bocox(VALL,km_all,im,jm,hx,hy)
                                                 call etim( boco_tim)
                                                 call btim( hfilt_tim)
    do j=1,jm
      call rbeta(km_all,hx,1,im,paspx,ssx,VALL(:,:,j))
    enddo
                                                 call etim( hfilt_tim)
  endif


                                                 call btim( boco_tim)
      call bocoy(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
  if(l_hgen)  then
                                                 call btim( hfilt_tim)
    do i=1,im
      call rbeta(km_all,hy,1,jm,paspy,ssy,HALL(:,i,:))
    enddo
                                                 call etim( hfilt_tim)
  endif

                                                 call btim( boco_tim)
      call bocox(HALL,km_all,im,jm,hx,hy,Fimax,Fjmax,2,gm)
                                                 call etim( boco_tim)
       
  if(l_hgen)  then
                                                 call btim( hfilt_tim)
    do j=1,jm
      call rbeta(km_all,hx,1,im,paspx,ssx,HALL(:,:,j))
    enddo
                                                 call etim( hfilt_tim)
  endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)
       if(lquart) then
         call downsending2_ens(HALL,VALL,km_all)
       else
         call downsending_ens(HALL,VALL,km_all)
       endif
                                                 call etim(   dnsend_tim)
!***
!*** Apply beta filter in vertical direction
!***
                                                 call btim(    vfilt_tim)
if(l_vertical_filter) then

        call sup_vrbeta1_ens(km3_all,hx,hy,hz,im,jm,lm, pasp1,ss1,VALL)

endif

                                                 call etim(    vfilt_tim)

!-----------------------------------------------------------------------
                        endsubroutine filtering_fast_ens


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad2_z_loc_g3
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure for localization                       !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - Apply vertical filter before and after horizontal              !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter                                               !
!     - 3 generatons                                                   !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
use mg_intstate, only: VALL,HALL
use mg_bocos, only: boco_2d_loc,bocoT_2d_loc
use mg_generations, only: upsending_loc,downsending_loc,weighting_loc
use mg_parameter, only: l_loc_vertical,n_ens,km_all,km3_all
use mg_parameter, only: km_4,km_16
implicit none


real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:):: H04
real(r_kind), allocatable, dimension(:,:,:):: H16

integer(i_kind) L,i,j
!-----------------------------------------------------------------------

allocate(VM3D(n_ens*km3,1-hx:im+hx,1-hy:jm+hy,lm))            ; VM3D=0.

allocate(H04(km_4 ,1-hx:im+hx,1-hy:jm+hy   ))                 ; H04=0.
allocate(H16(km_16,1-hx:im+hx,1-hy:jm+hy   ))                 ; H16=0.
!----------------------------------------------------------------------

!***
!*** Adjoint of beta filter in vertical direction  << Need modification >>>
!***

   if(l_loc_vertical) then
      call S2C_ens(VALL,VM3D,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
        call sup_vrbeta1T(km3_all,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call C2S_ens(VM3D,VALL,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
   endif


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>
   if(mype==0) then
       print *, 'Before upsending_loc'
   endif 
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>

       call upsending_loc(VALL,H04,H16,km,km_4,km_16)

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>
   if(mype==0) then
       print *, 'After upsending_loc'
   endif 
     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>
                                                 call etim( upsend_tim)
!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bocoT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!T      call rbetaT(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))

      call rbetaT(km_4 ,hx,1,im,hy,1,jm,pasp2,ss2,H04(:,:,:))
      call rbetaT(km_16,hx,1,im,hy,1,jm,pasp2,ss2,H16(:,:,:))

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

        call bocoT_2d_loc(H04,km_4 ,im,jm,hx,hy,Fimax(2),Fjmax(2),2)
        call bocoT_2d_loc(H16,km_16,im,jm,hx,hy,Fimax(2),Fjmax(2),3)


                                                 call etim(    bocoT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

!      call weighting_all(VALL,HALL,lhelm)

      call weighting_loc(VALL,H04,H16,km,km_4,km_16)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( boco_tim)

!T      call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
!T      call boco_2d(VALL,km,im,jm,hx,hy)

        call boco_2d_loc(H16,km_16,im,jm,hx,hy,Fimax(2),Fjmax(2),3)
        call boco_2d_loc(H04,km_4 ,im,jm,hx,hy,Fimax(2),Fjmax(2),2)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

!T      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
!  if(l_hgen)  then
!      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
!T  endif

      call rbeta(km_16,hx,1,im,hy,1,jm,pasp2,ss2,H16(:,:,:))
      call rbeta(km_4 ,hx,1,im,hy,1,jm,pasp2,ss2,H04(:,:,:))


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    boco_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)

       call downsending_loc(H16,H04,VALL,km,km_4,km_16)

                                                 call etim(   dnsend_tim)
!***
!*** Apply beta filter in vertical direction   << Need modification >>>
!***
   if(l_loc_vertical) then
      call S2C_ens(VALL,VM3D,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
      call sup_vrbeta1(km3_all,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call C2S_ens(VM3D,VALL,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
   endif


deallocate(VM3D) 

deallocate(H04)
deallocate(H16)

!-----------------------------------------------------------------------
                        endsubroutine filtering_rad2_z_loc_g3

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad2_z_loc_g4
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure for localization                       !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - Apply vertical filter before and after horizontal              !
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter                                               !
!     - 4 generatons                                                   !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
use mg_intstate, only: VALL,HALL
use mg_bocos, only: boco_2d_loc,bocoT_2d_loc
use mg_generations, only: upsending_loc,downsending_loc,weighting_loc
use mg_parameter, only: l_loc_vertical,n_ens,km_all,km3_all
use mg_parameter, only: km_4,km_16,km_64
implicit none


real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:):: H04
real(r_kind), allocatable, dimension(:,:,:):: H16
real(r_kind), allocatable, dimension(:,:,:):: H64

integer(i_kind) L,i,j
!-----------------------------------------------------------------------

allocate(VM3D(km3*n_ens,1-hx:im+hx,1-hy:jm+hy,lm))            ; VM3D=0.

allocate(H04(km_4 ,1-hx:im+hx,1-hy:jm+hy   ))                 ; H04=0.
allocate(H16(km_16,1-hx:im+hx,1-hy:jm+hy   ))                 ; H16=0.
allocate(H64(km_64,1-hx:im+hx,1-hy:jm+hy   ))                 ; H64=0.
!----------------------------------------------------------------------

!***
!*** Adjoint of beta filter in vertical direction  << Need modification >>>
!***

   if(l_loc_vertical) then
      call S2C_ens(VALL,VM3D,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
        call sup_vrbeta1T(km3_all,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call C2S_ens(VM3D,VALL,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
   endif


!==================== Adjoint (Conservative step) ======================

!***
!*** Adjoint interpolate and upsend 
!***
     
                                                 call btim( upsend_tim)

       call upsending_loc(VALL,H04,H16,H64,km,km_4,km_16,km_64)

                                                 call etim( upsend_tim)
!----------------------------------------------------------------------

!***
!*** Apply adjoint of Beta filter at all generations 
!***
                                                 call btim(    bocoT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!T      call rbetaT(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))

      call rbetaT(km_4 ,hx,1,im,hy,1,jm,pasp2,ss2,H04(:,:,:))
      call rbetaT(km_16,hx,1,im,hy,1,jm,pasp2,ss2,H16(:,:,:))
      call rbetaT(km_64,hx,1,im,hy,1,jm,pasp2,ss2,H64(:,:,:))

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

        call bocoT_2d_loc(H04,km_4 ,im,jm,hx,hy,Fimax(2),Fjmax(2),2)
        call bocoT_2d_loc(H16,km_16,im,jm,hx,hy,Fimax(2),Fjmax(2),3)
        call bocoT_2d_loc(H64,km_64,im,jm,hx,hy,Fimax(2),Fjmax(2),4)


                                                 call etim(    bocoT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

!      call weighting_all(VALL,HALL,lhelm)

      call weighting_loc(VALL,H04,H16,H64,km,km_4,km_16,km_64)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( boco_tim)

!T      call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)
!T      call boco_2d(VALL,km,im,jm,hx,hy)

        call boco_2d_loc(H64,km_64,im,jm,hx,hy,Fimax(2),Fjmax(2),4)
        call boco_2d_loc(H16,km_16,im,jm,hx,hy,Fimax(2),Fjmax(2),3)
        call boco_2d_loc(H04,km_4 ,im,jm,hx,hy,Fimax(2),Fjmax(2),2)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

!T      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,VALL(:,:,:))
!  if(l_hgen)  then
!      call rbeta(km,hx,1,im,hy,1,jm,pasp2,ss2,HALL(:,:,:))
!T  endif

      call rbeta(km_64,hx,1,im,hy,1,jm,pasp2,ss2,H64(:,:,:))
      call rbeta(km_16,hx,1,im,hy,1,jm,pasp2,ss2,H16(:,:,:))
      call rbeta(km_4 ,hx,1,im,hy,1,jm,pasp2,ss2,H04(:,:,:))


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    boco_tim)

!***
!***  Downsend, interpolate and add (Step 4)
!***  Then zero high generations (Step 5)
!***

                                                 call btim(   dnsend_tim)

       call downsending_loc(H64,H16,H04,VALL,km,km_4,km_16,km_64)

                                                 call etim(   dnsend_tim)
!***
!*** Apply beta filter in vertical direction   << Need modification >>>
!***

   if(l_loc_vertical) then
      call S2C_ens(VALL,VM3D,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
      call sup_vrbeta1(km3_all,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call C2S_ens(VM3D,VALL,1-hx,im+hx,1-hy,jm+hy,lm,km,km_all)
   endif


deallocate(VM3D) 

deallocate(H04)
deallocate(H16)
deallocate(H64)

!-----------------------------------------------------------------------
                        endsubroutine filtering_rad2_z_loc_g4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                          *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=1,jm
        do i=1,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L)=W(:,1+L)
            W(:,LM+L)=W(:,LM-L)
          end do
             call rbeta(kmax,hz,1,lm,  pasp,ss,W)
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

  
!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1T                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1T                                         *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm,  pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

        do j=1,jm
        do i=1,im
          do L=1,Lm
            W(:,L)=V(:,i,j,L)
          end do
          do L=1,hz
            W(:,1-L )=W(:,1+L )
            W(:,LM+L)=W(:,LM-L)
          end do
             call rbetaT(kmax,hz,1,lm, pasp,ss,W)
!
! Apply adjoint at the edges of domain
!
          do L=1,hz
            W(:,1+L)=W(:,1+L)+W(:,1-L)
            W(:,LM-L)=W(:,LM-L)+W(:,LM+L)
          enddo
          do l=1,Lm
            V(:,i,j,L)=W(:,L)
          end do
        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3                        &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss, V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,1:im,1:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(1:im,1:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,L

!----------------------------------------------------------------------

          do L=1,Lm
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j,1+L )
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call rbeta(kmax,hx,1,im, hy,1,jm, hz,1,lm, pasp,ss,W)

  
          do l=1,Lm
          do j=1,jm
          do i=1,im
            V(:,i,j,L)=W(:,i,j,L)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta3T                       &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta3                                           *
!                                                                     *
!**********************************************************************
(kmax,hx,hy,hz,im,jm,lm, pasp,ss,V)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: kmax,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1:lm),intent(inout):: V
real(r_kind),dimension(3,3,1:im,1:jm,1:lm), intent(in):: pasp
real(r_kind),dimension(1:im,1:jm,1:lm), intent(in):: ss

real(r_kind),dimension(1:kmax,1-hx:im+hx,1-hy:jm+hy,1-hz:lm+hz):: W

integer(i_kind):: i,j,l

!----------------------------------------------------------------------

          do L=1,Lm
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
            W(:,i,j,L)=V(:,i,j,L)
          end do
          end do
          end do

        do L=1,hz
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
              W(:,i,j,1-L )=W(:,i,j, 1+L)
              W(:,i,j,LM+L)=W(:,i,j,LM-L)
          end do
          end do
        end do
    
    
           call rbetaT(kmax,hx,1,im, hy,1,jm, hz,1,lm, pasp,ss,W)

!
! Apply adjoint at the edges of domain
!
        do L=1,hz
          do j=1-hy,jm+hy
          do i=1-hx,im+hx
              W(:,i,j,1+L )=W(:,i,j, 1+L)+W(:,i,j, 1-L)
              W(:,i,j,LM-L)=W(:,i,j,LM-L)+W(:,i,j,LM+L)
          end do
          end do
         end do
  
          do l=1,lm
          do j=1,jm
          do i=1,im
            V(:,i,j,l)=W(:,i,j,l)
          end do
          end do
          end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta3T

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1_ens                     &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                           *
!                                                                     *
!**********************************************************************
(km_en,hx,hy,hz,im,jm,lm, pasp,ss, VALL)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: km_en,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:km_en*lm,1-hx:im+hx,1-hy:jm+hy),intent(inout):: VALL
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:km_en,1-hz:lm+hz):: W

integer(i_kind):: i,j,L,k,k_ind,kloc

!----------------------------------------------------------------------

    do j=1,jm
    do i=1,im
      do k=1,km_en
        k_ind =(k-1)*Lm
          do L=1,Lm
            kloc=k_ind+L
            W(k,L)=VALL(kloc,i,j)
          end do
      enddo
          do L=1,hz
            W(:,1-L )=W(:,1+L )
            W(:,LM+L)=W(:,LM-L)
          end do

             call rbeta(km_en,hz,1,lm,  pasp,ss,W)

      do k=1,km_en
        k_ind =(k-1)*Lm
          do L=1,Lm
            kloc=k_ind+L
            VALL(kloc,i,j)= W(k,L)
          end do
      enddo
   enddo
   enddo

  
!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1_ens

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1T_ens                    &
!**********************************************************************
!                                                                     *
!     Adjoint of vrbeta1_ens                                          *
!                                                                     *
!**********************************************************************
(km_en,hx,hy,hz,im,jm,lm,  pasp,ss, VALL)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: km_en,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:km_en*lm,1-hx:im+hx,1-hy:jm+hy),intent(inout):: VALL
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:km_en,1-hz:lm+hz):: W

integer(i_kind):: i,j,L,k,k_ind,kloc

!----------------------------------------------------------------------

        do j=1,jm
        do i=1,im

          do k=1,km_en
            k_ind = (k-1)*Lm
              do L=1,Lm
                kloc=k_ind+L
                W(k,L)=VALL(kloc,i,j)
              end do
          enddo
             do L=1,hz
                W(:,1-L )=W(:,1+L )
                W(:,LM+L)=W(:,LM-L)
             end do

               call rbetaT(km_en,hz,1,lm, pasp,ss,W)
!
! Apply adjoint at the edges of domain
!
             do L=1,hz
               W(:,1+L )=W(:,1+L )+W(:,1-L)
               W(:,LM-L)=W(:,LM-L)+W(:,LM+L)
             enddo

          do k=1,km_en
            k_ind = (k-1)*Lm
              do l=1,Lm
                kloc=k_ind+L
                VALL(kloc,i,j)=W(k,L)
              enddo
          end do

        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T_ens

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1_bck                     &
!**********************************************************************
!                                                                     *
!     conversion of vrbeta1                                           *
!                                                                     *
!**********************************************************************
(km,km3,hx,hy,hz,im,jm,lm, pasp,ss, VALL)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: km,km3,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:km,1-hx:im+hx,1-hy:jm+hy),intent(inout):: VALL
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:km3,1-hz:lm+hz):: W

integer(i_kind):: i,j,L,k,k_ind,kloc

!----------------------------------------------------------------------

    do j=1,jm
    do i=1,im
      do k=1,km3
        k_ind =(k-1)*Lm
          do L=1,Lm
            kloc=k_ind+L
            W(k,L)=VALL(kloc,i,j)
          end do
      enddo
          do L=1,hz
            W(:,1-L )=W(:,1+L )
            W(:,LM+L)=W(:,LM-L)
          end do

             call rbeta(km3,hz,1,lm,  pasp,ss,W)

      do k=1,km3
        k_ind =(k-1)*Lm
          do L=1,Lm
            kloc=k_ind+L
            VALL(kloc,i,j)= W(k,L)
          end do
      enddo
   enddo
   enddo

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1_bck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine sup_vrbeta1T_bck                    &
!**********************************************************************
!                                                                     *
!     Adjoint of vrbeta1_bck                                          *
!                                                                     *
!**********************************************************************
(km,km3,hx,hy,hz,im,jm,lm,pasp,ss, VALL)
!----------------------------------------------------------------------
implicit none

integer(i_kind),intent(in):: km,km3,hx,hy,hz,im,jm,lm
real(r_kind),dimension(1:km,1-hx:im+hx,1-hy:jm+hy),intent(inout):: VALL
real(r_kind),dimension(1,1,1:lm), intent(in):: pasp
real(r_kind),dimension(1:lm), intent(in):: ss

real(r_kind),dimension(1:km3,1-hz:lm+hz):: W

integer(i_kind):: i,j,L,k,k_ind,kloc

!----------------------------------------------------------------------

        do j=1,jm
        do i=1,im

          do k=1,km3
            k_ind = (k-1)*Lm
              do L=1,Lm
                kloc=k_ind+L
                W(k,L)=VALL(kloc,i,j)
              end do
          enddo
             do L=1,hz
                W(:,1-L )=W(:,1+L )
                W(:,LM+L)=W(:,LM-L)
             end do

               call rbetaT(km3,hz,1,lm, pasp,ss,W)
!
! Apply adjoint at the edges of domain
!
             do L=1,hz
               W(:,1+L )=W(:,1+L )+W(:,1-L)
               W(:,LM-L)=W(:,LM-L)+W(:,LM+L)
             enddo

          do k=1,km3
            k_ind = (k-1)*Lm
              do l=1,Lm
                kloc=k_ind+L
                VALL(kloc,i,j)=W(k,L)
              enddo
          end do

        end do
        end do

!----------------------------------------------------------------------
                        endsubroutine sup_vrbeta1T_bck


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_filtering
