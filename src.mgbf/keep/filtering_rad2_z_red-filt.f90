!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filtering_rad2_z_rf
!***********************************************************************
!                                                                      !
! Multigrid filtering procedure 2:                                     !
!                                                                      !
!     - Multiple of 2D and 3D variables                                !
!     - Apply vertical filter before and after horizontal
!     - 1 upsending and downsending                                    !
!     - Applicaton of Helmholtz differential operator                  !
!     - 2d radial filter                                               !
!                                                                      !
!***********************************************************************
use mg_intstate, only: pasp1,pasp2,ss1,ss2
use mg_intstate, only: VALL,HALL
implicit none

real(r_kind), allocatable, dimension(:,:,:):: VM2D
real(r_kind), allocatable, dimension(:,:,:):: HM2D
real(r_kind), allocatable, dimension(:,:,:,:):: VM3D
real(r_kind), allocatable, dimension(:,:,:,:):: HM3D

integer(i_kind) L,i,j
!-----------------------------------------------------------------------

allocate(VM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; VM3D=0.
allocate(VM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; VM2D=0.
allocate(HM3D(km3,i0-hx:im+hx,j0-hy:jm+hy,lm))                  ; HM3D=0.
allocate(HM2D(km2,i0-hx:im+hx,j0-hy:jm+hy   ))                  ; HM2D=0.

!----------------------------------------------------------------------

!***
!*** Adjoint of beta filter in vertical direction
!***
      call stack_to_composite(VALL,VM2D,VM3D)
      call sup_vrbeta1T(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)


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
                                                 call btim(    bfiltT_tim)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!T      call rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbetaT(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!T        call bocoT_2d(VALL,km,im,jm,hx,hy)
        call bocoT_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


                                                 call etim(    bfiltT_tim)
!***
!*** Apply (a-b\nabla^2)
!***

                                                call btim( weight_tim)

      call weighting_all(VALL,HALL,lhelm)


                                                call etim( weight_tim)


!***
!*** Apply Beta filter at all generations (Step 7)
!***
                                                 call btim( bfilt_tim)

!T      call boco_2d(VALL,km,im,jm,hx,hy)
      call boco_2d(HALL,km,im,jm,hx,hy,Fimax,Fjmax,2,gm)


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! Filtering
!

!T      call rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,VALL(:,:,:))
  if(l_hgen)  then
      call rbeta(km,hx,i0,im,hy,j0,jm,pasp2,ss2,HALL(:,:,:))
  endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

                                                 call etim(    bfilt_tim)

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
      call stack_to_composite(VALL,VM2D,VM3D)
      call sup_vrbeta1(km3,hx,hy,hz,im,jm,lm,pasp1,ss1,VM3D)
      call composite_to_stack(VM2D,VM3D,VALL)




deallocate(VM3D) 
deallocate(VM2D)
deallocate(HM3D)
deallocate(HM2D)

!-----------------------------------------------------------------------
                        endsubroutine filtering_rad2_z_rf
