!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt_all
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
real(r_kind),allocatable,dimension(:,:,:):: A3D
real(r_kind),allocatable,dimension(:,:):: A2D

real(r_kind),allocatable,dimension(:,:,:):: F3D

!----------------------------------------------------------------------

allocate(A3D(km3,1:nm,1:mm,lm_a))
allocate(A2D(km2,1:nm,1:mm)

allocate(F3D(km3,1:nm,1:mm,lm))


     call stack_to_composite_spec(WORKA,A2D,A3D)

  if(lm_a>lm) then
     call lwq_vertical_adjoint_spec(km3,lm_a,lm,1,nm,1,mm,               &
                                    cvf1,cvf2,cvf3,cvf3,lref,A3D,F3D)
  else 

    do L=1,lm
      F3D(:,:,:,L)=A3D(:,:,:,L)
    enddo

  endif

     call composite_to_stack_spec(A2D,F3D,WORK)

     

     call anal_to_filt


!----------------------------------------------------------------------
                        endsubroutine anal_to_filt_all
