!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine filt_to_anal_all
!***********************************************************************
!                                                                      !
!  Transfer data from filter to analysis grid                          !
!                                                                      !
!***********************************************************************
real(r_kind),allocatable,dimension(:,:,:):: A3D
real(r_kind),allocatable,dimension(:,:):: A2D

real(r_kind),allocatable,dimension(:,:,:):: F3D

!----------------------------------------------------------------------

allocate(A3D(km3,1:nm,1:mm,lm_a))
allocate(A2D(km2,1:nm,1:mm)

allocate(F3D(km3,1:nm,1:mm,lm))


    
    call filt_to_anal
   
    call stack_to_composite_spec(WORK,A2D,F3D)
  
 if(lm_a>lm) then
    call lwq_vertical_direct_spec(km3,lm,lm_a,1,nm,1,mm,                 &
                                  cvf1,cvf2,cv43,cvf4,kref,F3D,A3D)
 else

   do L=1,lm
     A3D(:,:,:,L)=F3D(:,:,:,L)
   enddo
 
 endif  
 

    call composite_to_stack_spec(A2D,A3D,WORKA)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all
