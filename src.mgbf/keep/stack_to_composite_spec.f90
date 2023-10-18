!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite_spec              &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km_a,1:nm,1:mm)     ,intent(in):: ARR_ALL
real(r_kind),dimension(km3 ,1:nm,1:mm,lm_a),intent(out):: A3D
real(r_kind),dimension(km2 ,1:nm,1:mm)     ,intent(out):: A2D
integer(i_kind):: n,m,L,k
!----------------------------------------------------------------------
    do L=1,lm_a
      do m=1,mm
      do n=1,nm
        do k=1,km3
          A3D(k,n,m,L)=ARR_ALL( (k-1)*lm_a+L,n,m )
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          A2D(k,:,:)=ARR_ALL(km3*lm_a+k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine stack_to_composite_spec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack_spec              &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km2,1:nm,1:mm),   intent(in):: A2D
real(r_kind),dimension(km3,1:nm,1:mm,lm),intent(in):: A3D
real(r_kind),dimension(km ,1:nm,1:mm),   intent(out):: ARR_ALL
integer(i_kind):: n,m,L,k
!----------------------------------------------------------------------
    do L=1,lm_a
      do m=1,mm
      do n=1,nm
        do k=1,km3
          ARR_ALL( (k-1)*lm_a+L,n,m )=A3D(k,n,m,L)
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          ARR_ALL(km3*lm_a+k,:,:)=A2D(k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack_spec

