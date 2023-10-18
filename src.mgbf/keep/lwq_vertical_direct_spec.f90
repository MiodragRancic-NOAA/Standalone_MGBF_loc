!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_adjoint_spec            &
!***********************************************************************
!                                                                      !
!  Direct linerly weighted quadratic adjoint interpolation in vertical !
!  from reslution nm to resolution im                                  !
!                                                                      !
!              (  km <= nm )                                           !
!                                                                      !
!***********************************************************************
(km3,nm,km,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,W,F)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:nm), intent(in):: W
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:km), intent(out):: F
integer(i_kind):: k,n
!-----------------------------------------------------------------------

  f = 0.
do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    F(:,:,:,1) = F(:,:,:,1)+c2(n)*W(:,:,:,n)
    F(:,:,:,2) = F(:,:,:,2)+c3(n)*W(:,:,:,n)
    F(:,:,:,3) = F(:,:,:,3)+c4(n)*W(:,:,:,n)
  elseif &
    ( k==km-1) then
    F(:,:,:,km-2) = F(:,:,:,km-2)+c1(n)*W(:,:,:,n)
    F(:,:,:,km-1) = F(:,:,:,km-1)+c2(n)*W(:,:,:,n)
    F(:,:,:,km  ) = F(:,:,:,km  )+c3(n)*W(:,:,:,n)
  elseif( k==km) then
    F(:,:,:,k   ) = F(:,:,:,k   )+c2(n)*W(:,:,:,n)
  else
    F(:,:,:,k-1) = F(:,:,:,k-1)+c1(n)*W(:,:,:,n)
    F(:,:,:,k  ) = F(:,:,:,k  )+c2(n)*W(:,:,:,n)
    F(:,:,:,k+1) = F(:,:,:,k+1)+c3(n)*W(:,:,:,n)
    F(:,:,:,k+2) = F(:,:,:,k+2)+c4(n)*W(:,:,:,n)
  endif
enddo
    F(:,:,:,1 )=F(:,:,:,1 )+W(:,:,:,1 )
    F(:,:,:,km)=F(:,:,:,km)+W(:,:,:,nm)

!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_adjoint_spec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_direct_spec             &
!***********************************************************************
!                                                                      !
!  Linerly weighted direct quadratic interpolation in vertical         !
!  from reslouion im to resolution nm                                  !
!                                                                      !
!              (  km <= nm )                                           !
!                                                                      !
!***********************************************************************
(km3,km,nm,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,F,W)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,km,nm,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:km), intent(in):: F
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:nm), intent(out):: W
integer(i_kind):: k,n
!-----------------------------------------------------------------------

do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    W(:,:,:,n) =                   c2(n)*F(:,:,:,k)+c3(n)*F(:,:,:,k+1)+c4(n)*F(:,:,:,k+2)
  elseif &
    ( k==km-1) then
    W(:,:,:,n) =c1(n)*F(:,:,:,k-1)+c2(n)*F(:,:,:,k)+c3(n)*F(:,:,:,k+1)
  elseif &
    ( k==km)   then
    W(:,:,:,n) =                   c2(n)*F(:,:,:,k)
  else
    W(:,:,:,n) =c1(n)*F(:,:,:,k-1)+c2(n)*F(:,:,:,k)+c3(n)*F(:,:,:,k+1)+c4(n)*F(:,:,:,k+2)
  endif
enddo
    W(:,:,:,1 )=F(:,:,:,1 )
    W(:,:,:,nm)=F(:,:,:,km)
    
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct_spec
