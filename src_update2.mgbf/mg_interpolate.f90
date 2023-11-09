!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_interpolate
!***********************************************************************
!                                                                      !
!    general mapping between 2d arrays using linerly squared           !
!    interpolations                                                    !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds
use mg_parameter, only: xa0,ya0,xf0,yf0,dxa,dxf,dya,dyf                 &
                       ,nm,mm,lm                                        &
                       ,i0,j0,n0,m0                                     &
                       ,im,jm,ib,jb                                    
use mg_intstate, only: iref,jref                                        &
                      ,cx0,cx1,cx2,cx3                                  &
                      ,cy0,cy1,cy2,cy3                                  
use mg_intstate, only: irefq,jrefq                                      &
                      ,qx0,qx1,qx2                                      &
                      ,qy0,qy1,qy2
use mg_intstate, only: p_coef,q_coef
use mg_intstate, only: a_coef,b_coef

!use mpimod, only: mype
use mg_mppstuff, only: mype
use mg_mppstuff, only: finishMPI
implicit none

public lsqr_mg_coef

public l_vertical_direct_spec
public l_vertical_adjoint_spec

public l_vertical_direct_spec2
public l_vertical_adjoint_spec2

public lwq_vertical_coef
public lwq_vertical_direct
public lwq_vertical_adjoint

public lwq_vertical_direct_spec
public lwq_vertical_adjoint_spec

public def_offset_coef

public lsqr_direct_offset
public lsqr_adjoint_offset

public quad_direct_offset
public quad_adjoint_offset

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_offset_coef                         
!***********************************************************************
implicit none
                                                                      
real(r_kind):: r64,r32,r128
!-----------------------------------------------------------------------
 r64 = 1.0d0/64.0d0
 r32 = 1.0d0/32.0d0
 r128= 1.0d0/128.0d0

! p_coef =(/-3.,51,29,-3/)
! q_coef =(/-3.,19.0d0,51.0d0,-3.0d0/)
! p_coef = p_coef*r64
! q_coef = q_coef*r64

 p_coef =(/-9.,111.,29.,-3./)
 q_coef =(/-3.,29.,111.,-9./)
 p_coef = p_coef*r128 
 q_coef = q_coef*r128 

 a_coef =(/5.,30.,-3./)
 b_coef =(/-3.,30.,5./)
 a_coef=a_coef*r32 
 b_coef=b_coef*r32 
!-----------------------------------------------------------------------
                        endsubroutine def_offset_coef

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_mg_coef                            
!***********************************************************************
!                                                                      !
!   Prepare coeficients for mapping between:                           !
!   filter grid on analysis decomposition:  W(i0-ib:im+ib,j0-jb:jm+jb) !
!   and analysis grid:                      V(1:nm,1:mm)               !  
!                       - offset version -                             !
!                                                                      !
!              (  im < nm  and  jm < mm   )                            !
!                                                                      !
!***********************************************************************
implicit none
real(r_kind), dimension(1:nm):: xa
real(r_kind), dimension(1-ib:im+ib):: xf
real(r_kind), dimension(1:mm):: ya
real(r_kind), dimension(1-jb:jm+jb):: yf
integer(i_kind):: i,j,n,m
real(r_kind) x1,x2,x3,x4,x
real(r_kind) x1x,x2x,x3x,x4x
real(r_kind) rx2x1,rx3x1,rx4x1,rx3x2,rx4x2,rx4x3
real(r_kind) y1,y2,y3,y4,y
real(r_kind) y1y,y2y,y3y,y4y
real(r_kind) ry2y1,ry3y1,ry4y1,ry3y2,ry4y2,ry4y3
real(r_kind) cfl1,cfl2,cfl3,cll
real(r_kind) cfr1,cfr2,cfr3,crr
real(r_kind) x1_x,x2_x,x3_x
real(r_kind) y1_y,y2_y,y3_y
!-----------------------------------------------------------------------
!
! Initialize
!
 
   do n=1,nm
     xa(n)=xa0+dxa*(n-1)
   enddo

   do i=1-ib,im+ib
     xf(i)=xf0+dxf*(i-1)
   enddo

   do m=1,mm
     ya(m)=ya0+dya*(m-1)
   enddo

   do j=1-jb,jm+jb
     yf(j)=yf0+dyf*(j-1)
   enddo

!
! Find iref and jref
!
   do n=1,nm
     do i=1-ib,im+ib-1
       if( xa(n)< xf(i)) then
         iref(n)=i-2
         irefq(n)=i-1
         exit
       endif
     enddo
   enddo

!TEST
!   do n=1,nm
!     if(mype==0) then
!       print *,'irefq(n)=',irefq
!     endif
!   enddo
!     call finishMPI
!TEST

   do m=1,mm
     do j=1-jb,jm+jb-1
       if(ya(m) < yf(j)) then
         jref(m)=j-2
         jrefq(m)=j-1
         exit
       endif
     enddo
   enddo


   do n=1,nm
     i=iref(n)
     x1=xf(i)
     x2=xf(i+1)
     x3=xf(i+2)
     x4=xf(i+3)
     x = xa(n)
       x1x = x1-x   
       x2x = x2-x   
       x3x = x3-x   
       x4x = x4-x   
       rx2x1 = 1./(x2-x1)
       rx3x1 = 1./(x3-x1)
       rx4x1 = 1./(x4-x1)
       rx3x2 = 1./(x3-x2)
       rx4x2 = 1./(x4-x2)
       rx4x3 = 1./(x4-x3)
     CFL1 = x2x*x3x*rx2x1*rx3x1
     CFL2 =-x1x*x3x*rx2x1*rx3x2
     CFL3 = x1x*x2x*rx3x1*rx3x2
     CLL = x3x*rx3x2
     CFR1 = x3x*x4x*rx3x2*rx4x2
     CFR2 =-x2x*x4x*rx3x2*rx4x3
     CFR3 = x2x*x3x*rx4x2*rx4x3
     CRR =-x2x*rx3x2
       cx0(n)=CFL1*CLL
       cx1(n)=CFL2*CLL+CFR1*CRR
       cx2(n)=CFL3*CLL+CFR2*CRR
       cx3(n)=CFR3*CRR
   enddo

   do m=1,mm
     j=jref(m)
     y1=yf(j)
     y2=yf(j+1)
     y3=yf(j+2)
     y4=yf(j+3)
     y = ya(m)
       y1y = y1-y   
       y2y = y2-y   
       y3y = y3-y   
       y4y = y4-y   
       ry2y1 = 1./(y2-y1)
       ry3y1 = 1./(y3-y1)
       ry4y1 = 1./(y4-y1)
       ry3y2 = 1./(y3-y2)
       ry4y2 = 1./(y4-y2)
       ry4y3 = 1./(y4-y3)
     CFL1 = y2y*y3y*ry2y1*ry3y1
     CFL2 =-y1y*y3y*ry2y1*ry3y2
     CFL3 = y1y*y2y*ry3y1*ry3y2
     CLL = y3y*ry3y2
     CFR1 = y3y*y4y*ry3y2*ry4y2
     CFR2 =-y2y*y4y*ry3y2*ry4y3
     CFR3 = y2y*y3y*ry4y2*ry4y3
     CRR =-y2y*ry3y2
       cy0(m)=CFL1*CLL
       cy1(m)=CFL2*CLL+CFR1*CRR
       cy2(m)=CFL3*CLL+CFR2*CRR
       cy3(m)=CFR3*CRR
   enddo

!
! Quadratic interpolations
!
   do n=1,nm
     i=irefq(n)
     x1=xf(i)
     x2=xf(i+1)
     x3=xf(i+2)
     x = xa(n)
       x1_x = x1-x   
       x2_x = x2-x   
       x3_x = x3-x   
       rx2x1 = 1./(x2-x1)
       rx3x1 = 1./(x3-x1)
       rx3x2 = 1./(x3-x2)
       qx0(n) = x2_x*x3_x*rx2x1*rx3x1
       qx1(n) =-x1_x*x3_x*rx2x1*rx3x2
       qx2(n) = x1_x*x2_x*rx3x1*rx3x2
   enddo

   do m=1,mm
     i=jrefq(m)
     y1=yf(i)
     y2=yf(i+1)
     y3=yf(i+2)
     y = ya(m)
       y1_y = y1-y   
       y2_y = y2-y   
       y3_y = y3-y   
       ry2y1 = 1./(y2-y1)
       ry3y1 = 1./(y3-y1)
       ry3y2 = 1./(y3-y2)
       qy0(m) = y2_y*y3_y*ry2y1*ry3y1
       qy1(m) =-y1_y*y3_y*ry2y1*ry3y2
       qy2(m) = y1_y*y2_y*ry3y1*ry3y2
   enddo
 
!-----------------------------------------------------------------------
                        endsubroutine lsqr_mg_coef

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_coef                    &
!***********************************************************************
!                                                                      !
!  Prepare coeficients for vertical mapping between:                   !
!  analysis grid vertical resolution (nm) and                          !
!  generation one of filter grid vertical resoluition (im)             !
!                                                                      !
!              (  im <= nm )                                           !
!                                                                      !
!***********************************************************************
(nm,im,c1,c2,c3,c4,iref)
use mg_mppstuff, only: mype
implicit none

integer(i_kind), intent(in):: nm,im
real(r_kind), dimension(1:nm), intent(out):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(out):: iref

real(r_kind), dimension(1:nm):: y
real(r_kind), dimension(0:im+1):: x
real(r_kind):: dy,x1,x2,x3,x4,dx1,dx2,dx3,dx4 
real(r_kind):: dx13,dx23,dx24

integer(i_kind):: i,n
!-----------------------------------------------------------------------

   do i=0,im+1
     x(i)=(i-1)*1.
   enddo

    dy = 1.*(im-1)/(nm-1)
  do n=1,nm
    y(n)=(n-1)*dy
  enddo
    y(nm)=x(im)
 
  do n=2,nm-1
    i = y(n)+1
      x1 = x(i-1)
      x2 = x(i)
      x3 = x(i+1)
      x4 = x(i+2)
    iref(n)=i
      dx1 = y(n)-x1
      dx2 = y(n)-x2
      dx3 = y(n)-x3
      dx4 = y(n)-x4
      dx13 = dx1*dx3
      dx23 = 0.5*dx2*dx3
      dx24 = dx2*dx4
    c1(n) = -dx23*dx3
    c2(n) =  (    dx13+0.5*dx24)*dx3
    c3(n) = -(0.5*dx13+    dx24)*dx2
    c4(n) = dx23*dx2

    if(iref(n)==1) then
      c3(n)=c3(n)+c1(n)
      c1(n)=0.
    endif
    if(iref(n)==im-1) then
      c2(n)=c2(n)+c4(n)
      c4(n)=0.
    endif
  enddo
     iref(1)=1; c1(1)=0.; c2(1)=1.; c3(1)=0.; c4(1)=0.
     iref(nm)=im; c1(nm)=0.; c2(nm)=1.; c3(nm)=0.; c4(n)=0.
     
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_coef                            

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_adjoint                 &
!***********************************************************************
!                                                                      !
!  Direct linerly weighted quadratic adjoint interpolation in vertical !
!  from reslution nm to resolution im                                  !
!                                                                      !
!              (  km <= nm )                                           !
!                                                                      !
!***********************************************************************
(nm,km,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,w,f)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(in):: w
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(out):: f
integer(i_kind):: k,n
!-----------------------------------------------------------------------
  f = 0.
do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    f(1,:,:) = f(1,:,:)+c2(n)*w(n,:,:)
    f(2,:,:) = f(2,:,:)+c3(n)*w(n,:,:)
    f(3,:,:) = f(3,:,:)+c4(n)*w(n,:,:)
  elseif &
    ( k==km-1) then
    f(km-2,:,:) = f(km-2,:,:)+c1(n)*w(n,:,:)
    f(km-1,:,:) = f(km-1,:,:)+c2(n)*w(n,:,:)
    f(km  ,:,:) = f(km  ,:,:)+c3(n)*w(n,:,:)
  elseif( k==km) then
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
  else
    f(k-1,:,:) = f(k-1,:,:)+c1(n)*w(n,:,:)
    f(k  ,:,:) = f(k  ,:,:)+c2(n)*w(n,:,:)
    f(k+1,:,:) = f(k+1,:,:)+c3(n)*w(n,:,:)
    f(k+2,:,:) = f(k+2,:,:)+c4(n)*w(n,:,:)
  endif
enddo
    f(1,:,:)=f(1,:,:)+w(1,:,:)
    f(km,:,:)=f(km,:,:)+w(nm,:,:)

!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_adjoint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_direct                  &
!***********************************************************************
!                                                                      !
!  Linerly weighted direct quadratic interpolation in vertical         !
!  from reslouion im to resolution nm                                  !
!                                                                      !
!              (  km <= nm )                                           !
!                                                                      !
!***********************************************************************
(km,nm,imin,imax,jmin,jmax,c1,c2,c3,c4,kref,f,w)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,nm,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm), intent(in):: c1,c2,c3,c4
integer(i_kind), dimension(1:nm), intent(in):: kref
real(r_kind), dimension(1:km,imin:imax,jmin:jmax), intent(in):: f
real(r_kind), dimension(1:nm,imin:imax,jmin:jmax), intent(out):: w
integer(i_kind):: k,n
!-----------------------------------------------------------------------
do n=2,nm-1
  k = kref(n)
  if( k==1 ) then
    w(n,:,:) =             c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  elseif &
    ( k==km-1) then
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,:)+c3(n)*f(k+1,:,:)
  elseif &
    ( k==km)   then
    w(n,:,:) =                 c2(n)*f(k,:,:)
  else
    w(n,:,:) =c1(n)*f(k-1,:,:)+c2(n)*f(k,:,: )+c3(n)*f(k+1,:,:)+c4(n)*f(k+2,:,:)
  endif
enddo
    w(1,:,:)=f(1,:,:)
    w(nm,:,:)=f(km,:,:)
    
 
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lwq_vertical_adjoint_spec            &
!***********************************************************************
!                                                                      !
!  Direct linerly weighted quadratic adjoint interpolation in vertical !
!  from reslution nm to resolution km                                  !
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
!$OMP PARALLEL
!$OMP PRIVATE(n,k)
  F = 0.
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

!$OMP END PARALLEL
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

!$OMP PARALLEL
!$OMP PRIVATE(n,k)

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
    
!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine lwq_vertical_direct_spec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine l_vertical_adjoint_spec              &
!***********************************************************************
!                                                                      !
!  Adjoint of linear interpolations in vertical                        !
!  from reslution nm to resolution km                                  !
!                                                                      !
!              (  nm = 2*km-1 )                                        !
!                                                                      !
!***********************************************************************
(km3,nm,km,imin,imax,jmin,jmax,W,F)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:nm), intent(in):: W
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:km), intent(out):: F
integer(i_kind):: k,n
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP PRIVATE(n,k)
  F = 0.

      k=1
  do n=2,nm-1,2
    F(:,:,:,k  ) = F(:,:,:,k  )+0.5*W(:,:,:,n)
    F(:,:,:,k+1) = F(:,:,:,k+1)+0.5*W(:,:,:,n)
      k=k+1
  enddo

      k=1
  do n=1,nm,2
    F(:,:,:,k  ) = F(:,:,:,k  )+    W(:,:,:,n)
      k=k+1
  enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine l_vertical_adjoint_spec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine l_vertical_direct_spec               &
!***********************************************************************
!                                                                      !
!                                                                      !
!  Direct linear interpolations in vertical                            !
!  from reslution nm to resolution km                                  !
!                                                                      !
!              (  nm = 2*km-1 )                                        !
!                                                                      !
!***********************************************************************
(km3,km,nm,imin,imax,jmin,jmax,F,W)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,km,nm,imin,imax,jmin,jmax
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:km), intent(in):: F
real(r_kind), dimension(1:km3,imin:imax,jmin:jmax,1:nm), intent(out):: W
integer(i_kind):: k,n
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP PRIVATE(n,k)

      k=1
  do n=1,nm,2
    W(:,:,:,n) =F (:,:,:,k)
      k=k+1
  enddo

      k=1
  do n=2,nm-1,2
    W(:,:,:,n) = 0.5*(F(:,:,:,k)+F(:,:,:,k+1))
      k=k+1
  enddo
    
!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine l_vertical_direct_spec

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_direct_offset                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(km,i0-ib:im+ib,j0-jb:jm+jb) perform          !
! direct interpolations to get target array W(km,1:nm,1:mm)            !
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(V,W,km,ibm,jbm)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: km,ibm,jbm
real(r_kind), dimension(km,1-ibm:im+ibm,1-jbm:jm+jbm), intent(in):: V
real(r_kind), dimension(km,1:nm,1:mm),intent(out):: W  

real(r_kind), dimension(km,1:nm,1-jb:jm+jbm):: VX
integer(i_kind):: i,j,n,m
real(r_kind),dimension(km):: v0,v1,v2,v3     
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP PRIVATE(j,n,i,v0,v1,v2,v3)

   do j=1-jbm,jm+jbm
   do n=1,nm
       i = iref(n)
     v0(:)=V(:,i  ,j)
     v1(:)=V(:,i+1,j)
     v2(:)=V(:,i+2,j)
     v3(:)=V(:,i+3,j)
     VX(:,n,j) = cx0(n)*v0(:)+cx1(n)*v1(:)+cx2(n)*v2(:)+cx3(n)*v3(:)
   enddo
   enddo

   do m=1,mm
     j = jref(m)
   do n=1,nm
     v0(:)=VX(:,n,j  ) 
     v1(:)=VX(:,n,j+1) 
     v2(:)=VX(:,n,j+2) 
     v3(:)=VX(:,n,j+3) 
     W(:,n,m) =  cy0(m)*v0(:)+cy1(m)*v1(:)+cy2(m)*v2(:)+cy3(m)*v3(:)
   enddo
   enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine lsqr_direct_offset

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine lsqr_adjoint_offset                  &
!***********************************************************************
!                                                                      !
! Given a target array W(km,1:nm,1:mm) perform adjoint                 !
! interpolations to get source array V(km,i0-ib:im+ib,j0-jb:jm+jb)     !
! using two passes of 1d interpolator                                  !
!                      - offset version -                              !
!                                                                      !
!***********************************************************************
(W,V,km,ibm,jbm)
!-----------------------------------------------------------------------
implicit none
integer(i_kind):: km,ibm,jbm
real(r_kind), dimension(km,1:nm,1:mm),intent(in):: W  
real(r_kind), dimension(km,1-ibm:im+ibm,1-jbm:jm+jbm), intent(out):: V
real(r_kind), dimension(km,1:nm,1-jbm:jm+jbm):: VX
real(r_kind), dimension(km):: wk
real(r_kind), dimension(km):: vxk
integer(i_kind):: i,j,n,m,l,k
integer(i_kind):: ip1,ip2,ip3
integer(i_kind):: jp1,jp2,jp3
real(r_kind):: c0,c1,c2,c3
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP PRIVATE(j,m,n,i,c0,c1,c2,c3,wk,vxk)

   V(:,:,:) = 0.

   VX(:,:,:)=0.

   do m=1,mm
     j = jref(m)
     c0 = cy0(m)
     c1 = cy1(m)
     c2 = cy2(m)
     c3 = cy3(m)
   do n=1,nm
       wk(:)=W(:,n,m)
     VX(:,n,j  ) = VX(:,n,j  )+wk(:)*c0
     VX(:,n,j+1) = VX(:,n,j+1)+wk(:)*c1
     VX(:,n,j+2) = VX(:,n,j+2)+wk(:)*c2
     VX(:,n,j+3) = VX(:,n,j+3)+wk(:)*c3
   enddo
   enddo
 

   do n=1,nm
     i = iref(n)
     c0 = cx0(n)
     c1 = cx1(n)
     c2 = cx2(n)
     c3 = cx3(n)
   do j=1-jbm,jm+jbm
       vxk(:)=VX(:,n,j)
     V(:,i  ,j) = V(:,i  ,j)+vxk(:)*c0
     V(:,i+1,j) = V(:,i+1,j)+vxk(:)*c1
     V(:,i+2,j) = V(:,i+2,j)+vxk(:)*c2
     V(:,i+3,j) = V(:,i+3,j)+vxk(:)*c3
   enddo
   enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine lsqr_adjoint_offset

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine quad_direct_offset                   &
!***********************************************************************
!                                                                      !
! Given a source array  V(km,i0-ib:im+ib,j0-jb:jm+jb) perform          !
! direct interpolations to get target array W(km,1:nm,1:mm)            !
! using two passes of 1d interpolator                                  !
!                                                                      !
!***********************************************************************
(V,W,km,ibm,jbm)
!-----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: km,ibm,jbm
real(r_kind), dimension(km,1-ibm:im+ibm,1-jbm:jm+jbm), intent(in):: V
real(r_kind), dimension(km,1:nm,1:mm),intent(out):: W  

real(r_kind), dimension(km,1:nm,1-jbm:jm+jbm):: VX
integer(i_kind):: i,j,n,m
real(r_kind),dimension(km):: v0,v1,v2
!-----------------------------------------------------------------------

!$OMP PARALLEL
!$OMP PRIVATE(i,j,n,m,v0,v1,v2)

   do j=1-jbm,jm+jbm
   do n=1,nm
       i = irefq(n)
     v0(:)=V(:,i  ,j)
     v1(:)=V(:,i+1,j)
     v2(:)=V(:,i+2,j)
     VX(:,n,j) = qx0(n)*v0(:)+qx1(n)*v1(:)+qx2(n)*v2(:)
   enddo
   enddo

   do m=1,mm
     j = jrefq(m)
   do n=1,nm
     v0(:)=VX(:,n,j  ) 
     v1(:)=VX(:,n,j+1) 
     v2(:)=VX(:,n,j+2) 
     W(:,n,m) =  qy0(m)*v0(:)+qy1(m)*v1(:)+qy2(m)*v2(:)
   enddo
   enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine quad_direct_offset

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine quad_adjoint_offset                  &
!***********************************************************************
!                                                                      !
! Given a target array W(km,1:nm,1:mm) perform adjoint                 !
! interpolations to get source array V(km,i0-ib:im+ib,j0-jb:jm+jb)     !
! using two passes of 1d interpolator                                  !
!                      - offset version -                              !
!                                                                      !
!***********************************************************************
(W,V,km,ibm,jbm)
!-----------------------------------------------------------------------
implicit none
integer(i_kind):: km,ibm,jbm
real(r_kind), dimension(km,1:nm,1:mm),intent(in):: W  
real(r_kind), dimension(km,1-ibm:im+ibm,1-jbm:jm+jbm), intent(out):: V
real(r_kind), dimension(km,1:nm,1-jbm:jm+jbm):: VX
real(r_kind), dimension(km):: wk
real(r_kind), dimension(km):: vxk
integer(i_kind):: i,j,n,m,l,k
real(r_kind):: c0,c1,c2
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP PRIVATE(j,m,n,i,c0,c1,c2,wk,vxk)

   V(:,:,:) = 0.

   VX(:,:,:)=0.

   do m=1,mm
     j = jrefq(m)
     c0 = qy0(m)
     c1 = qy1(m)
     c2 = qy2(m)
   do n=1,nm
       wk(:)=W(:,n,m)
     VX(:,n,j  ) = VX(:,n,j  )+wk(:)*c0
     VX(:,n,j+1) = VX(:,n,j+1)+wk(:)*c1
     VX(:,n,j+2) = VX(:,n,j+2)+wk(:)*c2
   enddo
   enddo
 

   do n=1,nm
     i = irefq(n)
     c0 = qx0(n)
     c1 = qx1(n)
     c2 = qx2(n)
   do j=1-jbm,jm+jbm
       vxk(:)=VX(:,n,j)
     V(:,i  ,j) = V(:,i  ,j)+vxk(:)*c0
     V(:,i+1,j) = V(:,i+1,j)+vxk(:)*c1
     V(:,i+2,j) = V(:,i+2,j)+vxk(:)*c2
   enddo
   enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine quad_adjoint_offset

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine l_vertical_adjoint_spec2             &
!***********************************************************************
!                                                                      !
!  Adjoint of linear interpolations in vertical                        !
!  from reslution nm to resolution km                                  !
!                                                                      !
!              (  nm = 2*km-1 )                                        !
!                                                                      !
!***********************************************************************
(en,nm,km,imin,imax,jmin,jmax,W,F)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: en,nm,km,imin,imax,jmin,jmax
real(r_kind), dimension(1:nm*en,imin:imax,jmin:jmax), intent(in):: W
real(r_kind), dimension(1:km*en,imin:imax,jmin:jmax), intent(out):: F
integer(i_kind):: k,n,e,enm,ekm
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP PRIVATE(n,k)
  F = 0.

do e=0,en-1
  enm = e*nm
  ekm = e*km
      k=1
  do n=2,nm-1,2
    F(ekm+k  ,:,:) = F(ekm+k  ,:,:)+0.5*W(enm+n,:,:)
    F(ekm+k+1,:,:) = F(ekm+k+1,:,:)+0.5*W(enm+n,:,:)
      k=k+1
  enddo

      k=1
  do n=1,nm,2
    F(ekm+k,:,:) = F(ekm+k,:,:) +  W(enm+n,:,:)
      k=k+1
  enddo
enddo

!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine l_vertical_adjoint_spec2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine l_vertical_direct_spec2              &
!***********************************************************************
!                                                                      !
!                                                                      !
!  Direct linear interpolations in vertical                            !
!  from reslution nm to resolution km                                  !
!                                                                     !
!              (  nmax = 2*kmax-1 )                                    !
!                                                                      !
!***********************************************************************
(en,km,nm,imin,imax,jmin,jmax,F,W)
implicit none
!-----------------------------------------------------------------------
integer(i_kind), intent(in):: en,km,nm,imin,imax,jmin,jmax
real(r_kind), dimension(1:km*en,imin:imax,jmin:jmax), intent(in):: F
real(r_kind), dimension(1:nm*en,imin:imax,jmin:jmax), intent(out):: W
integer(i_kind):: k,n,e,enm,ekm
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP PRIVATE(n,k)

do e=0,en-1
  enm = e*nm
  ekm = e*km
      k=1
  do n=1,nm,2
    W(enm+n,:,:) =F (ekm+k,:,:)
      k=k+1
  enddo
      k=1
  do n=2,nm-1,2
    W(enm+n,:,:) = 0.5*(F(ekm+k,:,:)+F(ekm+k+1,:,:))
      k=k+1
  enddo
enddo
    
!$OMP END PARALLEL
!-----------------------------------------------------------------------
                        endsubroutine l_vertical_direct_spec2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_interpolate
