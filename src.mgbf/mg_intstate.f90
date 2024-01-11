!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_intstate
!***********************************************************************
!                                                                      !
! Contains declarations and allocations of internal state variables    !
! use for filtering                                                    !
!                       - offset version -                             !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use jp_pkind2, only: fpi
!GSI use mpimod, only: mype
use mg_mppstuff, only: mype
use mg_parameter, only: im,jm,nh,hx,hy,pasp01,pasp02,pasp03
use mg_parameter, only: lm,hz,p,nm,mm,ib,jb,nb,mb
use mg_parameter, only: l_loc,lm_a
!GSI use berror, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4
use mg_parameter, only: mg_weig1,mg_weig2,mg_weig3,mg_weig4,mg_weig5
use mg_mppstuff, only: my_hgen,finishMPI,barrierMPI
use jp_pbfil,only: cholaspect
use jp_pbfil,only: getlinesum
use jp_pbfil3, only: inimomtab,t22_to_3,tritform,t33_to_6,hextform
!TEST
!use gridmod, only: lat1,lon1
!TEST
implicit none

!T real(r_kind), allocatable,dimension(:,:,:), pointer :: V
real(r_kind), dimension(:,:,:), pointer :: V
!
! Composite control variable on first generation o filter grid
!
real(r_kind), target,allocatable,dimension(:,:,:):: VALL
real(r_kind), allocatable,dimension(:,:,:):: HALL
!
! Composite control variable on high generations of filter grid
!
!
!FOR ADJOINT TEST
!
!real(r_kind), allocatable,dimension(:,:):: A
!real(r_kind), allocatable,dimension(:,:):: B
!real(r_kind), allocatable,dimension(:,:):: A0
!real(r_kind), allocatable,dimension(:,:):: B0
!
real(r_kind), allocatable,dimension(:,:,:):: a_diff_f
real(r_kind), allocatable,dimension(:,:,:):: a_diff_h
real(r_kind), allocatable,dimension(:,:,:):: b_diff_f
real(r_kind), allocatable,dimension(:,:,:):: b_diff_h

!
! Localization weights
!
real(r_kind), allocatable,dimension(:,:,:):: w1_loc
real(r_kind), allocatable,dimension(:,:,:):: w2_loc
real(r_kind), allocatable,dimension(:,:,:):: w3_loc   
real(r_kind), allocatable,dimension(:,:,:):: w4_loc   

real(r_kind), allocatable,dimension(:,:):: p_eps
real(r_kind), allocatable,dimension(:,:):: p_del
real(r_kind), allocatable,dimension(:,:):: p_sig
real(r_kind), allocatable,dimension(:,:):: p_rho

real(r_kind), allocatable,dimension(:,:,:):: paspx
real(r_kind), allocatable,dimension(:,:,:):: paspy
real(r_kind), allocatable,dimension(:,:,:):: pasp1
real(r_kind), allocatable,dimension(:,:,:,:):: pasp2
real(r_kind), allocatable,dimension(:,:,:,:,:):: pasp3

real(r_kind), allocatable,dimension(:,:,:):: vpasp2
real(r_kind), allocatable,dimension(:,:,:):: hss2
real(r_kind), allocatable,dimension(:,:,:,:):: vpasp3
real(r_kind), allocatable,dimension(:,:,:,:):: hss3

real(r_kind), allocatable,dimension(:):: ssx
real(r_kind), allocatable,dimension(:):: ssy
real(r_kind), allocatable,dimension(:):: ss1
real(r_kind), allocatable,dimension(:,:):: ss2
real(r_kind), allocatable,dimension(:,:,:):: ss3

integer(fpi), allocatable,dimension(:,:,:):: dixs
integer(fpi), allocatable,dimension(:,:,:):: diys
integer(fpi), allocatable,dimension(:,:,:):: dizs

integer(fpi), allocatable,dimension(:,:,:,:):: dixs3
integer(fpi), allocatable,dimension(:,:,:,:):: diys3
integer(fpi), allocatable,dimension(:,:,:,:):: dizs3

integer(fpi), allocatable,dimension(:,:,:,:):: qcols

!real(r_kind), allocatable,dimension(:,:,:,:):: r_vol
!
!
! Composite stacked variable
!

real(r_kind), target,allocatable,dimension(:,:,:):: WORKA
real(r_kind), allocatable,dimension(:,:,:):: WORK 


integer(i_kind),allocatable,dimension(:):: iref,jref
integer(i_kind),allocatable,dimension(:):: irefq,jrefq
integer(i_kind),allocatable,dimension(:):: irefL,jrefL

integer(i_kind),allocatable,dimension(:):: Lref,Lref_h
real(r_kind),allocatable,dimension(:):: cvf1,cvf2,cvf3,cvf4
real(r_kind),allocatable,dimension(:):: cvh1,cvh2,cvh3,cvh4

real(r_kind),allocatable,dimension(:):: cx0,cx1,cx2,cx3
real(r_kind),allocatable,dimension(:):: cy0,cy1,cy2,cy3

real(r_kind),allocatable,dimension(:):: qx0,qx1,qx2
real(r_kind),allocatable,dimension(:):: qy0,qy1,qy2

real(r_kind),allocatable,dimension(:):: Lx0,Lx1
real(r_kind),allocatable,dimension(:):: Ly0,Ly1

real(r_kind),allocatable,dimension(:):: p_coef,q_coef
real(r_kind),allocatable,dimension(:):: a_coef,b_coef

real(r_kind),allocatable,dimension(:,:):: cf00,cf01,cf02,cf03           &
                                         ,cf10,cf11,cf12,cf13           &
                                         ,cf20,cf21,cf22,cf23           &
                                         ,cf30,cf31,cf32,cf33

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine allocate_mg_intstate(km_all,km_a_all)
!***********************************************************************
!                                                                      !
! Allocate internal state variables                                    !
!                                                                      !
!***********************************************************************
implicit none

integer(i_kind),intent(in):: km_all,km_a_all

integer(i_kind) :: km_all_4,km_all_16,km_all_64

if(l_loc) then
   km_all_4  = km_all/4
   km_all_16 = km_all/16
   km_all_64 = km_all/64
     allocate(w1_loc(km_all   ,1-hx:im+hx,1-hy:jm+hy))         ; w1_loc=0. 
     allocate(w2_loc(km_all_4 ,1-hx:im+hx,1-hy:jm+hy))         ; w2_loc=0. 
     allocate(w3_loc(km_all_16,1-hx:im+hx,1-hy:jm+hy))         ; w3_loc=0. 
     allocate(w4_loc(km_all_64,1-hx:im+hx,1-hy:jm+hy))         ; w4_loc=0. 
endif


allocate(V(1-hx:im+hx,1-hy:jm+hy,lm))                      ; V=0.
allocate(VALL(km_all,1-hx:im+hx,1-hy:jm+hy))                  ; VALL=0.
allocate(HALL(km_all,1-hx:im+hx,1-hy:jm+hy))                  ; HALL=0.


allocate(a_diff_f(km_all,1-hx:im+hx,1-hy:jm+hy))             ; a_diff_f=0. 
allocate(a_diff_h(km_all,1-hx:im+hx,1-hy:jm+hy))             ; a_diff_h=0. 
allocate(b_diff_f(km_all,1-hx:im+hx,1-hy:jm+hy))             ; b_diff_f=0. 
allocate(b_diff_h(km_all,1-hx:im+hx,1-hy:jm+hy))             ; b_diff_h=0. 

allocate(p_eps(1-hx:im+hx,1-hy:jm+hy))                            ; p_eps=0.
allocate(p_del(1-hx:im+hx,1-hy:jm+hy))                            ; p_del=0.
allocate(p_sig(1-hx:im+hx,1-hy:jm+hy))                            ; p_sig=0.
allocate(p_rho(1-hx:im+hx,1-hy:jm+hy))                            ; p_rho=0.

allocate(paspx(1,1,1:im))                                          ; paspx=0.
allocate(paspy(1,1,1:jm))                                          ; paspy=0.

allocate(pasp1(1,1,1:lm))                                           ; pasp1=0.
allocate(pasp2(2,2,1:im,1:jm))                                    ; pasp2=0.
allocate(pasp3(3,3,1:im,1:jm,1:lm))                               ; pasp3=0.

allocate(vpasp2(0:2,1:im,1:jm))                                   ; vpasp2=0.
allocate(hss2(1:im,1:jm,1:3))                                     ; hss2= 0.

allocate(vpasp3(1:6,1:im,1:jm,1:lm))                              ; vpasp3= 0.
allocate(hss3(1:im,1:jm,1:lm,1:6))                                ; hss3= 0.

allocate(ssx(1:im))                                             ; ssx=0.
allocate(ssy(1:jm))                                             ; ssy=0.
allocate(ss1(1:lm))                                             ; ss1=0.
allocate(ss2(1:im,1:jm))                                        ; ss2=0.
allocate(ss3(1:im,1:jm,1:lm))                                   ; ss3=0.

allocate(dixs(1:im,1:jm,3))                                     ; dixs=0
allocate(diys(1:im,1:jm,3))                                     ; diys=0

allocate(dixs3(1:im,1:jm,1:lm,6))                               ; dixs3=0
allocate(diys3(1:im,1:jm,1:lm,6))                               ; diys3=0
allocate(dizs3(1:im,1:jm,1:lm,6))                               ; dizs3=0

allocate(qcols(0:7,1:im,1:jm,1:lm))                             ; qcols=0

!
! In stnadalone version
!
!allocate(r_vol(km,0:nm,0:mm,2))                             ; r_vol=0.
!
! ... but in global version there will be 
!     r_vol2 and r_vol3 for 2d and 3d variables
! and r_vol3 will need to be given vertical structure
!

!
allocate(WORKA(km_a_all,1:nm,1:mm))                        ; WORKA=0.
allocate(WORK (km_all  ,1:nm,1:mm))                        ; WORK =0.

!
! for re-decomposition
!

allocate(iref(1:nm))                                     ; iref=0
allocate(jref(1:mm))                                     ; jref=0

allocate(irefq(1:nm))                                    ; irefq=0
allocate(jrefq(1:mm))                                    ; jrefq=0

allocate(irefL(1:nm))                                    ; irefq=0
allocate(jrefL(1:mm))                                    ; jrefq=0

allocate(cx0(1:nm))                                      ; cx0=0.
allocate(cx1(1:nm))                                      ; cx1=0.
allocate(cx2(1:nm))                                      ; cx2=0.
allocate(cx3(1:nm))                                      ; cx3=0.

allocate(cy0(1:mm))                                      ; cy0=0.
allocate(cy1(1:mm))                                      ; cy1=0.
allocate(cy2(1:mm))                                      ; cy2=0.
allocate(cy3(1:mm))                                      ; cy3=0.

allocate(qx0(1:nm))                                      ; qx0=0.
allocate(qx1(1:nm))                                      ; qx1=0.
allocate(qx2(1:nm))                                      ; qx2=0.

allocate(qy0(1:mm))                                      ; qy0=0.
allocate(qy1(1:mm))                                      ; cy1=0.
allocate(qy2(1:mm))                                      ; qy2=0.

allocate(Lx0(1:nm))                                      ; Lx0=0.
allocate(Lx1(1:nm))                                      ; Lx1=0.

allocate(Ly0(1:nm))                                      ; Ly0=0.
allocate(Ly1(1:nm))                                      ; Ly1=0.

!TEST
!       call finishMPI
!TEST

allocate(p_coef(4))                                      ; p_coef=0.
allocate(q_coef(4))                                      ; q_coef=0.

allocate(a_coef(3))                                      ; a_coef=0.
allocate(b_coef(3))                                      ; b_coef=0.


allocate(cf00(1:nm,1:mm))                            ; cf00=0.
allocate(cf01(1:nm,1:mm))                            ; cf01=0.
allocate(cf02(1:nm,1:mm))                            ; cf02=0.
allocate(cf03(1:nm,1:mm))                            ; cf03=0.
allocate(cf10(1:nm,1:mm))                            ; cf10=0.
allocate(cf11(1:nm,1:mm))                            ; cf11=0.
allocate(cf12(1:nm,1:mm))                            ; cf12=0.
allocate(cf13(1:nm,1:mm))                            ; cf13=0.
allocate(cf20(1:nm,1:mm))                            ; cf20=0.
allocate(cf21(1:nm,1:mm))                            ; cf21=0.
allocate(cf22(1:nm,1:mm))                            ; cf22=0.
allocate(cf23(1:nm,1:mm))                            ; cf23=0.
allocate(cf30(1:nm,1:mm))                            ; cf30=0.
allocate(cf31(1:nm,1:mm))                            ; cf31=0.
allocate(cf32(1:nm,1:mm))                            ; cf32=0.
allocate(cf33(1:nm,1:mm))                            ; cf33=0.

allocate(Lref(1:lm_a))                                  ; Lref=0
allocate(Lref_h(1:lm))                                 ; Lref_h=0

allocate(cvf1(1:lm_a))                                 ; cvf1=0.
allocate(cvf2(1:lm_a))                                 ; cvf2=0.
allocate(cvf3(1:lm_a))                                 ; cvf3=0.
allocate(cvf4(1:lm_a))                                 ; cvf4=0.

allocate(cvh1(1:lm))                                   ; cvh1=0.
allocate(cvh2(1:lm))                                   ; cvh2=0.
allocate(cvh3(1:lm))                                   ; cvh3=0.
allocate(cvh4(1:lm))                                   ; cvh4=0.


!-----------------------------------------------------------------------
                        endsubroutine allocate_mg_intstate

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine def_mg_weights
!***********************************************************************
!                                                                      !
! Define weights and scales                                            !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L
real(r_kind):: gen_fac
!-----------------------------------------------------------------------
if(mype.eq.0) then
  open(unit=10,file='a01.dat',status='new',action='write')
    write(10,'(f5.2)') pasp01
  close(10)

  open(unit=10,file='a02.dat',status='new',action='write')
    write(10,'(f5.2)') pasp02
  close(10)
endif
!-----------------------------------------------------------------------

      p_eps(:,:)=0.0
      p_del(:,:)=0.0
      p_sig(:,:)=0.0
      p_rho(:,:)=0.0

!--------------------------------------------------------
!
! For localization (for now)
!
    if(l_loc) then
      w1_loc(:,:,:)=mg_weig1
      w2_loc(:,:,:)=mg_weig2
      w3_loc(:,:,:)=mg_weig3
      w4_loc(:,:,:)=mg_weig4
    endif 
!--------------------------------------------------------
      gen_fac=1.
      a_diff_f(:,:,:)=mg_weig1 
      a_diff_h(:,:,:)=mg_weig1 

      b_diff_f(:,:,:)=0.
      b_diff_h(:,:,:)=0.

!      r_vol(:,:,:,1)=1.


      select case(my_hgen)
        case(2) 
          a_diff_h(:,:,:)=mg_weig2
        case(3) 
          a_diff_h(:,:,:)=mg_weig3 
        case(4) 
          a_diff_h(:,:,:)=mg_weig4
        case default 
          a_diff_h(:,:,:)=mg_weig5
      end select

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'In def mg_weights: After a_deff_h'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


          do L=1,lm
           pasp1(1,1,L)=pasp01
          enddo

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'In def mg_weights: After pasp1'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

          do i=1,im
            paspx(1,1,i)=pasp02
          enddo  
          do j=1,jm
            paspy(1,1,j)=pasp02
          enddo  

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'In def mg_weights: After paspx and paspy'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

          do j=1,jm
          do i=1,im
            pasp2(1,1,i,j)=pasp02*(1.+p_del(i,j))
            pasp2(2,2,i,j)=pasp02*(1.-p_del(i,j))
            pasp2(1,2,i,j)=pasp02*p_eps(i,j)     
            pasp2(2,1,i,j)=pasp02*p_eps(i,j)     
          end do
          end do

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'In def mg_weights: After pasp2 '
!      endif
!     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        do L=1,lm
          do j=1,jm
          do i=1,im
            pasp3(1,1,i,j,l)=pasp03*(1+p_del(i,j))
            pasp3(2,2,i,j,l)=pasp03
            pasp3(3,3,i,j,l)=pasp03*(1-p_del(i,j))
            pasp3(1,2,i,j,l)=pasp03*p_eps(i,j)
            pasp3(2,1,i,j,l)=pasp03*p_eps(i,j)
            pasp3(2,3,i,j,l)=pasp03*p_sig(i,j)
            pasp3(3,2,i,j,l)=pasp03*p_sig(i,j)
            pasp3(1,3,i,j,l)=pasp03*p_rho(i,j)
            pasp3(3,1,i,j,l)=pasp03*p_rho(i,j)
          end do
          end do
        end do
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'In def mg_weights: After pasp3 '
!      endif
!     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


        call cholaspect(1,lm,pasp1)
        call cholaspect(1,im,1,jm,pasp2)
        call cholaspect(1,im,1,jm,1,lm,pasp3)


        call getlinesum(hx,1,im,paspx,ssx)
        call getlinesum(hy,1,jm,paspy,ssy)
        call getlinesum(hz,1,lm,pasp1,ss1)
        call getlinesum(hx,1,im,hy,1,jm,pasp2,ss2)
        call getlinesum(hx,1,im,hy,1,jm,hz,1,lm,pasp3,ss3)
!-----------------------------------------------------------------------
                        endsubroutine def_mg_weights

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine init_mg_line
!***********************************************************************
!                                                                      !
! Inititate line filters                                               !
!                                                                      !
!***********************************************************************
integer(i_kind):: i,j,L,icol
logical:: ff
!-----------------------------------------------------------------------

  do j=1,jm
  do i=1,im
    call t22_to_3(pasp2(:,:,i,j),vpasp2(:,i,j))
  enddo
  enddo

  do l=1,lm
  do j=1,jm
  do i=1,im
    call t33_to_6(pasp3(:,:,i,j,l),vpasp3(:,i,j,l))
  enddo
  enddo
  enddo



  call inimomtab(p,nh,ff)

  call tritform(1,im,1,jm,vpasp2, dixs,diys, ff)

  do icol=1,3
    hss2(:,:,icol)=vpasp2(icol-1,:,:)
  enddo  


  call hextform(1,im,1,jm,1,lm,vpasp3,qcols,dixs3,diys3,dizs3, ff)


  do icol=1,6
    hss3(:,:,:,icol)=vpasp3(icol,:,:,:)
  enddo
 

!-----------------------------------------------------------------------
                        endsubroutine init_mg_line

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine deallocate_mg_intstate
!***********************************************************************
!                                                                      !
! Deallocate internal state variables                                  !
!                                                                      !
!***********************************************************************
implicit none

deallocate(V)

deallocate(HALL,VALL)

deallocate(a_diff_f,b_diff_f)
deallocate(a_diff_h,b_diff_h)
deallocate(p_eps,p_del,p_sig,p_rho,pasp1,pasp2,pasp3,ss1,ss2,ss3)
deallocate(dixs,diys)
deallocate(dixs3,diys3,dizs3)
deallocate(qcols)
!
! for testing
!
deallocate(WORKA)
deallocate(WORK )

!
! for re-decomposition
!
deallocate(iref,jref)
deallocate(irefq,jrefq)
deallocate(irefL,jrefL)

deallocate(cf00,cf01,cf02,cf03,cf10,cf11,cf12,cf13)
deallocate(cf20,cf21,cf22,cf23,cf30,cf31,cf32,cf33)

deallocate(Lref,Lref_h)

deallocate(cvf1,cvf2,cvf3,cvf4)

deallocate(cvh1,cvh2,cvh3,cvh4)

deallocate(cx0,cx1,cx2,cx3)
deallocate(cy0,cy1,cy2,cy3)

deallocate(qx0,qx1,qx2)
deallocate(qy0,qy1,qy2)

deallocate(Lx0,Lx1)
deallocate(Ly0,Ly1)

deallocate(p_coef,q_coef)
deallocate(a_coef,b_coef)
  
if(l_loc) then
  deallocate(w1_loc,w2_loc,w3_loc,w4_loc)
endif

!----------------------------------------------------------------------
                        endsubroutine deallocate_mg_intstate

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_intstate
