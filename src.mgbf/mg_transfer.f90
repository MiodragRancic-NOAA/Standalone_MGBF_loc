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
!TEST
!use mg_output, only: output_spec1_2dd
!TEST

implicit none
integer(i_kind):: n,m,l,k,i,j

public anal_to_filt_all
public filt_to_anal_all

public anal_to_filt
public filt_to_anal

public stack_to_composite
public composite_to_stack

public stack_to_composite_spec
public composite_to_stack_spec

public stack_to_composite_loc
public composite_to_stack_loc

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

allocate(A3D(km3,1:nm,1:mm,lm_a))
allocate(A2D(km2,1:nm,1:mm))

allocate(F3D(km3,1:nm,1:mm,lm))


!                                                 call btim(  s2c_tim)
     call stack_to_composite_spec(WORKA,A2D,A3D,lm_a,km_a)
!                                                 call etim(  s2c_tim)

!                                                 call btim(  vadj_tim)
  if(lm_a>lm) then
     call lwq_vertical_adjoint_spec(km3,lm_a,lm,1,nm,1,mm,               &
                                    cvf1,cvf2,cvf3,cvf4,lref,A3D,F3D)
  else 

    do L=1,lm
      F3D(:,:,:,L)=A3D(:,:,:,L)
    enddo

  endif
!                                                 call etim(  vadj_tim)

!                                                 call btim(  c2s_tim)
     call composite_to_stack_spec(A2D,F3D,WORK,lm,km)
!                                                 call etim(  c2s_tim)

!                                                 call btim(  a2f_tim)
     call anal_to_filt
!                                                 call etim(  a2f_tim)


deallocate(A3D,A2D,F3D)
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
real(r_kind),allocatable,dimension(:,:,:):: A2D

real(r_kind),allocatable,dimension(:,:,:,:):: F3D

!----------------------------------------------------------------------

allocate(A3D(km3,1:nm,1:mm,lm_a))
allocate(A2D(km2,1:nm,1:mm))

allocate(F3D(km3,1:nm,1:mm,lm))


    
    call filt_to_anal
   
    call stack_to_composite_spec(WORK,A2D,F3D,lm,km)
  
 if(lm_a>lm) then
    call lwq_vertical_direct_spec(km3,lm,lm_a,1,nm,1,mm,                 &
                                  cvf1,cvf2,cvf3,cvf4,lref,F3D,A3D)
 else

   do L=1,lm
     A3D(:,:,:,L)=F3D(:,:,:,L)
   enddo
 
 endif  
 

    call composite_to_stack_spec(A2D,A3D,WORKA,lm_a,km_a)

deallocate(A3D,A2D,F3D)
!----------------------------------------------------------------------
                        endsubroutine filt_to_anal_all

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine anal_to_filt
!***********************************************************************
!                                                                      !
!  Transfer data from analysis to first generaton of filter grid       !
!                                                                      !
!***********************************************************************
use mg_interpolate, only: lsqr_adjoint_offset
use mg_bocos, only:  bocoT_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: VLOC  

!----------------------------------------------------------------------

    allocate(VLOC(km,1-ib:im+ib,1-jb:jm+jb))                      


!T                                                 call btim(  aintp_tim)

      VLOC=0.
         call lsqr_adjoint_offset(WORK,VLOC,km)


!T                                                 call etim(  aintp_tim)


!***
!***  Apply adjoint lateral bc on PKF and WKF
!***
    

                                                 call btim(  bocoT_tim)
         call bocoT_2d(VLOC,km,im,jm,ib,jb)
                                                 call etim(  bocoT_tim)
 
       VALL=0.
       VALL(1:km,1:im,1:jm)=VLOC(1:km,1:im,1:jm)
      

    deallocate(VLOC)

!                                            call etim(   btrns1_tim)

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
use mg_bocos, only:  boco_2d
implicit none


real(r_kind),allocatable,dimension(:,:,:):: VLOC   
!TEST
!real(r_kind), allocatable, dimension(:,:):: PA
!TEST

!----------------------------------------------------------------------

!T                                            call btim(   btrns2_tim)

!***
!***  Define VLOC
!***

    allocate(VLOC(1:km,1-ib:im+ib,1-jb:jm+jb))                     

      VLOC=0.
      VLOC(1:km,1:im,1:jm)=VALL(1:km,1:im,1:jm)
        

!***
!***  Supply boundary conditions for VLOC
!***
                                                 call btim(  boco_tim)
         call boco_2d(VLOC,km,im,jm,ib,jb)

                                                 call etim(  boco_tim)


!T                                                 call btim(   intp_tim)

         call lsqr_direct_offset(VLOC,WORK,km)

!T                                                 call etim(   intp_tim)
    deallocate(VLOC)


!T                                                 call etim(   btrns2_tim)

!----------------------------------------------------------------------
                        endsubroutine filt_to_anal


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
                        subroutine stack_to_composite_loc               &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(in):: ARR_ALL
real(r_kind),dimension(km3*n_ens,1-hx:im+hx,1-hy:jm+hy,lm),intent(out):: A3D
real(r_kind),dimension(km2*n_ens,1-hx:im+hx,1-hy:jm+hy)   ,intent(out):: A2D
integer(i_kind):: k,L
integer(i_kind):: n,n_inc
!----------------------------------------------------------------------
  do n=1,n_ens
    n_inc = km*(n-1)

    do k=1,km3
      do L=1,lm
        A3D(n_inc+k,:,:,L)=ARR_ALL(n_inc+(k-1)*lm+L,:,:)
      enddo
    enddo

      do k=1,km2
        A2D(n_inc+k,:,:)  =ARR_ALL(n_inc+km3*lm+k,:,:)
      enddo 

   enddo
 
!----------------------------------------------------------------------
                        endsubroutine stack_to_composite_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine composite_to_stack_loc               &
!***********************************************************************
!                                                                      !
!  Transfer data from composite to stack variables                     !
!                                                                      !
!***********************************************************************
(A2D,A3D,ARR_ALL)
!----------------------------------------------------------------------
implicit none
real(r_kind),dimension(km2*n_ens,1-hx:im+hx,1-hy:jm+hy),   intent(in):: A2D
real(r_kind),dimension(km3*n_ens,1-hx:im+hx,1-hy:jm+hy,lm),intent(in):: A3D
real(r_kind),dimension(km ,1-hx:im+hx,1-hy:jm+hy),   intent(out):: ARR_ALL
integer(i_kind):: k,L
integer(i_kind):: n,n_inc
!----------------------------------------------------------------------

  do n=1,n_ens
    n_inc = km*(n-1)

    do k=1,km3
      do L=1,lm
        ARR_ALL(n_inc+(k-1)*lm+L,:,: )=A3D(n_inc+k,:,:,L)
      enddo
    enddo

    do k=1,km2
      ARR_ALL(n_inc+km3*lm+k,:,:)=A2D(n_inc+k,:,:)
    enddo 

  enddo
!----------------------------------------------------------------------
                        endsubroutine composite_to_stack_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine stack_to_composite_spec              &
!***********************************************************************
!                                                                      !
!  Transfer data from stack to composite variables                     !
!                                                                      !
!***********************************************************************
(ARR_ALL,A2D,A3D,lmax,kmax)
!----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: lmax,kmax
real(r_kind),dimension(kmax,1:nm,1:mm)     ,intent(in):: ARR_ALL
real(r_kind),dimension(km3 ,1:nm,1:mm,lmax),intent(out):: A3D
real(r_kind),dimension(km2 ,1:nm,1:mm)     ,intent(out):: A2D
integer(i_kind):: n,m,L,k
!----------------------------------------------------------------------
    do L=1,lmax
      do m=1,mm
      do n=1,nm
        do k=1,km3
          A3D(k,n,m,L)=ARR_ALL( (k-1)*lmax+L,n,m )
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          A2D(k,:,:)=ARR_ALL(km3*lmax+k,:,:)
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
(A2D,A3D,ARR_ALL,lmax,kmax)
!----------------------------------------------------------------------
implicit none
integer(i_kind),intent(in):: lmax,kmax
real(r_kind),dimension(km2 ,1:nm,1:mm),     intent(in):: A2D
real(r_kind),dimension(km3 ,1:nm,1:mm,lmax),intent(in):: A3D
real(r_kind),dimension(kmax,1:nm,1:mm),     intent(out):: ARR_ALL
integer(i_kind):: n,m,L,k
!----------------------------------------------------------------------
    do L=1,lmax
      do m=1,mm
      do n=1,nm
        do k=1,km3
          ARR_ALL( (k-1)*lmax+L,n,m )=A3D(k,n,m,L)
        enddo
      enddo
      enddo
    enddo

        do k=1,km2
          ARR_ALL(km3*lmax+k,:,:)=A2D(k,:,:)
        enddo 

!----------------------------------------------------------------------
                        endsubroutine composite_to_stack_spec


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_transfer
