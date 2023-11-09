!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_entrymod
!***********************************************************************
!                                                                      !
!   Initialize and finialize multigrid Beta filter for modeling of     !
!   background error covariance                                        !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_parameter
!use mpimod, only: mype
use mg_mppstuff, only: mype
use mg_mppstuff, only: init_mg_MPI,finishMPI,barrierMPI
use mg_domain, only: init_mg_domain
use mg_domain_loc, only: init_domain_loc
use mg_intstate, only: allocate_mg_intstate,def_mg_weights              & 
                      ,init_mg_line                                     &
                      ,deallocate_mg_intstate                           &
                      ,cvf1,cvf2,cvf3,cvf4,lref                         &
                      ,WORKA
use mg_interpolate,only: lsqr_mg_coef,lwq_vertical_coef,def_offset_coef
use mg_input, only: input_2d,input_3d,input_spec1_2d
use mg_output, only: output_spec1_2d,output_vertical_2d

implicit none

public mg_initialize
public mg_finalize

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_initialize
!**********************************************************************!
!                                                                      !
!   Initialization subroutine                                          !
!                                                     M. Rancic (2020) !
!***********************************************************************

real(r_kind), allocatable, dimension(:,:):: PA

!---------------------------------------------------------------------------
!
!               Firs set of subroutines is called only once and serves to 
!               initialte the MGBF run                                 
! 
!---------------------------------------------------------------------------

!****
!**** Initialize run multigrid Beta filter parameters
!****

      call init_mg_parameter


!****
!**** Initialize MPI
!****

      call init_mg_MPI

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish init_mg_MPI'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!***
!*** Initialize integration domain
!***

      call init_mg_domain

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish init_mg_domain'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      if(l_loc) then
        call init_domain_loc
      endif

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish init_mg_loc'
!      endif

!      call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


!---------------------------------------------------------------------------
!
!               All others are function of km2,km3,km,nm,mm,im,jm
!               and needs to be called separately for each application
! 
!---------------------------------------------------------------------------
!***
!*** Define km and WORKA array based on input from mg_parameters and
!*** depending on specific application
!***
 

!***
!*** Allocate variables, define weights, prepare mapping 
!*** between analysis and filter grid
!***

      call allocate_mg_intstate(km_all,km_a_all)

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish allocate_mg_intstate'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call def_offset_coef

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish def_offset_coef'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      call def_mg_weights

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!        print *,'Finish def_mg_weights'
!      endif
!      call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      if(mgbf_line) then
         call init_mg_line
      endif

      call lsqr_mg_coef 

      call lwq_vertical_coef(lm_a,lm,cvf1,cvf2,cvf3,cvf4,lref)

!***
!*** Just for testing of standalone version. In GSI WORKA will be given
!*** through a separate subroutine 
!***

!    call input_3d(WORKA(     1:  lm,:,:),1,1,     1,mm,nm,  lm,mm0,4,3)
!    call input_3d(WORKA(  lm+1:2*lm,:,:),1,1,  lm+1,mm,nm,2*lm,mm0,6,5)
!    call input_3d(WORKA(2*lm+1:3*lm,:,:),1,1,2*lm+1,mm,nm,3*lm,mm0,2,1)
!    call input_3d(WORKA(3*lm+1:4*lm,:,:),1,1,3*lm+1,mm,nm,4*lm,mm0,3,2)
!    call input_3d(WORKA(4*lm+1:5*lm,:,:),1,1,4*lm+1,mm,nm,5*lm,mm0,7,3)
!    call input_3d(WORKA(5*lm+1:6*lm,:,:),1,1,5*lm+1,mm,nm,6*lm,mm0,4,5)

!    call input_3d(WORKA(6*lm+1:6*lm+1,:,:),1,1,6*lm+1,mm,nm,6*lm+1,mm0,2,1)
!    call input_3d(WORKA(6*lm+2:6*lm+2,:,:),1,1,6*lm+2,mm,nm,6*lm+2,mm0,4,1)
!    call input_3d(WORKA(6*lm+3:6*lm+3,:,:),1,1,6*lm+3,mm,nm,6*lm+3,mm0,5,1)
!    call input_3d(WORKA(6*lm+4:6*lm+4,:,:),1,1,6*lm+4,mm,nm,6*lm+4,mm0,7,1)

     WORKA(:,:,:)=0.

if(ldelta) then

allocate(PA(1:nm,1:mm))

    PA = 0.
    call input_spec1_2d(PA,nxm/2,nym/2,'md')

!    WORKA(3*lm+1:4*lm,:,:)=0.
    WORKA(3*lm_a+lm_a/2,:,:)=PA(:,:)


deallocate(PA)

endif
!TEST
!             call finishMPI
!TEST

!-----------------------------------------------------------------------
                        endsubroutine mg_initialize

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine mg_finalize
!**********************************************************************!
!                                                                      !
!   Finalize multigrid Beta Function                                   !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mg_parameter, only: nm,mm

real(r_kind), allocatable, dimension(:,:):: PA, VA
integer(i_kind):: n,m,L
!-----------------------------------------------------------------------

if(ldelta) then

!
! Horizontal cross-section
!

allocate(PA(1:nm,1:mm))

     PA(:,:)=WORKA(3*lm_a+lm_a/2,:,:)

          call output_spec1_2d(PA)

deallocate(PA)

!
! Vertical cross-section
!
     
allocate(VA(1:nm,1:lm_a))


     do l=1,lm_a
     do n=1,nm
        VA(n,l)=WORKA(3*lm_a+l,n,mm/2)
     enddo
     enddo

          call output_vertical_2d(VA,nym/2)

deallocate(VA)

endif

     call barrierMPI


          call deallocate_mg_intstate          

!-----------------------------------------------------------------------
                        endsubroutine mg_finalize
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_entrymod
