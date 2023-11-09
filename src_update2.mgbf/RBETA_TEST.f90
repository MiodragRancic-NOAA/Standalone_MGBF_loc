!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        program RBETA_TEST 
!***********************************************************************
!                                                                      !
!   Multigrid Beta filter for modeling background error covariance     !
!                                                                      !
!                                                     M. Rancic (2020) !
!***********************************************************************
use mpi
use kinds, only: r_kind,i_kind
use mg_entrymod, only: mg_initialize,mg_finalize
use mg_mppstuff, only: finishMPI,mype
use mg_filtering, only: filtering_procedure
use mg_transfer, only: anal_to_filt_all,filt_to_anal_all 
use mg_transfer, only: anal_to_filt_all2,filt_to_anal_all2 
use mg_parameter, only: mgbf_proc,l_new_map
use mg_timers

implicit none

!-----------------------------------------------------------------------

                                                   call btim(   total_tim)
                                                   call btim(    init_tim)
!***
!*** Initialzie multigrid Beta filter                                   
!***

          call mg_initialize
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'finish mg_initialize'
!      endif
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                                                   call etim(    init_tim)
!***
!*** From the analysis to first generation of filter grid
!***
                                                   call btim(    an2filt_tim)

        if(l_new_map) then
          call anal_to_filt_all2
        else
          call anal_to_filt_all
        endif 
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'finish anal_to_filter_all'  
!     endif
!     call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

                                                   call etim(    an2filt_tim)


!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!***
!*** Filtering
!***
!======================================================================

       call filtering_procedure(mgbf_proc)

!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!      if(mype==0) then
!         print *,'finish filtering_procdure'
!      endif
!         call finishMPI
!TEST->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!======================================================================

!***
!*** From first generation of filter grid to analysis grid
!***

                                                   call btim(   filt2an_tim)
        if(l_new_map) then
          call filt_to_anal_all2
        else 
          call filt_to_anal_all
        endif

                                                   call etim(   filt2an_tim)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!***
!*** Adjoint test if needed
!***
      
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
   
      
!==================== Forward (Smoothing step) ========================
!***
!*** DONE! Deallocate variables
!***
                                                   call btim(   output_tim)
       call mg_finalize

                                                   call etim(   output_tim)
                                                   call etim(   total_tim)


!***
!*** Print wall clock and cpu timing
!***
      call print_mg_timers("timing_cpu.csv", print_cpu)



      call finishMPI


!-----------------------------------------------------------------------
                        endprogram RBETA_TEST 
