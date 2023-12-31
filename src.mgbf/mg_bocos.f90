!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        module mg_bocos
!***********************************************************************
!                                                                      !
!  Provide communication between subdomains and supply halos on        !
!  filter grid                                                         !
!                       - offset version -                             !
!                                                                      !
! Libraries: mpi                                                       !
! Modules: kinds, mg_mppstuff, mg_parameter, mg_domain                 !
!                                                     M. Rancic (2022) !
!***********************************************************************
!use mpi
use kinds, only: r_kind,i_kind
!use mpimod, only: mype,mpi_comm_world
use mg_mppstuff, only: mype,mpi_comm_comp,mpi_comm_work
use mg_mppstuff, only: itype,rtype,dtype,my_hgen                        &
                      ,barrierMPI,finishMPI,l_hgen,mype_hgen
use mg_parameter, only: gm

implicit none

interface boco_2d
  module procedure boco_2d_g1 
  module procedure boco_2d_gh 
endinterface

interface bocoT_2d
  module procedure bocoT_2d_g1 
  module procedure bocoT_2d_gh 
endinterface

interface boco_3d
  module procedure boco_3d_g1 
  module procedure boco_3d_gh 
endinterface

interface bocoT_3d
  module procedure bocoT_3d_g1 
  module procedure bocoT_3d_gh 
endinterface

interface upsend_all
  module procedure upsend_all_g1
  module procedure upsend_all_gh
endinterface

interface downsend_all     
  module procedure downsend_all_gh
  module procedure downsend_all_g2
endinterface

interface bocox
  module procedure bocox_2d_g1
  module procedure bocox_2d_gh
endinterface

interface bocoy
  module procedure bocoy_2d_g1
  module procedure bocoy_2d_gh
endinterface

interface bocoTx
  module procedure bocoTx_2d_g1
  module procedure bocoTx_2d_gh
endinterface

interface bocoTy
  module procedure bocoTy_2d_g1
  module procedure bocoTy_2d_gh
endinterface

public boco_2d_loc
public bocoT_2d_loc

public upsend_loc_g12
public upsend_loc_g23
public upsend_loc_g34

public downsend_loc_g43
public downsend_loc_g32
public downsend_loc_g21

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_2d_g1                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for generation 1                 !
!                                                                      !
!                        - offset version -                            !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                     

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay,nbxy
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatay = km*imax*nby
      ndatax = km*(jmax+2*nby)*nbx

!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!
! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif


! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

!----------------------------------------------------------------------
!
! SEND extended boundaries toward WEST and EAST 
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST 
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if


!
! Assign received values from EAST and WEST
!

! From west

   if(lwest) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_2d_gh                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for high generations             !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*imax*nby
      ndatax = km*(jmax+2*nby)*nbx


!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!


! From south

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif


! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

!
!  SEND extended boundaries to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_2d_g1                          &
!***********************************************************************
!                                                                      !
! Adjoint of side sending subroutine:                                  !
! Supplies (nbx,nby) lines of halos in (x,y) directions, including     !
! values at the edges of the subdomains and assuming mirror boundary   !
! conditions just for generation 1                                     !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatax =km*(jmax+2*nby)*nbx
      ndatay =km*imax*nby

!
! SEND extended halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j) = W(:,-nbx+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!

! From west

   if(lwest) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do
   else 
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND boundaries SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby), stat = iaerr )

              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_comp, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! ASSIGN received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!-----------------------------------------------------------------------
                        endsubroutine bocoT_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_2d_gh                          &
!***********************************************************************
!                                                                      !
!  Supply n-lines inside of domains, including edges, with halos from  !
!  the surrounding domains.  Assume mirror boundary conditions at the  !
!  boundaries of the domain. For high multigrid generations.           !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+2*nby)*nbx
      ndatay =km*imax*nby

!
! SEND extended halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j) = W(:,-nbx+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do
   else 
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby  
              do i=1,imax
                sBuf_S(:,i,j) = W(:,i,-nby+j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoT_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_3d_g1                           &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for generation 1                 !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km3,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,im,jm,Lm,nbx,nby,nbz
real(r_kind),dimension(km3,1-nbx:im+nbx,1-nby:jm+nby,1-nbz:Lm+nbz)    &
                      ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to generation one
!
          g_ind=1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im  
          jmax = jm

!-----------------------------------------------------------------------
      ndatay = km3*imax*nby*Lm
      ndatax = km3*(jmax+2*nby)*nbx*Lm


!
! SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km3,1:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km3,1:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j,L)=W(:,i,jmax-nby+j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km3,1:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km3,1:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j,L)=W(:,i,jmax+1-j,L)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j,L)=rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j,L)=W(:,i,nby+1-j,L)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j,L)=rBuf_S(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!
! SEND extended boundaries toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j,L) = W(:,imax-nbx+i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries WEST and EAST
!

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

!
! Assign received values from  EAST and WEST
!
! From west

   if(lwest) then

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j,L)= W(:,nbx+1-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j,L)= rBuf_W(:,i,j,L)
     enddo
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=W(:,imax-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=rBuf_E(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff


!-----------------------------------------------------------------------
                        endsubroutine boco_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_3d_gh                           &
!**********************************************************************!

! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for high generations             !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km3,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                &
                    ,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,im,jm,Lm,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind),dimension(km3,1-nbx:im+nbx,1-nby:jm+nby,1-nbz:Lm+nbz)    &
                      ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km3*imax*nby*Lm
      ndatax = km3*(jmax+2*nby)*nbx*Lm

!
! SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km3,1:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km3,1:imax,nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j,L)=W(:,i,jmax-nby+j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
! RECEIVE boundaries from SOUTH and NORTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km3,1:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km3,1:imax,nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!TEST
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
!TEST

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j,L)=W(:,i,jmax+1-j,L)
     enddo
     enddo
     enddo

   else

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j,L)=rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j,L)=W(:,i,nby+1-j,L)
     end do
     end do
     end do

   else

     do L=1,Lm
     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j,L)=rBuf_S(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!TEST
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
!TEST


!
! SEND extended boundaries to WEST and EAST   
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j,L) = W(:,i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=1,Lm
                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j,L) = W(:,imax-nbx+i,j,L)
                  enddo
                enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km3,nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

!
! Deallocate send bufferes from EAST and WEST
!
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if

!
! Assign received values from WEST and EAST
!
! From west

   if(lwest) then

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j,L)= W(:,nbx+1-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j,L)= rBuf_W(:,i,j,L)
     enddo
     enddo
     enddo


   endif

! From east

   if(least) then

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=W(:,imax+1-i,j,L)
     end do
     end do
     end do

   else 

     do L=1,Lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j,L)=rBuf_E(:,i,j,L)
     enddo
     enddo
     enddo

   endif

!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1-L )=W(:,:,:, 1+L)
          W(:,:,:,LM+L)=W(:,:,:,LM-L)
        end do


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_3d_g1                          &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain                                            *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(W,km3,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km3,im,jm,Lm,nbx,nby,nbz
real(r_kind), dimension(km3,1-nbx:im+nbx,1-nby:jm+nby,1-nbz:Lm+nbz)   &
                       ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

         g_ind=1

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          imax = im
          jmax = jm

!----------------------------------------------------------------------
      ndatax =km3*(jmax+2*nby)*nbx *Lm
      ndatay =km3*imax*nby *Lm

!
! SEND extended halos toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km3,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j,L) = W(:,-nbx+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km3,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j,L) = W(:,imax+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(2), isend)

      end if
!
! RECEIVE extended halos from EAST and WEST
!
! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km3,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


          allocate( rBuf_W(1:km3,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if
!
! Assign received extended halos from WEST and EAST to interior of domains
!

! From west

   if(lwest) then
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j,L)= W(:,i,j,L)+W(:,1-i,j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
     end do
     end do
     end do
   endif

! From east

   if(least) then
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+nbx-i,j,L)
     end do
     end do
     end do
   else 
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
     end do
     end do
     end do
   endif

!
! Send halos SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km3,1:imax,1:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km3,1:imax,1:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j,L)=W(:,i,jmax+j,L)
              enddo
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_comp, sHandle(1), isend)

      end if


!
! RECEIVE boundaries from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km3,1:imax,1:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km3,1:imax,1:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+W(:,i,1-j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
     end do
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+nby-j,L)
     enddo
     enddo
     enddo
   else
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo
   endif

!----------------------------------------------------------------------
!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1+L )=W(:,:,:, 1+L)+W(:,:,:, 1-L)
          W(:,:,:,LM-L)=W(:,:,:,LM-L)+W(:,:,:,LM+L)
        end do


!----------------------------------------------------------------------
!
!                           DEALLOCATE sBufferes
!


      if( itarg_w >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif 
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif 
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif 
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif


!-----------------------------------------------------------------------
                        endsubroutine bocoT_3d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_3d_gh                          &
!***********************************************************************
!                                                                      *
!  Supply n-lines inside of domains, including edges, with halos from  *
!  the surrounding domains.  Assume mirror boundary conditions at the  *
!  boundaries of the domain                                            *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(W,km,im,jm,Lm,nbx,nby,nbz,Fimax,Fjmax,mygen_min,mygen_max)

!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flsouth,Flnorth                      &
                    ,Fitarg_n,Fitarg_s,Fitarg_w,Fitarg_e                
use mpi
implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,Lm,nbx,nby,nbz,mygen_min,mygen_max
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby,1-nbz:Lm+nbz)    &
                       ,intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:,:)::                         &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+2*nby)*nbx *Lm
      ndatay =km*imax*nby *Lm

!
! SEND extended halos toward WEST and EAST
!
! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j,L) = W(:,-nbx+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j,L) = W(:,imax+i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)
      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


          allocate( rBuf_W(1:km,1:nbx,1-nby:jmax+nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )


      end if

!
! Assign received extended halos from WEST and EAST
!

! From west

   if(lwest) then
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j,L)= W(:,i,j,L)+W(:,1-i,j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j,L)= W(:,i,j,L)+rBuf_W(:,i,j,L)
     end do
     end do
     end do
   endif

! From east

   if(least) then
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+W(:,imax+1+nbx-i,j,L)
     end do
     end do
     end do
   else 
     do L=1,lm
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j,L)= W(:,imax-nbx+i,j,L)+rBuf_E(:,i,j,L)
     end do
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby,L) = W(:,i,j,L)
              enddo
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby,1:Lm), stat = iaerr )

              do L=Lm,1,-1
              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j,L)=W(:,i,jmax+j,L)
              enddo
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby,1:Lm), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if


!-----------------------------------------------------------------------
!
! Assign received halos from SOUTH and NORTH
!

   if(lsouth) then
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+W(:,i,1-j,L)
     end do
     end do
     end do
   else
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,j,L)= W(:,i,j,L)+rBuf_S(:,i,j,L)
     end do
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+W(:,i,jmax+1+nby-j,L)
     enddo
     enddo
     enddo
   else
     do L=1,lm
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j,L)= W(:,i,jmax-nby+j,L)+rBuf_N(:,i,j,L)
     enddo
     enddo
     enddo
   endif


!
! Set up mirror b.c. at the bottom and top of domain 
!
        do L=1,nbz
          W(:,:,:,1+L )=W(:,:,:, 1+L)+W(:,:,:, 1-L)
          W(:,:,:,LM-L)=W(:,:,:,LM-L)+W(:,:,:,LM+L)
        end do


!-----------------------------------------------------------------------
!
!                           DEALLOCATE sBufferes
!

      if( itarg_w >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
         deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
         deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
         deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
         deallocate( sBuf_N, stat = ierr )
      end if
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoT_3d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_g1                        &
!***********************************************************************
!                                                                      !
!         Upsend data from generation one to generation two            !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(Harray,Warray,km)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km
real(r_kind), dimension(km,1:imL,1:jmL),intent(in):: Harray
real(r_kind), dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(out):: Warray

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!
       g_ind=1

       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)


       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:) = 0.0d0
   endif

     ndata =km*imL*jmL

!
! --- Send data to SW portion of processors at higher generation
!

      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      endif
!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
                Warray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      endif

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then

        nebpe = itargdn_se

        if(nebpe /= mype) then

          call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

        endif
             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               Warray(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                      mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,          &
                         mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

      endif


!-----------------------------------------------------------------------
                        endsubroutine upsend_all_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_all_gh                        &
!***********************************************************************
!                                                                      *
!         Upsend data from one grid generation to another              *
!         (Just for high grid generations)                             *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(Harray,Warray,km,mygen_dn,mygen_up)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne    &
                    ,Fitarg_up                                          &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km
real(r_kind), dimension(km,1:imL,1:jmL),intent(in):: Harray
real(r_kind), dimension(km,1-hx:im+hx,1-hy:jm+hy),intent(out):: Warray
integer(i_kind),intent(in):: mygen_dn,mygen_up

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne,flag_up
integer(i_kind):: itarg_up
integer:: g_ind

!-----------------------------------------------------------------------
!
! Define generational flags
!
 
       g_ind=2 

       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)                                          


!-----------------------------------------------------------------------

   if(my_hgen==mygen_up) then
      Warray(:,:,:)=0.0d0
   endif

     ndata =km*imL*jmL
!TEST
!    if(mype==0) then 
!       write(0,*) 'From upsend_all_gh.f90: ndata=',ndata
!    endif
!TEST


      if(  lsendup_sw ) then

        nebpe = itarg_up
    

        allocate( sBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = Harray(:,i,j)
             enddo
             enddo
!TEST
!       write(0,*) 'UPSEND_ALL_GH SW: ndata,mype,nepbe=',ndata,mype,nebpe
!TEST

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )


      end if

!
! --- Receive SW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        allocate( rBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Warray(:,i,j)=Rbuf_SW(:,i,j)
             enddo
             enddo

      endif
!TEST
!         call finishMPI
!TEST


!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( sBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = Harray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_work, sHandle(2), isend)

        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

      end if

!
! --- Receive SE portion of data at higher generation


      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
        nebpe = itargdn_se


        allocate( rBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

      endif


!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        allocate( sBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NW(:,i,j) = Harray(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )


    end if

!
! --- Receive NW portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw
 

        allocate( rBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)

        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Warray(:,i,jmL+j)=rBuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)

      end if

!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( sBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NE(:,i,j) = Harray(:,i,j)
             enddo
             enddo
!TEST
!       write(0,*) 'UPSEND_ALL_GH NE: ndata,mype,nepbe,mygen_up',ndata,mype,nebpe,mygen_up
!TEST

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                       mpi_comm_work, sHandle(4), isend)

         call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

      end if

!
! --- Receive NE portion of data at higher generation
!

      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne

        allocate( rBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)

        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Warray(:,imL+i,jmL+j)=rBuf_NE(:,i,j)
             enddo
             enddo

          deallocate( rBuf_NE, stat = iderr)

      endif

!TEST
!         call finishMPI
!TEST

!-----------------------------------------------------------------------
                        endsubroutine upsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_all_gh                      &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(Warray,Harray,km,mygen_up,mygen_dn)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne     &
                    ,Fitarg_up                                           &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw       
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km
real(r_kind), dimension(km,1:im,1:jm),intent(in):: Warray
real(r_kind), dimension(km,1:imL,1:jmL),intent(out):: Harray
integer, intent(in):: mygen_up,mygen_dn
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
!-----------------------------------------------------------------------

     Harray(:,:,:) = 0.0d0
!
! Define generational flags
!

         g_ind=2
       lsendup_sw=Flsendup_sw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_se=Flsendup_se(g_ind).and.(my_hgen==mygen_dn)
       lsendup_nw=Flsendup_nw(g_ind).and.(my_hgen==mygen_dn)
       lsendup_ne=Flsendup_ne(g_ind).and.(my_hgen==mygen_dn)

       itarg_up=Fitarg_up(g_ind)

       ndata =km*imL*jmL

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation

 
  if(my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
        nebpe = itargdn_sw


        allocate( sBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )


  endif
!
! --- Receive SW portion of data at lower generation


      if( lsendup_sw ) then

        nebpe = itarg_up


        allocate( rBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

      endif

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(my_hgen==mygen_up .and.  itargdn_se >= 0 ) then
        nebpe = itargdn_se

        allocate( sBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )


  endif
!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( rBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
     end if

!
! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

  if(my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
        nebpe = itargdn_nw


        allocate( sBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )


  endif
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        allocate( rBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)


      end if


! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if(my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
        nebpe = itargdn_ne


        allocate( sBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )


  endif
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( rBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

      end if

!-----------------------------------------------------------------------
                        endsubroutine downsend_all_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_all_g2                      &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      *
!                       - offset version -                             *
!                                                                      *
!***********************************************************************
(Warray,Harray,km)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain, only: Flsendup_sw,Flsendup_se,Flsendup_nw,Flsendup_ne     &
                    ,Fitarg_up                                           &
                    ,itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw       
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km
real(r_kind), dimension(km,1:im,1:jm),intent(in):: Warray
real(r_kind), dimension(km,1:imL,1:jmL),intent(out):: Harray
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE             

real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe

logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne  
integer:: mygen_up,mygen_dn
integer(i_kind):: itarg_up                                           
integer(i_kind):: g_ind
!-----------------------------------------------------------------------
!
! Define generational flags
!
    mygen_up=2
    mygen_dn=1

         g_ind=1
       lsendup_sw=Flsendup_sw(g_ind)
       lsendup_se=Flsendup_se(g_ind)
       lsendup_nw=Flsendup_nw(g_ind)
       lsendup_ne=Flsendup_ne(g_ind)

       itarg_up=Fitarg_up(g_ind)


      ndata =km*imL*jmL


!
! Send data down to generation 1
!
LSEND:  if(my_hgen==mygen_up) then
!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation
 
        nebpe = itargdn_sw

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = Warray(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )

        endif
!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

        nebpe = itargdn_se

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = Warray(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

        endif

! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

        nebpe = itargdn_nw

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NW(:,i,j) = Warray(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )

        endif

!
! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

        nebpe = itargdn_ne
        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:km,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NE(:,i,j) = Warray(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

        endif


    endif LSEND   

!
! --- Receive SW portion of data at lower generation
!

      if( lsendup_sw .and. mype /= itarg_up ) then

        nebpe = itarg_up


        call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )


      else &

!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se .and. mype /= itarg_up) then

        nebpe = itarg_up

        call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )


      else &


!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw .and. mype /= itarg_up) then

        nebpe = itarg_up

        call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_comp, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )


      else &


!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne .and. mype /= itarg_up) then
        nebpe = itarg_up

        call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_comp, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )


      end if
   
!
! Assign received and prescribed values
!     
      if( lsendup_sw ) then

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_se ) then

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_nw ) then

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NW(:,i,j)
             enddo
             enddo

      else &
      if( lsendup_ne ) then

             do j=1,jmL
             do i=1,imL
               Harray(:,i,j)=dBuf_NE(:,i,j)
             enddo
             enddo

       endif


!-----------------------------------------------------------------------
                        endsubroutine downsend_all_g2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocox_2d_g1                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction assuming mirror boundary  !
! conditions at the end of domain. Version for generation 1            !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_w,Fitarg_e,Flwest,Fleast

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatax = km*jmax*nbx

!----------------------------------------------------------------------
!
! SEND extended boundaries toward WEST and EAST 
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE boundaries from EAST and WEST 
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if


!
! Assign received values from EAST and WEST
!

! From west

   if(lwest) then

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif


!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocox_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocox_2d_gh                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines of halos in x direction assuming mirror boundary  !
! conditions at the end of domain. Version for high generations        !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_w,Fitarg_e,Flwest,Fleast,Flsouth,Flnorth

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W           

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatax = km*jmax*nbx

!
!  SEND halos to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,1:jmax), stat = iaerr )

                do j=1,jmax
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax-i,j)
     end do
     end do

   else 

     do j=1,jmax
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocox_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_g1                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction assuming mirror boundary  !
! conditions at the end of domain. Version for generation 1            !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flsouth,Flnorth                     

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          g_ind = 1

          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 

          imax = im       
          jmax = jm


!-----------------------------------------------------------------------
      ndatay = km*imax*nby


!
!  SEND boundaries toward SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_comp, sHandle(1), isend)

      end if

!
! RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if
!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From SOUTH

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif



!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoy_2d_gh                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines of halos in y direction assuming mirror boundary  !
! conditions at the end of domain. Version for high generations        !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fitarg_n,Fitarg_s,Flwest,Fleast,Flsouth,Flnorth                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
integer(i_kind) g_ind,g
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!
 
       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2 
         g = my_hgen
         l_sidesend=.true.
       else
         l_sidesend=.false.
       endif

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)                 


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*imax*nby

!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!

! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

! From south

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      endif
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoy_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_g1                         &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines close to edges of the subdomins from neighboring  !
! halos in x direction assuming mirror boundary conditions             !
! Version for generation 1                                             !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   

integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=1
!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatax =km*jmax*nbx

!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(1), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e


          allocate( rBuf_E(1:km,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w


         allocate( rBuf_W(1:km,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received halos from WEST and EAST to interrior of domains
!

! From west

   if(lwest) then
     do j=1,jmax
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1,jmax
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1,jmax
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx+1-i,j)
     end do
     end do
   else 
     do j=1,jmax
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_w  >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      endif
      if( itarg_e  >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      endif

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTx_2d_gh                         &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nbx lines close to edges of the subdomins from neighboring  !
! halos in x direction assuming mirror boundary conditions             !
! Version for high generations                                         ! 
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Flwest,Fleast,Flnorth,Fitarg_w,Fitarg_e                
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_E,sBuf_W             &
                                             ,rBuf_E,rBuf_W   
integer(i_kind) itarg_w,itarg_e,imax,jmax
logical lwest,least,lnorth

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_w = Fitarg_w(g_ind)
          itarg_e = Fitarg_e(g_ind)

          lwest   = Flwest(g_ind)
          least   = Fleast(g_ind)

          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*jmax*nbx
!
! SEND halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1:jmax), stat = iaerr )

              do j=1,jmax
              do i=1-nbx,0
                sBuf_W(:,i+nbx,j) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,0:nbx,0:jmax), stat = iaerr )

              do j=1,jmax
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,1:nbx,1:jmax), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=1,jmax
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1,jmax
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1,jmax
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+nbx+1-i,j)
     end do
     end do
   else 
     do j=1,jmax
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

      if( itarg_w  >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e  >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTx_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTy_2d_g1                         &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines close to edges of the subdomins from neighboring  !
! halos in y direction assuming mirror boundary conditions             !
! Version for generation 1                                             !                                             
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby)
!-----------------------------------------------------------------------
use mg_domain, only: Flsouth,Flnorth,Fitarg_n,Fitarg_s
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:):: sBuf_N,sBuf_S             &
                                             ,rBuf_N,rBuf_S

integer(i_kind) itarg_n,itarg_s,imax,jmax
logical lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!

         g_ind=1
!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)

          imax = im    
          jmax = jm


!----------------------------------------------------------------------
      ndatay =km*imax*nby

!
! SEND SOUTH and NORTH halos
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby), stat = iaerr )

              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_comp, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_comp, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )


      end if

!
! ASSIGN received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby+1-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!

      if( itarg_s  >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n  >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_g1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoTy_2d_gh                         &
!***********************************************************************
!                                                                      !
! Side sending subroutine:                                             !
! Supplies nby lines close to edges of the subdomins from neighboring  !
! halos in y direction assuming mirror boundary conditions             !
! Version for high generations                                         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,mygen_min,mygen_max)
!-----------------------------------------------------------------------
use mg_domain, only: Fleast,Flsouth,Flnorth                             &
                    ,Fitarg_n,Fitarg_s
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,mygen_min,mygen_max
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::  sBuf_N,sBuf_S            &
                                              ,rBuf_N,rBuf_S
integer(i_kind) itarg_n,itarg_s,itarg_e,imax,jmax
logical least,lsouth,lnorth                                       

integer(i_kind) sHandle(2),rHandle(2),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatay
logical l_sidesend
integer(i_kind) g_ind,g,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


       if(mygen_min <= my_hgen .and. my_hgen <= mygen_max) then
         g_ind=2
         g = my_hgen
         l_sidesend=.true.
       else 
         l_sidesend=.false.
       endif


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n(g_ind)
          itarg_s = Fitarg_s(g_ind)

          least   = Fleast(g_ind)
          lsouth  = Flsouth(g_ind)
          lnorth  = Flnorth(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------

      ndatay =km*imax*nby
!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby), stat = iaerr )

              do j=1-nby,0
              do i=1,imax
                sBuf_S(:,i,j+nby) = W(:,i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(1), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+nby+1-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

      if( itarg_s  >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n  >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!                           DEALLOCATE sBufferes

      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoTy_2d_gh

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine boco_2d_loc                          &
!**********************************************************************!
!                                                                      !
! Side sending subroutine:                                             !
! Supplies (nbx,nby) lines of halos in (x,y) directions assuming       !
! mirror boundary conditions. Version for localiztion                  ! 
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!**********************************************************************!
(W,km,im,jm,nbx,nby,Fimax,Fjmax,g)
!-----------------------------------------------------------------------
use mg_domain_loc, only: Fitarg_n_loc,Fitarg_s_loc,Fitarg_w_loc,Fitarg_e_loc &
                        ,Flwest_loc,Fleast_loc,Flsouth_loc,Flnorth_loc                      

use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,g
real(r_kind),dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                  sBuf_N,sBuf_E,sBuf_S,sBuf_W           &
                                 ,rBuf_N,rBuf_E,rBuf_S,rBuf_W           

integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical:: lwest,least,lsouth,lnorth                                      

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,l,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
integer(i_kind) g_ind
logical l_sidesend
!-----------------------------------------------------------------------
!
! Limit communications to selected number of generations
!

         l_sidesend=.true.

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff
!
! from mg_domain      
! 
          itarg_n = Fitarg_n_loc(g)
          itarg_s = Fitarg_s_loc(g)
          itarg_w = Fitarg_w_loc(g)
          itarg_e = Fitarg_e_loc(g)

          lwest   = Flwest_loc(g)
          least   = Fleast_loc(g)
          lsouth  = Flsouth_loc(g)
          lnorth  = Flnorth_loc(g)                 


!
!  Keep this for now but use only Mod(nxm,8)=Mod(nym,8)=0
!

          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!-----------------------------------------------------------------------
      ndatay = km*imax*nby
      ndatax = km*(jmax+2*nby)*nbx


!
!  SEND boundaries to SOUTH and NORTH
!

! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

            allocate( sBuf_S(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_S(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

            allocate( sBuf_N(1:km,1:imax,nby), stat = iaerr )

                do j=1,nby
                  do i=1,imax
                    sBuf_N(:,i,j)=W(:,i,jmax-nby+j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                              mpi_comm_work, sHandle(1), isend)

      end if
!
!     RECEIVE boundaries from NORTH and SOUTH
!

! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

          allocate( rBuf_N(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe, &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

          allocate( rBuf_S(1:km,1:imax,nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )

      end if

!
! Assign received values from NORTH and SOUTH
!


! From south

   if(lsouth) then

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=W(:,i,nby+1-j)
     end do
     end do

   else

     do j=1,nby
     do i=1,imax
       W(:,i,-nby+j)=rBuf_S(:,i,j)
     enddo
     enddo

   endif


! --- from NORTH ---

   if( lnorth) then

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=W(:,i,jmax+1-j)
     enddo
     enddo

   else

     do j=1,nby
     do i=1,imax
       W(:,i,jmax+j)=rBuf_N(:,i,j)
     enddo
     enddo

   endif

!
!  SEND extended boundaries to WEST and EASTH
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_W(:,i,j) = W(:,i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )

                do j=1-nby,jmax+nby
                  do i=1,nbx
                    sBuf_E(:,i,j) = W(:,imax-nbx+i,j)
                  enddo
                enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype, &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
!     RECEIVE extended boundaries from EAST and WEST
!

! --- from EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if( itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,  &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if

!
! Assign received values from  WEST and EAST
!

! From west

   if(lwest) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= W(:,nbx+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,-nbx+i,j)= rBuf_W(:,i,j)
     enddo
     enddo


   endif

! From east

   if(least) then

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=W(:,imax+1-i,j)
     end do
     end do

   else 

     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax+i,j)=rBuf_E(:,i,j)
     enddo
     enddo

   endif

!-----------------------------------------------------------------------
!
!                           DEALLOCATE rBufferes
!
      if( itarg_w >= 0 ) then
        deallocate( rBuf_W, stat = iderr)
      end if
      if( itarg_e >= 0 ) then
        deallocate( rBuf_E, stat = iderr)
      end if
      if( itarg_s >= 0 ) then
        deallocate( rBuf_S, stat = iderr)
      end if
      if( itarg_n >= 0 ) then
        deallocate( rBuf_N, stat = iderr)
      end if

!
!                           DEALLOCATE sBufferes
!

      if( itarg_w >= 0 ) then
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_W, stat = ierr )
      end if
      if( itarg_e >= 0 ) then
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_E, stat = ierr )
      end if
      if( itarg_s >= 0 ) then
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_S, stat = ierr )
      end if
      if( itarg_n >= 0 ) then
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_N, stat = ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine boco_2d_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine bocoT_2d_loc                         &
!***********************************************************************
!                                                                      !
!  Supply n-lines inside of domains, including edges, with halos from  !
!  the surrounding domains.  Assume mirror boundary conditions at the  !
!  boundaries of the domain. Vesrion for localization.                 !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(W,km,im,jm,nbx,nby,Fimax,Fjmax,g)
!-----------------------------------------------------------------------
use mg_domain_loc, only: Flwest_loc,Fleast_loc,Flsouth_loc,Flnorth_loc      &
                    ,Fitarg_w_loc,Fitarg_e_loc,Fitarg_s_loc,Fitarg_n_loc             
use mpi

implicit none

!-----------------------------------------------------------------------
integer(i_kind), intent(in):: km,im,jm,nbx,nby,g
real(r_kind), dimension(km,1-nbx:im+nbx,1-nby:jm+nby),intent(inout):: W
integer(i_kind), dimension(gm), intent(in):: Fimax,Fjmax
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                           &
                                        sBuf_N,sBuf_E,sBuf_S,sBuf_W     &
                                       ,rBuf_N,rBuf_E,rBuf_S,rBuf_W   
integer(i_kind) itarg_n,itarg_s,itarg_w,itarg_e,imax,jmax
logical lwest,least,lsouth,lnorth                                       

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,L,i,j
integer(i_kind) isend,irecv,nebpe
integer(i_kind) ndatax,ndatay
logical l_sidesend
integer(i_kind) g_ind,k
!-----------------------------------------------------------------------
!
! Limit comminications to selected number of generations
!


         g_ind=g
         l_sidesend=.true.


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

FILT_GRID:    if(l_sidesend) then

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!
! from mg_domain
!
          itarg_n = Fitarg_n_loc(g_ind)
          itarg_s = Fitarg_s_loc(g_ind)
          itarg_w = Fitarg_w_loc(g_ind)
          itarg_e = Fitarg_e_loc(g_ind)

          lwest   = Flwest_loc(g_ind)
          least   = Fleast_loc(g_ind)
          lsouth  = Flsouth_loc(g_ind)
          lnorth  = Flnorth_loc(g_ind)


          if(least) then
            imax = Fimax(g)
          else 
            imax = im       !   << Note that is not necesseraly im from
          endif             !      mg_parameter.  Could be also imL >>>
          if(lnorth) then
            jmax = Fjmax(g)
          else  
            jmax = jm
          endif


!----------------------------------------------------------------------
      ndatax =km*(jmax+2*nby)*nbx
      ndatay =km*imax*nby

!
! SEND extended halos toward WEST and EAST
!

! --- toward WEST ---

      if( itarg_w >= 0) then
        nebpe = itarg_w

              allocate( sBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_W(:,i,j) = W(:,-nbx+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_W, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(4), isend)

      end if

! --- toward EAST ---

      if( itarg_e >= 0 ) then
        nebpe = itarg_e

              allocate( sBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )

              do j=1-nby,jmax+nby
              do i=1,nbx
                sBuf_E(:,i,j) = W(:,imax+i,j)
              enddo
              enddo

              call MPI_ISEND( sBuf_E, ndatax, dtype, nebpe, mype,       &
                              mpi_comm_work, sHandle(2), isend)

      end if

!
! RECEIVE extended halos from EAST and WEST
!

! --- from EAST ---

      if(  itarg_e >= 0 ) then
        nebpe = itarg_e

          allocate( rBuf_E(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_E, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

      end if

! --- from WEST ---

      if(  itarg_w >= 0 ) then
        nebpe = itarg_w

          allocate( rBuf_W(1:km,1:nbx,1-nby:jmax+nby), stat = iaerr )
          call MPI_IRECV( rBuf_W, ndatax, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )

      end if
!
! Assign received values from WEST and EAST
!

! From west

   if(lwest) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,i,j)= W(:,i,j)+W(:,1-i,j)
     end do
     end do
   else
     do j=1-nby,jmax+nby
     do i=1,nbx
      W(:,i,j)= W(:,i,j)+rBuf_W(:,i,j)
     end do
     end do
   endif

! From east

   if(least) then
     do j=1-nby,jmax+nby
     do i=1,nbx
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+W(:,imax+1+nbx-i,j)
     end do
     end do
   else 
     do j=1-nby,jmax+nby
     do i=1,nbx  
       W(:,imax-nbx+i,j)= W(:,imax-nbx+i,j)+rBuf_E(:,i,j)
     end do
     end do
   endif

!
! SEND halos toward SOUTH and NORTH
!
! --- toward SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s

              allocate( sBuf_S(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby  
              do i=1,imax
                sBuf_S(:,i,j) = W(:,i,-nby+j)
              enddo
              enddo

              call MPI_ISEND( sBuf_S, ndatay, dtype, nebpe, mype,  &
                              mpi_comm_work, sHandle(3), isend)
      end if

! --- toward NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n

             allocate( sBuf_N(1:km,1:imax,1:nby), stat = iaerr )

              do j=1,nby
              do i=1,imax
                sBuf_N(:,i,j)=W(:,i,jmax+j)
              enddo
              enddo

             call MPI_ISEND( sBuf_N, ndatay, dtype, nebpe, mype,        &
                             mpi_comm_work, sHandle(1), isend)

      end if

!
! RECEIVE halos from NORTH and SOUTH
!
!
! --- from NORTH ---

      if( itarg_n >= 0 ) then
        nebpe = itarg_n


          allocate( rBuf_N(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_N, ndatay, dtype, nebpe, nebpe,          &
                      mpi_comm_work, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )

      end if

! --- from SOUTH ---

      if( itarg_s >= 0 ) then
        nebpe = itarg_s


          allocate( rBuf_S(1:km,1:imax,1:nby), stat = iaerr )
          call MPI_IRECV( rBuf_S, ndatay, dtype, nebpe, nebpe,          &
                       mpi_comm_work, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )


      end if

!
! Assign received values from SOUTH and NORTH
!

! From south

   if(lsouth) then
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+W(:,i,1-j)
     end do
     end do
   else
     do j=1,nby
     do i=1,imax
       W(:,i,j)= W(:,i,j)+rBuf_S(:,i,j)
     end do
     end do
   endif

!  From north

   if(lnorth) then
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+W(:,i,jmax+1+nby-j)
     enddo
     enddo
   else
     do j=1,nby
     do i=1,imax
       W(:,i,jmax-nby+j)= W(:,i,jmax-nby+j)+rBuf_N(:,i,j)
     enddo
     enddo
   endif

!-----------------------------------------------------------------------

!                           DEALLOCATE rBufferes

        deallocate( rBuf_W, stat = iderr)
        deallocate( rBuf_E, stat = iderr)
        deallocate( rBuf_S, stat = iderr)
        deallocate( rBuf_N, stat = iderr)

!                           DEALLOCATE sBufferes

      if( itarg_w  >= 0 ) then
         call MPI_WAIT( sHandle(4), istat, ierr )
      end if
      if( itarg_e  >= 0 ) then
         call MPI_WAIT( sHandle(2), istat, ierr )
      end if
      if( itarg_s  >= 0 ) then
         call MPI_WAIT( sHandle(3), istat, ierr )
      end if
      if( itarg_n  >= 0 ) then
         call MPI_WAIT( sHandle(1), istat, ierr )
      end if


!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

     endif FILT_GRID

!fffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff

!-----------------------------------------------------------------------
                        endsubroutine bocoT_2d_loc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_loc_g12                       &
!***********************************************************************
!                                                                      !
!         Upsend data from generation one to generation two            !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(V,H,km_4,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: lsendup_sw_loc,lsendup_se_loc                  &
                        ,lsendup_nw_loc,lsendup_ne_loc                  &
                        ,Fitargup_loc12                                &
                        ,itargdn_sw_loc21,itargdn_se_loc21              &
                        ,itargdn_ne_loc21,itargdn_nw_loc21
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_4,flag
real(r_kind), dimension(km_4,1:imL,1:jmL),intent(in):: V
real(r_kind), dimension(km_4,1-hx:im+hx,1-hy:jm+hy),intent(out):: H

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
integer(i_kind):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne

!-----------------------------------------------------------------------
   mygen_dn=1
   mygen_up=2
!
! Define generational flags
!

       itarg_up=Fitargup_loc12(flag)


       itargdn_sw = itargdn_sw_loc21
       itargdn_se = itargdn_se_loc21                  
       itargdn_ne = itargdn_ne_loc21
       itargdn_nw = itargdn_nw_loc21

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

!N   if(my_hgen==mygen_up) then
      H(:,:,:) = 0.0d0
!N   endif

     ndata =km_4*imL*jmL

!
! --- Send data to SW portion of processors at higher generation
!

      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      endif
!
! --- Receive SW portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
      if( itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      endif

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
      if( itargdn_se >= 0 ) then

        nebpe = itargdn_se

        if(nebpe /= mype) then

          call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

        endif
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

!      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
      if( itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                      mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
      if( itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,          &
                         mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

      endif


!-----------------------------------------------------------------------
                        endsubroutine upsend_loc_g12

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_loc_g23                       &
!***********************************************************************
!                                                                      !
!         Upsend data from generation three to generation four         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(V,H,km_16,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: lsendup_sw_loc,lsendup_se_loc                  &
                        ,lsendup_nw_loc,lsendup_ne_loc                  &
                        ,Fitargup_loc23                                &
                        ,itargdn_sw_loc32,itargdn_se_loc32              &
                        ,itargdn_ne_loc32,itargdn_nw_loc32    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_16,flag
real(r_kind), dimension(km_16,1:imL,1:jmL),intent(in):: V
real(r_kind), dimension(km_16,1-hx:im+hx,1-hy:jm+hy),intent(out):: H

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
integer(i_kind):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne

!-----------------------------------------------------------------------
   mygen_dn=2
   mygen_up=3
!
! Define generational flags
!

       itarg_up=Fitargup_loc23(flag)

       itargdn_sw = itargdn_sw_loc32
       itargdn_se = itargdn_se_loc32                  
       itargdn_ne = itargdn_ne_loc32
       itargdn_nw = itargdn_nw_loc32

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

!N   if(my_hgen==mygen_up) then
      H(:,:,:) = 0.0d0
!N   endif

     ndata =km_16*imL*jmL

!
! --- Send data to SW portion of processors at higher generation
!

      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      endif
!
! --- Receive SW portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_sw >= 0 ) then
      if( itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      endif

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_se >= 0 ) then
      if( itargdn_se >= 0 ) then

        nebpe = itargdn_se

        if(nebpe /= mype) then

          call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

        endif
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

!      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
      if( itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                      mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
      if( itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,          &
                         mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

      endif


!-----------------------------------------------------------------------
                        endsubroutine upsend_loc_g23

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine upsend_loc_g34                       &
!***********************************************************************
!                                                                      !
!         Upsend data from generation three to generation four         !
!                                                                      !
!                       - offset version -                             !
!                                                                      !
!***********************************************************************
(V,H,km_64,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: lsendup_sw_loc,lsendup_se_loc                  &
                        ,lsendup_nw_loc,lsendup_ne_loc                  &
                        ,Fitargup_loc34                                &
                        ,itargdn_sw_loc43,itargdn_se_loc43              &
                        ,itargdn_ne_loc43,itargdn_nw_loc43    
use mpi

implicit none

!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_64,flag
real(r_kind), dimension(km_64,1:imL,1:jmL),intent(in):: V
real(r_kind), dimension(km_64,1-hx:im+hx,1-hy:jm+hy),intent(out):: H

!-----------------------------------------------------------------------
real(r_kind), allocatable, dimension(:,:,:)::                            &
                                         sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE &
                                        ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j
integer(i_kind) isend,irecv,nebpe

integer(i_kind):: mygen_dn,mygen_up
integer(i_kind):: itarg_up
integer(i_kind):: itargdn_sw,itargdn_se,itargdn_ne,itargdn_nw
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne

!-----------------------------------------------------------------------
   mygen_dn=3
   mygen_up=4
!
! Define generational flags
!

       itarg_up=Fitargup_loc34(flag)

       itargdn_sw = itargdn_sw_loc43
       itargdn_se = itargdn_se_loc43                  
       itargdn_ne = itargdn_ne_loc43
       itargdn_nw = itargdn_nw_loc43

       lsendup_sw = lsendup_sw_loc
       lsendup_se = lsendup_se_loc
       lsendup_nw = lsendup_nw_loc
       lsendup_ne = lsendup_ne_loc
!-----------------------------------------------------------------------

!N   if(my_hgen==mygen_up) then
      H(:,:,:) = 0.0d0
!N   endif

     ndata =km_64*imL*jmL

!
! --- Send data to SW portion of processors at higher generation
!

      if(  lsendup_sw ) then

        nebpe = itarg_up
    
        if(nebpe == mype) then
           
             do j=1,jmL
             do i=1,imL
                dBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SW(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                       mpi_comm_comp, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )

        deallocate( sBuf_SW, stat = ierr )

        endif

      endif
!
! --- Receive SW portion of data at higher generation
!

      if( itargdn_sw >= 0 ) then

        nebpe = itargdn_sw

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_SW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(1), irecv)
          call MPI_WAIT( rHandle(1), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
                H(:,i,j)=dBuf_SW(:,i,j)
             enddo
             enddo

      endif

!
! --- Send data to SE portion of processors at higher generation
!

      if( lsendup_se ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
                dBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_SE(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype, &
                       mpi_comm_comp, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )

        deallocate( sBuf_SE, stat = ierr )

        endif

      end if

!
! --- Receive SE portion of data at higher generation
!

      if( itargdn_se >= 0 ) then

        nebpe = itargdn_se

        if(nebpe /= mype) then

          call MPI_IRECV( dBuf_SE, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(2), irecv)
          call MPI_WAIT( rHandle(2), istat, ierr )

        endif
             do j=1,jmL
             do i=1,imL
               H(:,imL+i,j)=dBuf_SE(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NW portion of processors at higher generation
!

      if( lsendup_nw ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NW(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NW(:,i,j) = V(:,i,j)
             enddo
             enddo

         call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_comp, sHandle(3), isend)

         call MPI_WAIT( sHandle(3), istat, ierr )

         deallocate( sBuf_NW, stat = ierr )

      end if

    end if

!
! --- Receive NW portion of data at higher generation
!

!      if( my_hgen==mygen_up .and. itargdn_nw >= 0 ) then
      if( itargdn_nw >= 0 ) then

        nebpe = itargdn_nw
 
        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NW, ndata, dtype, nebpe, nebpe,          &
                          mpi_comm_comp, rHandle(3), irecv)
          call MPI_WAIT( rHandle(3), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,i,jmL+j)=dBuf_NW(:,i,j)
             enddo
             enddo

      endif
!
! --- Send data to NE portion of processors at higher generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        if(nebpe == mype) then

             do j=1,jmL
             do i=1,imL
               dBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        else

        allocate( sBuf_NE(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_NE(:,i,j) = V(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype, &
                      mpi_comm_comp, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )

         deallocate( sBuf_NE, stat = ierr )

        endif

      end if

!
! --- Receive NE portion of data at higher generation
!

!N      if( my_hgen==mygen_up .and. itargdn_ne >= 0 ) then
      if( itargdn_ne >= 0 ) then

        nebpe = itargdn_ne

        if(nebpe /= mype) then
          call MPI_IRECV( dBuf_NE, ndata, dtype, nebpe, nebpe,          &
                         mpi_comm_comp, rHandle(4), irecv)
          call MPI_WAIT( rHandle(4), istat, ierr )
        endif

             do j=1,jmL
             do i=1,imL
               H(:,imL+i,jmL+j)=dBuf_NE(:,i,j)
             enddo
             enddo

      endif


!-----------------------------------------------------------------------
                        endsubroutine upsend_loc_g34

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_loc_g43                     &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(W,Z,km_64,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: Fitargup_loc34                                  &
                    ,itargdn_sw_loc43,itargdn_se_loc43                   &
                    ,itargdn_ne_loc43,itargdn_nw_loc43                   &
                    ,lsendup_sw_loc,lsendup_se_loc                       &
                    ,lsendup_nw_loc,lsendup_ne_loc
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_64,flag
real(r_kind), dimension(km_64,1:im,1:jm),intent(in):: W
real(r_kind), dimension(km_64,1:imL,1:jmL),intent(out):: Z
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_64,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe
integer(i_kind) itargdn_sw,itargdn_se,itargdn_nw,itargdn_ne

integer(i_kind):: itarg_up                                           
!-----------------------------------------------------------------------

     Z(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc34(flag)

       itargdn_sw = itargdn_sw_loc43
       itargdn_se = itargdn_se_loc43
       itargdn_nw = itargdn_nw_loc43
       itargdn_ne = itargdn_ne_loc43

       ndata =km_64*imL*jmL

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation

     if(itargdn_sw >= 0) then

        nebpe = itargdn_sw


        allocate( sBuf_SW(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = W(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )

     endif

!
! --- Receive SW portion of data at lower generation


      if( lsendup_sw ) then

        nebpe = itarg_up


        allocate( rBuf_SW(1:km_64,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

      endif

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

     if(itargdn_se >= 0) then

        nebpe = itargdn_se

        allocate( sBuf_SE(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = W(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )

     endif
!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( rBuf_SE(1:km_64,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
     end if

!
! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

     if(itargdn_nw >= 0) then

        nebpe = itargdn_nw


        allocate( sBuf_NW(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NW(:,i,j) = W(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )

     endif

!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        allocate( rBuf_NW(1:km_64,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)


      end if


! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

     if(itargdn_ne >= 0) then

        nebpe = itargdn_ne

        allocate( sBuf_NE(1:km_64,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NE(:,i,j) = W(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

     endif

!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( rBuf_NE(1:km_64,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=1,jmL
             do i=1,imL
               Z(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

      end if

!-----------------------------------------------------------------------
                        endsubroutine downsend_loc_g43

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_loc_g32                     &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(Z,H,km_16,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: Fitargup_loc23                                  &
                    ,itargdn_sw_loc32,itargdn_se_loc32                   &
                    ,itargdn_ne_loc32,itargdn_nw_loc32                   &
                    ,lsendup_sw_loc,lsendup_se_loc                       &
                    ,lsendup_nw_loc,lsendup_ne_loc
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_16,flag
real(r_kind), dimension(km_16,1:im,1:jm),intent(in):: Z
real(r_kind), dimension(km_16,1:imL,1:jmL),intent(out):: H
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_16,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe
integer(i_kind) itargdn_sw,itargdn_se,itargdn_nw,itargdn_ne

integer(i_kind):: itarg_up                                           
!-----------------------------------------------------------------------

     H(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc23(flag)

       itargdn_sw = itargdn_sw_loc32
       itargdn_se = itargdn_se_loc32
       itargdn_nw = itargdn_nw_loc32
       itargdn_ne = itargdn_ne_loc32

       ndata =km_16*imL*jmL

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation

 
  if( itargdn_sw >= 0 ) then

        nebpe = itargdn_sw


        allocate( sBuf_SW(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = Z(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )

  endif

!
! --- Receive SW portion of data at lower generation


      if( lsendup_sw ) then

        nebpe = itarg_up


        allocate( rBuf_SW(1:km_16,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=1,jmL
             do i=1,imL
               H(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

      endif

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if( itargdn_se >= 0 ) then

        nebpe = itargdn_se

        allocate( sBuf_SE(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = Z(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )


  endif
!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( rBuf_SE(1:km_16,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=1,jmL
             do i=1,imL
               H(:,i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
     end if

!
! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

  if( itargdn_nw >= 0 ) then

        nebpe = itargdn_nw


        allocate( sBuf_NW(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NW(:,i,j) = Z(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )


  endif
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        allocate( rBuf_NW(1:km_16,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=1,jmL
             do i=1,imL
               H(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)


      end if


! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if( itargdn_ne >= 0 ) then
        nebpe = itargdn_ne


        allocate( sBuf_NE(1:km_16,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NE(:,i,j) = Z(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )

  endif
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then
        nebpe = itarg_up

        allocate( rBuf_NE(1:km_16,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=1,jmL
             do i=1,imL
               H(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

      end if

!-----------------------------------------------------------------------
                        endsubroutine downsend_loc_g32

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        subroutine downsend_loc_g21                     &
!***********************************************************************
!                                                                      *
!         Downsending data from low resolution pes    (mygen_up)       *
!         to the concurent high-resolution pes        (mygen_dn)       *
!         and add the existing and the recevied values                 *
!                                                                      !
!                       - offset version -                             !
!                                                                      *
!***********************************************************************
(H,V,km_4,flag)
!-----------------------------------------------------------------------
use mg_parameter, only: im,jm,imL,jmL,hx,hy
use mg_domain_loc, only: Fitargup_loc12                                  &
                    ,itargdn_sw_loc21,itargdn_se_loc21                   &
                    ,itargdn_nw_loc21,itargdn_ne_loc21                   &
                    ,lsendup_sw_loc,lsendup_se_loc                       &
                    ,lsendup_nw_loc,lsendup_ne_loc
use mpi

implicit none
!-----------------------------------------------------------------------

integer(i_kind), intent(in):: km_4,flag
real(r_kind), dimension(km_4,1:im,1:jm),intent(in):: H
real(r_kind), dimension(km_4,1:imL,1:jmL),intent(out):: V
!-----------------------------------------------------------------------

real(r_kind), allocatable, dimension(:,:,:)::                            &
                            sBuf_SW,sBuf_SE,sBuf_NW,sBuf_NE              &
                           ,rBuf_SW,rBuf_SE,rBuf_NW,rBuf_NE              

real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_SW
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_SE
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_NW
real(r_kind),dimension(1:km_4,1:imL,1:jmL):: dBuf_NE

integer(i_kind) sHandle(4),rHandle(4),ISTAT(MPI_STATUS_SIZE)
integer(i_kind) iaerr,ierr,iderr,ndata,i,j,L
integer(i_kind) isend,irecv,nebpe
integer(i_kind) itargdn_sw,itargdn_se,itargdn_nw,itargdn_ne
logical:: lsendup_sw,lsendup_se,lsendup_nw,lsendup_ne

integer(i_kind):: itarg_up                                           
!-----------------------------------------------------------------------

     V(:,:,:) = 0.0d0
!
! Define generational flags
!

       itarg_up=Fitargup_loc12(flag)

       itargdn_sw = itargdn_sw_loc21
       itargdn_se = itargdn_se_loc21
       itargdn_nw = itargdn_nw_loc21
       itargdn_ne = itargdn_ne_loc21

       ndata =km_4*imL*jmL

!
! --- Send data from SW portion of processors at the higher generation
!     to corresponding  PE's at lower generation
  
 
  if( itargdn_sw >= 0 ) then
        nebpe = itargdn_sw


        allocate( sBuf_SW(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_SW(:,i,j) = H(:,i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(1), isend)
        call MPI_WAIT( sHandle(1), istat, ierr )
        deallocate( sBuf_SW, stat = ierr )
  
  endif

!
! --- Receive SW portion of data at lower generation
!


      if( lsendup_sw ) then

        nebpe = itarg_up


        allocate( rBuf_SW(1:km_4,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SW, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(1), irecv)
        call MPI_WAIT( rHandle(1), istat, ierr )

             do j=1,jmL
             do i=1,imL
               V(:,i,j)=rBuf_SW(:,i,j)  
             enddo
             enddo

        deallocate( rBuf_SW, stat = iderr)

      endif

!
! --- Send data from SE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if( itargdn_se >= 0 ) then
        nebpe = itargdn_se

        allocate( sBuf_SE(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
               sBuf_SE(:,i,j) = H(:,imL+i,j)
             enddo
             enddo

        call MPI_ISEND( sBuf_SE, ndata, dtype, nebpe, mype,  &
                       mpi_comm_work, sHandle(2), isend)
        call MPI_WAIT( sHandle(2), istat, ierr )
        deallocate( sBuf_SE, stat = ierr )


  endif
!
! --- Receive SE portion of data at lower generation

 
      if( lsendup_se ) then
        nebpe = itarg_up


        allocate( rBuf_SE(1:km_4,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_SE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(2), irecv)
        call MPI_WAIT( rHandle(2), istat, ierr )

             do j=1,jmL
             do i=1,imL
               V(:,i,j)=Rbuf_SE(:,i,j)
             enddo
             enddo

       deallocate( rBuf_SE, stat = iderr)
  
     end if

!
! --- Send data from NW portion of processors at the higher generation
!     to corresponding  PE's at lower generantion

  if( itargdn_nw >= 0 ) then

        nebpe = itargdn_nw


        allocate( sBuf_NW(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NW(:,i,j) = H(:,i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NW, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(3), isend)
        call MPI_WAIT( sHandle(3), istat, ierr )
        deallocate( sBuf_NW, stat = ierr )


  endif
!
! --- Receive NW portion of data at lower generation


      if( lsendup_nw ) then

        nebpe = itarg_up

        allocate( rBuf_NW(1:km_4,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NW, ndata, dtype, nebpe, nebpe, &
                       mpi_comm_work, rHandle(3), irecv)
        call MPI_WAIT( rHandle(3), istat, ierr )

             do j=1,jmL
             do i=1,imL
               V(:,i,j)=Rbuf_NW(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NW, stat = iderr)


      end if


! --- Send data from NE portion of processors at the higher generation
!     to corresponding  PE's at lower generation

  if( itargdn_ne >= 0 ) then

        nebpe = itargdn_ne


        allocate( sBuf_NE(1:km_4,1:imL,1:jmL), stat = iaerr )

             do j=1,jmL
             do i=1,imL
                sBuf_NE(:,i,j) = H(:,imL+i,jmL+j)
             enddo
             enddo

        call MPI_ISEND( sBuf_NE, ndata, dtype, nebpe, mype,  &
                        mpi_comm_work, sHandle(4), isend)
        call MPI_WAIT( sHandle(4), istat, ierr )
        deallocate( sBuf_NE, stat = ierr )


  endif
!
! --- Receive NE portion of data at lower generation
!

      if( lsendup_ne ) then

        nebpe = itarg_up

        allocate( rBuf_NE(1:km_4,1:imL,1:jmL), stat = iaerr )

        call MPI_IRECV( rBuf_NE, ndata, dtype, nebpe, nebpe, &
                        mpi_comm_work, rHandle(4), irecv)
        call MPI_WAIT( rHandle(4), istat, ierr )

             do j=1,jmL
             do i=1,imL
               V(:,i,j)=rBuf_NE(:,i,j)
             enddo
             enddo

        deallocate( rBuf_NE, stat = iderr)

      end if

!-----------------------------------------------------------------------
                        endsubroutine downsend_loc_g21

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        endmodule mg_bocos
