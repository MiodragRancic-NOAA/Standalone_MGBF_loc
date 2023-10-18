!!,nx&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                        program convert
!***********************************************************************
!                                                                      !
!   Convert small arrays into a global array for generation 1          !
!                                                                      !
!***********************************************************************
use def_kind
implicit none
integer(kint), parameter:: nm=72,mm=48
integer(kint), parameter:: nxm=20,mym=20
integer(kint), parameter:: LM=65
integer(kint), parameter:: nm0=nxm*nm,mm0=mym*mm,npem1=nxm*mym-1
integer(kint), parameter:: n0=1,m0=1

real(kfpt),dimension(1:nm,1:mm,0:npem1):: v1
real(kfpt),dimension(1:nm0,1:mm0):: v1_glb
character(len=3):: c_mype
character(len=2):: c_mylev
character(len=11) ime
character(len=11) ime_glb
integer:: mype,n_ind,m_ind,mype_loc,mylev
integer(kint):: L,n,m,nx,my,nglob
integer(kint):: n0_min,n0_max,m0_min,m0_max

real(kfpt),dimension(1:nm,1:Lm,nxm):: Z
real(kfpt),dimension(1:nm0,1:Lm):: Z_glb

logical:: lvert=.true.     !  false for scale decomposition tests
logical:: lxz =.false.
!----------------------------------------------------------------------
!
! G1
!

   do mype=0,npem1

     do m=1,mm
     do n=1,nm
        read(100+mype,'(f12.6)') v1(n,m,mype)
     end do
     end do

     print *,'read fort.',100+mype

     my = mype/nxm
     nx = mype - my*nxm

       m_ind = my*mm+m0
     do m=1,mm
       n_ind = nx*nm+n0
     do n=1,nm
       v1_glb(n_ind,m_ind)=v1(n,m,mype)
       n_ind = n_ind+1
     end do
       m_ind = m_ind+1
     end do
   
  
   enddo 


!!   n0_min=2*nm+n0
!!   n0_max=5*nm
 
!!   m0_min=2*mm+m0
!!   m0_max=5*mm

   n0_min=4*nm+n0
   n0_max=15*nm

   m0_min=4*mm+m0
   m0_max=15*mm

!
!  For src.decomp
!
   
!   n0_min=0
!   n0_max=nm0
!  
!   m0_min=0
!   m0_max=mm0



     open(11,file='v_xy.dat',form='formatted')

     do m=m0_min,m0_max
     do n=n0_min,n0_max
        write(11,'(f12.6)') v1_glb(n,m)
     end do
     end do

     close(11)

!
! Output 1D crossesction
!
     open(13,file='v_x.dat',form='formatted')

     do n=n0_min,n0_max
        write(13,'(f12.6)') v1_glb(n,(m0_max+m0_min)/2)
     enddo

     close(13)

!TEST

!     stop
 

!
! Put together verical cross sections
!
 if(lvert) then

   do nx=1,nxm

     do L=1,Lm
     do n=1,nm
        read(500+nx,'(f12.6)') Z(n,L,nx)
!if(Z(n,L,nx)>0.) then
!  print *,'n,L,nx,Z(n,L,nx)=',n,L,nx,Z(n,L,nx)
!endif
     end do
     end do

     print *,'read fort.',500+nx

   enddo

   do nx=1,nxm
     do L=1,Lm
     do n=1,nm
       nglob=(nx-1)*nm+n
       Z_glb(nglob,L)=Z(n,L,nx)
!TEST
!if(nx==3.and.n==nm-3.and.Z(n,L,nx)>0.) then
!if(Z(n,L,nx)>0.) then
!  print *,'nglob,L,Z_glb(nglob,L)=',nglob,L,Z_glb(nglob,L)
!endif
!TEST
     enddo
     enddo
   enddo

   if(lxz) then
        
     open(12,file='v_xz.dat',form='formatted')

!     do l=1,lm
!     do n=1,nm0
     do l=15,35
     do n=2*nm+1,5*nm
        write(12,'(f12.6)') Z_glb(n,L)
!TEST
!  if(n==141) then
!     print*,'WRITE: n,L,Z_glb(n,L)=',n,L,Z_glb(n,L)
!  endif
!TEST
     end do
     end do

     close(12)

   endif

     open(18,file='v_z.dat',form='formatted')
       do l=1,lm
!         write(18,'(f12.6)') Z_glb(3*nm+nm/2,l)
         write(18,'(f12.6)') Z_glb((n0_min+n0_max)/2,l)
       enddo
     close(18)

 endif
  

   
!----------------------------------------------------------------------
                        endprogram convert
