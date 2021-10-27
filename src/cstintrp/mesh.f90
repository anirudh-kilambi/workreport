module mesh

use grids

implicit none

 integer, pointer :: nn,ne,nf,nfbd,fmesh,nftrdim,nconnect
 INTEGER, pointer :: in(:,:) => NULL() &
                  ,tseg(:,:) => NULL() &
                  ,pseg(:,:) => NULL() &
                  , itm(:,:) => NULL() &
                  ,  nmit(:) => NULL() &
                  ,  tmit(:) => NULL() &
                  ,  imit(:) => NULL() &
                  ,   sbd(:) => NULL()
 real(kind=8), pointer :: xgr(:) => NULL() &
                         ,ygr(:) => NULL()
 real(kind=4), pointer :: longr(:) => NULL() &
                         ,latgr(:) => NULL()
 
 contains
!
!
! --------------------------------------- 
!    link module defined pointer to mesh pointer arrays
! --------------------------------------- 
!
      subroutine link_mesh( mesh_p )
      implicit none
      ! arguments
      type(type_mesh), target :: mesh_p
      ! locals
      integer i,j,k
      character (len=20) :: format1
      character (len=80) :: nodname,elename
      double precision :: xdum,ydum
      logical message

! start transfer
      nn  => mesh_p % nn
      ne  => mesh_p % ne
      nf  => mesh_p % nf
      nconnect  => mesh_p % nconnect
      in  => mesh_p % in
      xgr => mesh_p % xgr
      ygr => mesh_p % ygr
      tseg=> mesh_p % tseg
      pseg=> mesh_p % pseg
      itm => mesh_p % itm
      nmit=> mesh_p % nmit
      imit=> mesh_p % imit
      tmit=> mesh_p % tmit
      sbd => mesh_p % sbd
      nfbd=> mesh_p % nfbd
      nftrdim => mesh_p % nftrdim
! end transfer

      end subroutine link_mesh
!
! --------------------------------------- 
!    deallocate mesh arrays
! --------------------------------------- 
!
      subroutine deallocate_mesh( mesh_p )
      implicit none
      ! arguments
      type(type_mesh), target :: mesh_p
      ! locals
      integer i,j,k
      character (len=20) :: format1
      character (len=80) :: nodname,elename
      double precision :: xdum,ydum
      logical message

! start transfer
      mesh_p % nn = 0
      mesh_p % ne = 0
      mesh_p % nf = 0
      mesh_p % nconnect = 0
      if( associated( mesh_p % in ) )    nullify( mesh_p % in )
      if( associated( mesh_p % xgr ) )   nullify( mesh_p % xgr )
      if( associated( mesh_p % ygr ) )   nullify( mesh_p % ygr )
      if( associated( mesh_p % longr ) ) nullify( mesh_p % longr )
      if( associated( mesh_p % latgr ) ) nullify( mesh_p % latgr )
      if( associated( mesh_p % tseg ) )  nullify( mesh_p % tseg )
      if( associated( mesh_p % pseg ) )  nullify( mesh_p % pseg )
      if( associated( mesh_p % itm ) )   nullify( mesh_p % itm )
      if( associated( mesh_p % nmit ) )  nullify( mesh_p % nmit )
      if( associated( mesh_p % imit ) )  nullify( mesh_p % imit )
      if( associated( mesh_p % tmit ) )  nullify( mesh_p % tmit )
      if( associated( mesh_p % sbd ) )   nullify( mesh_p % sbd )
      mesh_p % nfbd  = 0
      mesh_p % nftrdim = 10
! end transfer

      end subroutine deallocate_mesh
        
!
! --------------------------------------- 
!                                          
!    compute some pointers between node to cell, face to cell,nodes, etc....                         
!                                          
! --------------------------------------- 
!
      subroutine pointer_mesh( mesh_p, message, debug )
!
      implicit none

      ! arguments
      type(type_mesh), target :: mesh_p

      ! locals
      integer pouet
      integer, allocatable :: ndeb(:),nfin(:)
      logical message
      logical, optional :: debug
!
! entiers de service
      integer i,j,ind,k,tk,ik,lk,pj,it,ne1,ne2,ne3,i1,i2,nc,nd,buf,j1,j2
      logical trouve
      character :: format1*20
!
! associate scalar
      nn      => mesh_p % nn
      ne      => mesh_p % ne
      nf      => mesh_p % nf
      nfbd    => mesh_p % nfbd
      nftrdim => mesh_p % nftrdim      
      nconnect=> mesh_p % nconnect
      nftrdim = 10

! allocate grid arrays
      allocate(mesh_p % itm(nconnect,ne))
      allocate(mesh_p % tseg(2,nconnect*nn))
      allocate(mesh_p % pseg(2,nconnect*nn))
      allocate(mesh_p % nmit(nn+1))
      allocate(mesh_p % imit(nn*nftrdim))
      allocate(mesh_p % tmit(nn*nftrdim))
      allocate(ndeb(nn+1))
      allocate(nfin(nn+1))
! associate arrays
      in  => mesh_p % in
      tseg=> mesh_p % tseg
      pseg=> mesh_p % pseg
      itm => mesh_p % itm
      nmit=> mesh_p % nmit
      imit=> mesh_p % imit
      tmit=> mesh_p % tmit

if (message) &
      write(*,*) 'compute mesh-related pointers'
!
      tseg=0
!
!
! triangles (elements) mitoyens d'un noeud
! on suppose que chaque noeud est entoure de nftrdim quadrilateres
! et on remplit les vecteurs tmit et imit
!
if (present(debug)) then
write(*,*) nn,ne,nconnect
do i=1,ne
  write(*,*) in(:,i)
enddo
endif
      ndeb(1)=1
      nfin(1)=1
      do i=1,nn-1
        ndeb(i+1)=ndeb(i)+nftrdim
        nfin(i+1)=ndeb(i+1)
      enddo
!
        do i=1,ne
           do k=1,nconnect
              ne1=in(k,i)
              ind=nfin(ne1)
              tmit(ind)=i
              imit(ind)=k
              nfin(ne1)=ind+1
              if (ne1.lt.nn) then
              do while (nfin(ne1).eq.ndeb(ne1+1))
                 ne1=ne1+1
                 do j1=nfin(ne1)-1,ndeb(ne1),-1
                    tmit(j1+1)=tmit(j1)
                    imit(j1+1)=imit(j1)
                 enddo
                 ndeb(ne1)=ndeb(ne1)+1
                 nfin(ne1)=nfin(ne1)+1
              enddo
              endif
           enddo
         enddo
!
! elimination des vides laisses dans tmit et imit
!
         ind=0
         do i=1,nn
            do j=ndeb(i),nfin(i)-1
               ind=ind+1
               tmit(ind)=tmit(j)
               imit(ind)=imit(j)
            enddo
        enddo
        nmit(1)=1
        do i=1,nn
           nmit(i+1)=nmit(i)+nfin(i)-ndeb(i)
        enddo
if (present(debug)) then
do i=1,nn
  write(*,*) 'tmit',tmit(nmit(i):nmit(i+1)-1)
  write(*,*) 'imit',imit(nmit(i):nmit(i+1)-1)
enddo
endif
!
!
!
      nf=0
      do i=1,ne
         do j=1,nconnect
            ne1=in(j,i)
            ne2=in(mod(j,nconnect)+1,i)
!
!   checker les anciennes faces des elements mitoyens
!
            k=nmit(ne2)
            tk=tmit(k)
            trouve=.true.
            ik=imit(k)            
            do while(tk.lt.i.and.trouve)
               ne3=in(mod(ik,nconnect)+1,tk)
!if (present(debug)) then
!               write(*,*) 'tk',i,tk,ne1,ne2,ne3
!endif
               if (ne1.eq.ne3) then
                  lk=itm(ik,tk)
                  tseg(2,lk)=i
                  itm(j,i)=lk
                  trouve=.false.
               endif
               k=k+1
               tk=tmit(k)
               ik=imit(k)            
            enddo
!
! la face est nouvelle!
!
            if (trouve) then
               nf=nf+1
               itm(j,i)=nf
               tseg(1,nf)=i
               pseg(1,nf)=ne1
               pseg(2,nf)=ne2
            endif
         enddo      
      enddo

if (present(debug)) then
do i=1,nf
  write(*,*) tseg(:,i)
enddo
endif
      deallocate (ndeb)
      deallocate (nfin)


! --------- creation of the boundary pointer

      nfbd=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            nfbd=nfbd+1
         endif
      enddo
      allocate(mesh_p % sbd(nfbd))
      sbd => mesh_p % sbd
      k=0
      do i=1,nf
         if (tseg(2,i).eq.0) then
            k=k+1
            sbd(k)=i
         endif
      enddo

if (message) then      
 format1='(1x,a10,i7)'
      write(*,fmt=format1) 'faces    :',nf
endif

      return
      end subroutine pointer_mesh

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


      SUBROUTINE ARC(XX1,YY1,XX2,YY2,ARCL)
!
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

!----------------------------------------------------------------------
!     solve the arc length through the earth center ("grand arc")
!
!     INPUT : (XX1,YY1) position in E-longitude, N-latitude (degrees)
!             (XX2,YY2) position in E-longitude, N-latitude (degrees)
!     OUTPUT: ARCL "grand arc" length (metres)
!
!----------------------------------------------------------------------

!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
! arguments
!----------------------------------------------------------------------
      double precision XX1,YY1,XX2,YY2 
!----------------------------------------------------------------------
! locals
!----------------------------------------------------------------------
      double precision xa,ya,za,xb,yb,zb,ab,aob,ARCL
      double precision, parameter :: rearth = 6371.d3
!----------------------------------------------------------------------

      xa=cos(yy1)*cos(xx1)
      ya=cos(yy1)*sin(xx1)
      za=sin(yy1)

      xb=cos(yy2)*cos(xx2)
      yb=cos(yy2)*sin(xx2)
      zb=sin(yy2)

      AB=SQRT((xb-xa)**2+(yb-ya)**2+(zb-za)**2)
      AOB=2.d0*asin(AB*0.5d0)
      ARCL=rearth*AOB

      END subroutine arc

! --------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------
  subroutine get_latlon_xyz_centroid( mesh_p )
    implicit none
    ! arguments
    type(type_mesh), target :: mesh_p
    ! locals
    real(kind=4), dimension(:,:), allocatable :: c3d, n3d
    integer nn, ne, nc
    integer j, inode, icell

    nn = mesh_p % nn; ne = mesh_p % ne; nc = mesh_p % nconnect
    if (.not.associated( mesh_p % lonc )) then
       allocate( mesh_p % lonc(ne) )
       allocate( mesh_p % latc(ne) )
    endif

    allocate(c3d(3,ne), n3d(3,nn))
       call ez_lac(n3d, mesh_p % longr, mesh_p % latgr, nn)
    c3d(:,:) = zerd
    do icell = 1, ne
       do j = 1, nc
          inode = mesh_p % in(j,icell)
          c3d(:, icell) = c3d(:, icell) + n3d(:,inode)
       enddo
    enddo
    c3d(:,:) = c3d(:,:) / real(nc)
    call ez_cal( mesh_p % lonc, mesh_p % latc, c3d, ne)

    deallocate(c3d,n3d)

    end subroutine get_latlon_xyz_centroid

end module mesh
