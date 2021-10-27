module mesh_grid
!-----------------------------------------------------------------------
! modules
!-----------------------------------------------------------------------
      use grids
      use mesh

contains
!-------------------------------------------------------------------------!
!--------------------- generate mesh info --------------------------------!
!-------------------------------------------------------------------------!

      subroutine init_mesh( grid )

      implicit none

      ! arguments
      type(type_grid), intent(inout), target :: grid

      ! locals
      integer :: nx, ny ! grid dimensions
      real(kind=4), dimension(:,:), pointer :: lon_corn, lat_corn
      logical,      dimension(:,:), pointer :: mask
      type(type_mesh), pointer :: mesh_p

!-----------------------------------------------------------------------
! iniitalize local dimensions
!-----------------------------------------------------------------------

      nx = grid % nx
      ny = grid % ny
      call remove_periodic_points( grid, grid % lld % mask )

      mesh_p => grid % meshd
      nn  => mesh_p % nn
      ne  => mesh_p % ne
nconnect  => mesh_p % nconnect
nconnect=4

!-----------------------------------------------------------------------
!
!     compute connectivity based on mask
!     between structured grid and unstructured mesh
!     and find and remove redundant nodal points from the mesh info
!
!-----------------------------------------------------------------------

      call get_connectivity( mesh_p, grid, .true. )

!-----------------------------------------------------------------------
!
!     construct unstructured mesh pointers
!
!-----------------------------------------------------------------------

!      call pointer_mesh( mesh_p, .true. )
!      call construct_tree_mesh( mesh_p, .true. )

!-----------------------------------------------------------------------
!
!     allocate other grid fields
!
!-----------------------------------------------------------------------

!      allocate( mesh_p % geo_p % area        (ne), &
!     &          mesh_p % geo_p % frac        (ne), &
!     &          mesh_p % geo_p % centroid_lat(ne), &
!     &          mesh_p % geo_p % centroid_lon(ne))
!
!      mesh_p % geo_p % area = zerd
!      mesh_p % geo_p % frac = zerd
!      mesh_p % geo_p % centroid_lat = zerd
!      mesh_p % geo_p % centroid_lon = zerd

end subroutine init_mesh

!***********************************************************************

!-----------------------------------------------------------------------
!
!     this routine computes the
!     mesh connectivity based on mask
!     present_grid already set
!     nn,ne, nconnect already linked
!     assumes that lon/lat corners correspond to NEMO F-point
!
!-----------------------------------------------------------------------
      subroutine get_connectivity( mesh_p, grid, message )

      implicit none
!-----------------------------------------------------------------------
! arguments
!-----------------------------------------------------------------------
      type(type_mesh), target :: mesh_p
      type(type_grid), target :: grid
      logical message

!-----------------------------------------------------------------------
! locals
!-----------------------------------------------------------------------
      integer nx, ny
      real(kind=4), dimension(:,:), pointer :: lon_corn, lat_corn
      logical,      dimension(:,:), pointer :: mask
      integer i, j, ele, node, ik, jk, node2
      integer, pointer, dimension(:,:) :: pstr
      integer, pointer, dimension(:) :: pele
      integer, allocatable :: inele(:,:),pnode(:)

!-----------------------------------------------------------------------
! iniitalize local dimensions
!-----------------------------------------------------------------------

      nx = grid % nx
      ny = grid % ny
      mask => grid % lld % mask
      nn  => mesh_p % nn
      ne  => mesh_p % ne
nconnect  => mesh_p % nconnect
      allocate( inele(nconnect,nx*ny), pnode((nx+1)*(ny+1)) )

!-----------------------------------------------------------------------
!
!     link arrays
!
!-----------------------------------------------------------------------

      mask     => grid % lld % mask
      lon_corn => grid % lld % lonc
      lat_corn => grid % lld % latc

!-----------------------------------------------------------------------
!
!     adjust north pole latitude
!
!-----------------------------------------------------------------------

      if (grid % std_grd % grtyp == 'U' .or. &
          grid % std_grd % grtyp == 'B' .or. &
          grid % std_grd % grtyp == 'G' .or. &
        ( grid % std_grd % grtyp == 'Z'.and. grid % nbdy == 'N') ) then
         lat_corn(:,ny) = 90.
      endif

!-----------------------------------------------------------------------
!
!     adjust southern pole latitude
!
!-----------------------------------------------------------------------

      if (grid % std_grd % grtyp == 'U' .or. &
          grid % std_grd % grtyp == 'B' .or. &
          grid % std_grd % grtyp == 'G') then
         lat_corn(:,0) = -90.
      endif

!-----------------------------------------------------------------------
!
! find list of active element
!
!-----------------------------------------------------------------------

      allocate( mesh_p % pm2ps(nx,ny) )
      pstr   => mesh_p % pm2ps ; pstr(:,:) = 0
      allocate( mesh_p % ps2pm(nx*ny) )
      pele   => mesh_p % ps2pm; pele(:) = 0

      ele=0
      pnode(:)=0
      do j=1,ny
         do i=1,nx
            if (mask(i,j)) then
               ele = ele + 1
               pstr(i,j) = ele
               pele(ele) = nx*(j-1) + i
! go through the 4 vertices of the quadrangle (nconnect=4)
               do jk=0,1
                  do ik=0,1
                     node = (j - 1 + jk) * (nx+1) + i + ik
                     pnode(node) = 1
                     inele(2*jk+ik+1,ele) = node
                   enddo
               enddo
               ! switch order for anticloxkwise convention
               node = inele(3,ele)
               inele(3,ele) = inele(4,ele)
               inele(4,ele) = node
            endif
         enddo
      enddo

      ! pass T-point grid lat/lon to mesh centroids
      ne = ele
      allocate( mesh_p % lonc(ne), mesh_p % latc(ne) )
      ele=0
      do j=1,ny
         do i=1,nx
            if (mask(i,j)) then
               ele = ele + 1
               mesh_p % lonc(ele) = grid % lld % lont(i,j)
               mesh_p % latc(ele) = grid % lld % latt(i,j)
            endif
         enddo
      enddo

!-----------------------------------------------------------------------
!
! identify redundant corners
!
!-----------------------------------------------------------------------

      if ( grid % webdy ) then

         select case( grid % werepeat )

         case (1) ! GEM style redundancy

           ! identify the east point
           i = nx + 1
           do j = 0 , ny
              node = j * (nx+1) + i
              node2 = j * (nx+1) + 1
              if ( pnode(node2) > 0 ) pnode(node) = - node2
           enddo

         case (2) ! NEMO style redundancy

           ! identify the east point (the cell 1 and nx are already duplicate and should be eliminated)
           i = nx
           do j = 0 , ny
              node = j * (nx+1) + i
              node2 = j * (nx+1) + 2
              if ( pnode(node2) > 0 ) pnode(node) = - node2
           enddo

         end select
      endif

      select case ( grid % nbdy )

! will do that when we can have connectivity (i.e. mix of triangle and quadrangle)
!      case('N')
!
!           j = ny
!           do i = 1, nx
!              node = j * (nx+1) + i + 1
!              pnode(node) = - ( j * (nx+1) + 1 )
!           enddo

      case('T')
           j = ny - grid % nrepeat
           do i = 1, nx/2
              node = j * (nx+1) + i + 1
              node2 = (j-1) * (nx+1) + nx - i + 2
              if ( pnode(node2) > 0 ) pnode(node) = - node2
           enddo

      case('F')
           j = ny - grid % nrepeat
           do i = nx/2+1, nx
              node = j * (nx+1) + i + 1
              node2 = j * (nx+1) + nx - i + 1
              if ( pnode(node2) > 0 ) pnode(node) = - node2
           enddo

      end select

!-----------------------------------------------------------------------
!
! reduce the list of nodes
!
!-----------------------------------------------------------------------

      nn = 0
      do j = 0 , ny
         do i = 0 , nx
            node = j * (nx+1) + i + 1
            if ( pnode(node) == 1 ) then
               nn = nn + 1
            endif
         enddo
      enddo

      allocate(mesh_p % longr(nn))
      allocate(mesh_p % latgr(nn))
      longr => mesh_p % longr
      latgr => mesh_p % latgr

      nn = 0
      do j = 0 , ny
         do i = 0 , nx
            node = j * (nx+1) + i + 1
            if ( pnode(node) == 1 ) then
               nn = nn + 1
               pnode(node) = nn
               longr(nn) = lon_corn(i,j)
               latgr(nn) = lat_corn(i,j)
            endif
         enddo
      enddo

!-----------------------------------------------------------------------
!
! fill the final connectivity table
!
!-----------------------------------------------------------------------

      allocate(mesh_p % in(nconnect,ne))
      in => mesh_p % in

      do ele = 1, ne
         do ik = 1, nconnect
            node = pnode(inele(ik ,ele))
            if ( node < 0 ) then ! redundant point, find the associated point
               node = - node
               j = ( node - 1 ) / ( nx + 1 )
               i = node - j * ( nx + 1 ) - 1
               node = j * (nx+1) + i + 1
!               write(*,*) 'find node',i,j,ele, node, pnode(node)
               node = pnode(node)
               if (node==0) stop
            endif
           in(ik,ele) = node
        enddo
      enddo
!      write(*,*) 'min/max in',minval(in),maxval(in)
      deallocate( inele, pnode )

!-----------------------------------------------------------------------
!
! print if needed
!
!-----------------------------------------------------------------------

      if (message) then      
         write(*,fmt='(1x,a10,i7)') 'nn       :',nn
         write(*,fmt='(1x,a10,i7)') 'ne       :',ne
      endif

      end subroutine get_connectivity

!! --------------------------------------------------------------------------------------
!! --------------------------------------------------------------------------------------
      subroutine construct_tree_mesh_center(mesh_p, message)
!! -----------------------------------------------------------------------------------------------
!!  *** construct a kd-tree for the unstructured mesh based on centroids
!! -----------------------------------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
! arguments
!-----------------------------------------------------------------------
      type(type_mesh), target :: mesh_p
      logical message
       
!-----------------------------------------------------------------------
! locals
!-----------------------------------------------------------------------
      integer, pointer :: ne
      REAL(KIND=4) :: start, finish
      INTEGER ntot, ji, jj, i1, i2
      REAL(kind=4), ALLOCATABLE, DIMENSION(:,:) :: data3d
      TYPE (kdtree2), pointer :: my_kdtree
      REAL(KIND=4), DIMENSION(:), POINTER :: lon, lat

      call cpu_time(start)

! create the KD-tree for the source grid
      ne     => mesh_p % ne
      lon    => mesh_p % lonc
      lat    => mesh_p % latc

      allocate( mesh_p % kdtree )
      my_kdtree => mesh_p % kdtree

      allocate(data3d(3,ne))
      call ez_lac(data3d, lon, lat, ne)
      !
      my_kdtree = kdtree2_create(data3d)
      deallocate(data3d)
      my_kdtree % the_data => NULL()

      if (message) &
          write(*,*) 'kdtree',mesh_p % kdtree % n, mesh_p % kdtree % dimen
      call cpu_time(finish)
      if (message) &
          print '("KD-tree creation Time = ",f6.3," seconds.")',finish-start

      end subroutine construct_tree_mesh_center

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module
