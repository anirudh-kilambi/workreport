MODULE grids

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE grids ***
  !!
  !! history :
  !!     Original code : F. Roy (November 2009)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!
  USE kdtree2_module
  USE std
  USE mympp

  IMPLICIT NONE

  ! parameters
  real(kind=8), parameter :: zerd = 0.0_8, oned = 1.0_8, epsilon = 1e-15_8

!********************************************************************
!
! structured weight
!
!********************************************************************

TYPE type_weights
 integer :: nx, ny
 integer :: numwgt                                                ! can be the max number of weights across all points
 real(kind=8), pointer, dimension(:,:,:) ::  wgt => NULL()        ! structured weights
 integer, pointer, dimension(:,:,:)      :: iwgt => NULL() &
                                          , jwgt => NULL()        ! structured indices of weights
 integer, pointer, dimension(:,:)        :: nwgt => NULL()        ! number of weights/points for each destination point
 logical, pointer, dimension(:,:)        :: mwgt => NULL()        ! mask for the interpolation
END TYPE type_weights

!********************************************************************
!
! structured data
!
!********************************************************************

TYPE type_data
 logical     , pointer, dimension(:,:)   ::  m => NULL()   ! mask for cells
 real(kind=8), pointer, dimension(:,:)   ::  u => NULL() & ! main field or first vector component
                                          ,  v => NULL()   ! second vector component
END TYPE type_data

!********************************************************************
!
! structured location arrays
!
!********************************************************************

TYPE type_location
 real(kind=8), pointer, dimension(:,:)   :: ag => NULL() &  ! fractional position along the i-direction of a cell
                                          , bg => NULL()    ! fractional position along the j-direction of a cell
 real(kind=4), pointer, dimension(:,:)   :: pg => NULL()    ! differential rotation angle
 integer, pointer, dimension(:,:)        :: ig => NULL() &  ! i-index along the grid
                                          , jg => NULL()    ! j-index along the grid
 logical, pointer, dimension(:,:)        :: mg => NULL()    ! points located outside the grid are tagged as false
END TYPE type_location

!********************************************************************
!
! latlon structured
!
!********************************************************************

TYPE type_lonlat
 logical     , pointer ::     mask(:,:) => NULL()   ! mask for cells
 real(kind=4), pointer ::     lonc(:,:) => NULL() & ! corner points
                             ,latc(:,:) => NULL() &
                             ,lont(:,:) => NULL() & ! centered T-point
                             ,latt(:,:) => NULL() &
                             ,vect(:,:) => NULL() & ! dummy vector
                             ,tic (:)   => NULL() & 
                             ,tac (:)   => NULL() &
                             ,t3d(:,:,:) => NULL() ! 3d cartesian coordinates at centered T-points
END TYPE type_lonlat

!********************************************************************
!
! Geometry declaration
!
!********************************************************************

TYPE type_geometry

 real(kind=8),pointer, dimension (:)  ::  &
      centroid_lon => NULL() &
     ,centroid_lat => NULL() &
     ,area  => NULL() &
     ,frac  => NULL() &
     ,art1  => NULL() &
     ,art2  => NULL() &
     ,tanearth => NULL() &
     , coscell => NULL()

END TYPE type_geometry

TYPE type_weight_unstr

  integer :: num_links_map, max_num_links
  integer, pointer, dimension (:)  :: src_cell, dst_cell
  real(kind=8), pointer, dimension (:,:)  :: wts_map
  integer, pointer, dimension(:,:) :: link_add

  logical :: alloc_wgts = .false.

END TYPE type_weight_unstr

!********************************************************************
!
! Unstructured Mesh declaration
!
!********************************************************************

TYPE type_mesh
 integer :: nn,ne,nf,nfbd,fmesh,nconnect
 integer :: nftrdim=15
 INTEGER, pointer :: in(:,:) => NULL() &
                  ,tseg(:,:) => NULL() &
                  ,pseg(:,:) => NULL() &
                  , itm(:,:) => NULL() &
                  ,  nmit(:) => NULL() &
                  ,  tmit(:) => NULL() &
                  ,  imit(:) => NULL() &
                  ,   sbd(:) => NULL()
 real(kind=4), pointer :: longr(:) => NULL() &
                         ,latgr(:) => NULL() &
                          ,lonc(:) => NULL() &
                          ,latc(:) => NULL()
 real(kind=8), pointer ::   xgr(:) => NULL() &
                           ,ygr(:) => NULL() &
                            ,xc(:) => NULL() &
                            ,yc(:) => NULL()
 integer, pointer ::         pm2ps(:,:) => NULL()
 integer, pointer ::         ps2pm(:)   => NULL()
 real(kind=8) :: scx,scy
 TYPE (kdtree2), pointer :: kdtree => NULL()
 TYPE(type_geometry)     :: geo_p
END TYPE type_mesh


!********************************************************************
!
! grid data type
!
!********************************************************************

TYPE type_grid
 integer :: nx, ny
 integer           :: level
 character(len=20) :: ctype
 ! management of periodicity
 INTEGER werepeat, nrepeat
 LOGICAL webdy
 CHARACTER(LEN=1) :: nbdy = ''
 TYPE(type_lonlat)   :: lld
 TYPE(type_data)     :: vard
 TYPE(type_location) :: locd
 TYPE(type_weights)  :: wgtd
 TYPE (kdtree2), pointer :: kdtree => NULL()
 integer dxtree, dytree   ! stride in creating the tree (the more structured is the grid, the larger these values can be)
 integer nxred, nyred     ! final size of the reduced grid
 integer, pointer, dimension(:)          :: ired => NULL() &         ! pointer to reduced index on kd-tree
                                          , jred => NULL()
 TYPE (std_grid) :: std_grd ! regular grid std descriptor
 TYPE (std_grid) :: aux_grd ! auxiliary grid std for axis descriptor
 integer         :: ng
 TYPE (type_grid), ALLOCATABLE, DIMENSION(:) :: subgrd
 TYPE (type_mesh) :: meshd
END TYPE type_grid

!********************************************************************
!
! variable meta data
!
!********************************************************************

TYPE variable
   TYPE(std_variable) :: varstd
   ! additional info for mask
   logical :: ln_mask = .false.
   TYPE(std_variable) :: mask
   ! additional info for vector
   logical :: ln_vector = .false.
   TYPE(std_variable) :: varstdv
END TYPE variable


PUBLIC :: init_grid_parameter, init_sub_grid, get_sub_axes, getll, allocate_grid_intrp, construct_tree, &
          get_delta_i_xyz, get_dlatlon_xyz


CONTAINS

  SUBROUTINE init_grid_parameter( grid, unitf )
    !! -----------------------------------------------------------------------------------------------
    !!  *** Initialize cstintrp options, dimensions, parameters, validate with 
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE

    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    integer             :: unitf ! RPN file unit

    ! * Local variables
    type (std_grid), pointer &
                        :: grd   ! grid type 
    integer             :: ierr
    integer             :: nx, ny, gid, ig
    character(len=1)    :: gridtype
    integer, pointer    :: ng
    ! librmn functions
    integer, external    :: ezqkdef, ezget_nsubgrids


    grd => grid % std_grd
    ng  => grid % ng
    nx = grd % ni
    ny = grd % nj

    gid = ezqkdef(nx, ny, grd % grtyp, grd % ig1, grd % ig2, grd % ig3, grd % ig4, unitf) ! defines the grid
    grd % grid_id = gid
    write(*,*) 'grid info',grd

    ng=1

    select case ( grd % grtyp )
    case ('U')
      ng = ezget_nsubgrids(gid)
      if ( ng /= 2 ) then
         write(*,*) 'CSTINTRP: no support for such U grid',grd % grid_id
         write(*,*) 'CSTINTRP: number of subgrids',ng
         call my_abort
      endif
    case default
    end select

    allocate( grid % subgrd(ng) )

  END SUBROUTINE init_grid_parameter

  SUBROUTINE init_sub_grid( grid )
    !! -----------------------------------------------------------------------------------------------
    !!  *** Initialize subgrids if U super grid or just pass info to subgrid
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE

    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 

    ! * Local variables
    type (std_grid), pointer &
                        :: grd   ! grid type 
    type (type_grid), dimension(:), pointer &
                        :: subg  ! grid type 
    integer, pointer    :: ng    ! number of subgrid
    integer             :: nx, ny
    integer             :: ierr
    integer             :: gid, ig
    integer, allocatable, dimension(:) :: &
                           subgid
    type(std_grid), pointer :: &
                           gr, g1 ! pointer grid for subgrid and its auxiliary descriptor
   ! librmn functions
    integer, external   :: ezget_subgridids, ezgxprm

    grd  => grid % std_grd
    subg => grid % subgrd
    ng   => grid % ng
    nx = grd % ni
    ny = grd % nj
    grid % nx = nx
    grid % ny = ny

    select case ( ng )
    case(1)
      subg(1) % std_grd = grd
      subg(1) % nx      = grid % nx
      subg(1) % ny      = grid % ny
      subg(1) % dxtree  = grid % dxtree
      subg(1) % dytree  = grid % dytree
      subg(1) % nbdy    = grid % nbdy
      subg(1) % nrepeat = grid % nrepeat
      subg(1) % webdy   = grid % webdy
      subg(1) % werepeat= grid % werepeat
    case(2)
      ! get info about all subgrida
      allocate(subgid(ng))
      gid = grd % grid_id
      ierr = ezget_subgridids(gid,subgid)
      do ig = 1, ng
        gr => subg(ig) % std_grd
        g1 => subg(ig) % aux_grd
        gid = subgid(ig)
        gr % grid_id = gid
        ierr = ezgxprm(gid, gr % ni, gr % nj, gr % grtyp, gr % ig1, gr % ig2, gr % ig3, gr % ig4, &
                                              g1 % grtyp, g1 % ig1, g1 % ig2, g1 % ig3, g1 % ig4) ! returns the grid info and auxiliary grid info
        g1 % ni = gr % ni
        g1 % nj = gr % nj
        if ( gr % grtyp /= 'Z' )  then; print *,'CSTINTRP: no support for such grid type subgrid', gid         ; call my_abort; endif
        if ( gr % ni /= nx   )    then; print *,'CSTINTRP: no support for such grid type subgrid dimensions NI'; call my_abort; endif
        if ( gr % nj /= ny/2 )    then; print *,'CSTINTRP: no support for such grid type subgrid dimensions NJ'; call my_abort; endif
        subg(ig) % nx = gr % ni
        subg(ig) % ny = gr % nj
      enddo
      deallocate(subgid)
    end select

    ! get 1d axes if RPN grids
    call get_sub_axes( grid )

  END SUBROUTINE init_sub_grid

  SUBROUTINE get_sub_axes( grid )
    !! -----------------------------------------------------------------------------------------------
    !!  *** get 1D axes from subgrid
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    ! locals
    type (std_grid), pointer &
                        :: grd   ! grid type 
    type (type_grid), dimension(:), pointer &
                        :: subg  ! grid type 
    integer, pointer    :: ng    ! number of subgrid
    real(kind=4), pointer, dimension(:) &
                        :: tic, tac
    integer             :: nx, ny, ig, ierr, gid
   ! librmn functions
    integer, external   :: gdgaxes

    grd  => grid % std_grd
    subg => grid % subgrd
    ng   => grid % ng

    select case( grd % grtyp )
    case( 'Z','U' )

      do ig = 1, ng
        nx = subg(ig) % nx
        ny = subg(ig) % ny
        allocate(subg(ig) % lld % tic(nx),subg(ig) % lld % tac(ny))
        tic => subg(ig) % lld % tic
        tac => subg(ig) % lld % tac
        gid = subg(ig) % std_grd % grid_id
        ierr = gdgaxes( gid, tic, tac)
      enddo
    case default
      write(*,*) 'no sub-axes to be found'
    end select

  END SUBROUTINE get_sub_axes

  SUBROUTINE getll(grid, unitf)
    !! -----------------------------------------------------------------------------------------------
    !!  *** read the lat/lon
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    integer             :: unitf ! RPN file unit

    ! * Local variables
    type (std_grid), pointer &
                        :: grd   ! grid type 
    integer, pointer    :: ng    ! number of subgrid
    type (type_grid), dimension(:), pointer &
                        :: subg  ! grid type 
    integer             :: nx, ny, ig, ierr, status_st90
    integer             :: gid, nv
    type(type_grid), pointer &
                        :: gr ! pointer grid
    real(kind=4), pointer, dimension(:,:) &
                        :: wrkr4
    integer nlist, liste(10), l, nk
   ! librmn functions
    integer, external   :: gdll, fstinl, fstluk

    grd  => grid % std_grd
    subg => grid % subgrd
    ng   => grid % ng

    ! FD not sure this is needed
    if ( ( grd % grtyp == 'X' .or. grd % grtyp == 'O') .and. ng > 1 ) then
       write(*,*) 'CSTINTRP: NEMO SUPER GRID ... IMPOSSIBLE'
       call my_abort
    endif

    do ig = 1, ng

      gr => subg(ig)
      nx = gr % nx
      ny = gr % ny
      allocate( gr % lld % latt(nx,ny), gr % lld % lont(nx,ny) )

      select case ( grd % grtyp )
      case('Z','B','U','N','A','S','Y','L','G')

        gid = gr % std_grd % grid_id
        ierr = gdll(gid, gr % lld % latt, gr % lld % lont)
        if (ierr.lt.0) then
          print*, 'erreur gdll', ierr, gr % std_grd
          call my_abort
        endif

      case('X','O','M')

         nx = gr % nx
         ny = gr % ny

       !! Lecture, Lat-Lon nemo

         ! longitude
         wrkr4 => gr % lld % lont
         ierr = fstinl(unitf, nx, ny, nk, -1, ' ', -1, -1, -1, ' ', '>>',liste, nlist, 10)
         if ( ierr /= 0 ) then
            write(*,*) 'found no record for >>'
            call my_abort
         endif
         l = nlist
         ierr = fstluk(wrkr4,liste(l),nx,ny,nk)
         if ( ierr < 0 ) then
            write(*,*) 'problem reading >>'
            call my_abort
         endif
         ! latitude
         wrkr4 => gr % lld % latt
         ierr = fstinl(unitf, nx, ny, nk, -1, ' ', -1, -1, -1, ' ', '^^',liste, nlist, 10)
         if ( ierr /= 0 ) then
            write(*,*) 'found no record for ^^'
            call my_abort
         endif
         l = nlist
         ierr = fstluk(wrkr4,liste(l),nx,ny,nk)
         if ( ierr < 0 ) then
            write(*,*) 'problem reading ^^'
            call my_abort
         endif

         call construct_tree( gr )
      case default
      end select

      write(*,*) 'latlon read successfully'
      write(*,*) 'min/max lon',minval(gr % lld % lont),maxval(gr % lld % lont)
      write(*,*) 'min/max lat',minval(gr % lld % latt),maxval(gr % lld % latt)

      ! special treatment for unstructured grid
      if ( grd % grtyp == 'M' ) call read_mesh_connectivity( gr, unitf )
 
    enddo
  END SUBROUTINE getll

  subroutine read_mesh_connectivity( grid, unitf )
    !! --------------------------------------------------------------------------------------
    !!  *** allocate and read the connectivity of the mesh
    !! --------------------------------------------------------------------------------------
    ! arguments
    integer, intent(in) :: unitf ! RPN file unit
    type(type_grid), intent(inout), target :: grid
    ! locals
    type(type_mesh), pointer :: mesh_p
    integer, dimension(:,:), pointer :: in
    integer nn, ne, nconnect
    integer ni, nj, nk
    integer, parameter :: nmax = 10
    integer nlist, liste(nmax), l
    integer ierr
    integer i,j,n
    real(kind=4), dimension(:), pointer :: lon,lat
   ! librmn functions
    integer, external   :: fstinl, fstluk

    mesh_p => grid % meshd

    ierr = fstinl(unitf, ni, nj, nk, -1, ' ', -1, -1, -1, ' ', '##',liste, nlist, nmax)
    if ( ierr /= 0 ) then
       write(*,*) 'found no record for ##'
       call my_abort
    endif

    ! allocate the connectivity
    !
    nconnect = 3 
    ne = ni * nj / nconnect
    if ( ne * nconnect /= ni * nj ) call my_abort

    allocate( mesh_p % in(nconnect,ne) )
    in => mesh_p % in
    mesh_p % ne = ne
    mesh_p % nconnect = nconnect

    ! read the connectivity
    !
    l = nlist
    ierr = fstluk(in,liste(l),ni,nj,nk)
    if ( ierr < 0 ) then
     write(*,*) 'problem reading ##'
       call my_abort
    endif

    nn = grid % nx * grid % ny
    mesh_p % nn = nn
    allocate( mesh_p % longr(nn) )
    allocate( mesh_p % latgr(nn) )
    lon => mesh_p % longr
    lat => mesh_p % latgr
    n = 0
    do j = 1, grid % ny
      do i = 1, grid % nx
         n = n + 1
         lon(n) = grid % lld % lont(i,j)
         lat(n) = grid % lld % latt(i,j)
      enddo
    enddo

    ! make sure that the connecticvity is between 1 and nn
    in = in + 1

    write(*,*) 'connectivity read successfully'
    write(*,*) 'min/max in',minval(in),maxval(in)

  end subroutine read_mesh_connectivity

  SUBROUTINE allocate_grid_intrp( grid )
    !! -----------------------------------------------------------------------------------------------
    !!  *** allocate or reallocate interpolation pointers
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    ! locals
    integer nx,ny

    nx = grid % std_grd % ni
    ny = grid % std_grd % nj

    if ( .not.associated( grid % lld % mask) ) &
    allocate( grid % lld % mask ( nx, ny ) )

    if ( associated( grid % locd % Ag) ) return
    allocate( grid % locd % Ag ( nx, ny ) )
    allocate( grid % locd % Bg ( nx, ny ) )
    allocate( grid % locd % ig ( nx, ny ) )
    allocate( grid % locd % jg ( nx, ny ) )

  END SUBROUTINE allocate_grid_intrp

  SUBROUTINE allocate_grid_weights( grid )
    !! -----------------------------------------------------------------------------------------------
    !!  *** allocate or reallocate interpolation pointers
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    ! locals
    integer nx,ny, nw

    if ( associated( grid % wgtd % wgt) ) return
    nw = grid % wgtd % numwgt
    nx = grid % std_grd % ni
    ny = grid % std_grd % nj
    allocate( grid % wgtd % nwgt ( nx, ny ) )
    allocate( grid % wgtd %  wgt ( nx, ny, nw ) )
    allocate( grid % wgtd % iwgt ( nx, ny, nw ) )
    allocate( grid % wgtd % jwgt ( nx, ny, nw ) )

  END SUBROUTINE allocate_grid_weights

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  SUBROUTINE construct_tree(grid)
    !! -----------------------------------------------------------------------------------------------
    !!  *** construct a reduced tree, especially useful when the grid is quasi-regular
    !! -----------------------------------------------------------------------------------------------
    IMPLICIT NONE
    ! * Arguments
    type (type_grid), target &
                        :: grid  ! grid type 
    ! locals
    integer nx,ny
    REAL(KIND=4) :: start, finish
    INTEGER ntot, ji, jj, i1, i2
    REAL(kind=4), ALLOCATABLE, DIMENSION(:,:) :: data3d
    INTEGER dxtree, dytree
    INTEGER, POINTER ::  nxred, nyred
    TYPE (kdtree2), pointer :: my_kdtree
    INTEGER, DIMENSION(:), POINTER :: iil, jjl
    REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: lonred, latred
    REAL(KIND=4), DIMENSION(:,:), POINTER :: lon, lat

    call cpu_time(start)

! create the KD-tree for the source grid
    nx     =  grid % std_grd % ni
    ny     =  grid % std_grd % nj
    lon    => grid % lld % lont
    lat    => grid % lld % latt
    nxred  => grid % nxred
    nyred  => grid % nyred
    dxtree =  grid % dxtree
    dytree =  grid % dytree

    allocate( grid % kdtree )
    my_kdtree => grid % kdtree

    ! create a coarsened grid
    nxred = (nx-1) / dxtree + 1
    nyred = (ny-1) / dytree + 1
    dxtree = grid % dxtree
    dytree = grid % dytree

    ALLOCATE( lonred(nxred,nyred), &
              latred(nxred,nyred))
    ! allocate reduced index array
    ALLOCATE( grid % ired(nxred), grid % jred(nyred) )
    iil    => grid % ired
    jjl    => grid % jred

    lonred(:,:) = lon(1:nx:dxtree,1:ny:dytree)
    latred(:,:) = lat(1:nx:dxtree,1:ny:dytree)
    i2 = 0
    DO jj=1, ny, dytree
       i2 = i2 + 1
       jjl(i2) = jj
    ENDDO
    i1 = 0
    DO ji=1, nx, dxtree
       i1 = i1 + 1
       iil(i1) = ji
    ENDDO

    ntot = nxred * nyred
    allocate(data3d(3,ntot))
    call ez_lac(data3d, lonred, latred, ntot)
    deallocate( lonred, latred )
    !
    my_kdtree = kdtree2_create(data3d)
    deallocate(data3d)
    my_kdtree % the_data => NULL()
    write(*,*) 'kdtree',grid % kdtree % n, grid % kdtree % dimen
    call cpu_time(finish)
    print '("KD-tree creation Time = ",f6.3," seconds.")',finish-start

  END SUBROUTINE construct_tree

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  subroutine init_grid_mask(grid,nomvar,unitf,ln_mask)
    implicit none
    ! arguments
    type(type_grid), target :: grid
    character(len=*), intent(in) :: nomvar
    integer, intent(in) :: unitf
    logical, intent(in) :: ln_mask
    ! locals
    type(type_grid), dimension(:), pointer :: subg
    integer ierr, status_st90
    integer nx, ny
    integer nlist, liste(10), l, nk
    integer, allocatable, dimension(:,:) :: wrki4
    integer ig, ng, i1
    type(type_grid), pointer :: gr ! pointer grid
   ! librmn functions
    integer, external   :: fstinl, fstluk
    logical, dimension(:,:), pointer :: mask

    nx = grid % nx
    ny = grid % ny
    nk = 1
    subg => grid % subgrd
    ng   =  grid % ng

    if (associated(grid % lld % mask)) deallocate(grid % lld % mask)
    allocate( grid % lld % mask(nx,ny) )
    mask => grid % lld % mask
    mask(:,:) = .true. ! all points are good!

    do ig = 1, ng
       gr => subg(ig)
       nx = gr % nx
       ny = gr % ny
       if (associated(gr % lld % mask)) deallocate(gr % lld % mask)
       allocate( gr % lld % mask(nx,ny) )
       gr % lld % mask(:,:) = .true.
    enddo

    if (.not.ln_mask) then
       write(*,*) 'setting grid mask to T', grid % std_grd
       return ! if no request for mask, return and continue without
    endif

    ierr = fstinl(unitf, nx, ny, nk, -1, ' ', -1, -1, -1, '@@', nomvar,liste, nlist, 10)
    if ( ierr /= 0 .or. nlist == 0) then
       write(*,*) 'found no @@ record for mask of ',TRIM(nomvar)
       write(*,*) 'of grid ',grid % std_grd
       write(*,*) 'setting grid mask to T'
       return
    endif

    l = nlist
    allocate(wrki4(nx,ny))
    ierr = fstluk(wrki4,liste(l),nx,ny,nk)
    where(wrki4 == 0 ) mask = .false.
    deallocate(wrki4)

    ! distribute to sub grid
    select case(ng)
    case(1)
       gr => subg(ng)
       gr % lld % mask = grid % lld % mask
    case(2) ! Yin-Yang grid, cut the mask in 2
       do ig = 1, ng
          gr => subg(ig)
          ny = gr % ny
          i1 = ny * ( ig - 1 )
          gr % lld % mask(:,:) = grid % lld % mask(:,i1+1:i1+ny)
       enddo
    end select

  end subroutine init_grid_mask

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  subroutine remove_periodic_points(grd, M)
  ! remove duplicate points
  implicit none
  ! arguments
  type(type_grid), intent(in) :: grd
  logical, dimension(:,:), intent(inout) :: M
  ! locals
  integer nx,ny,i,j

    nx = SIZE( M, 1 )
    ny = SIZE( M, 2 )

    ! start with west-east
    if( grd % webdy ) then ! w-e periodicity active
      select case( grd % werepeat )
      case(1)
        M(nx,:) = .false.
      case(2)
        M( 1,:) = .false.
        M(nx,:) = .false.
      end select
    endif

    select case( grd % nbdy )
    case('F') ! F-folding
      do j = 1, grd % nrepeat
        M(:,ny-j+1) = .false.
      enddo
    case('T') ! T-folding
      do j = 1, grd % nrepeat
        M(:,ny-j+1) = .false.
      enddo
      j = ny - grd % nrepeat
        M(nx/2+1:nx,j) = .false.
    end select

  end subroutine remove_periodic_points

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  subroutine create_corner_mesh( grid )
    implicit none
    ! arguments
    type(type_grid), target :: grid
    ! locals
    real(kind=4), dimension(:,:,:), allocatable :: data3d, data3m
    integer ntot, ntot2, nx, ny, nx2, ny2
    integer i,j

    nx = grid % nx; ny = grid % ny
    ntot = nx * ny
    allocate(data3d(3,nx,ny))
    call ez_lac(data3d, grid % lld % lont, grid % lld % latt, ntot)

    nx2 = nx + 1; ny2 = ny + 1
    allocate(data3m(3,nx,0:ny))

    ! south border extrapolation
    j = 1
    do i = 1, nx
          data3m(:, i , j-1) = 1.5 *   data3d(:, i , j  ) &
                             - 0.5 *   data3d(:, i , j+1)
    enddo

    ! central averaging
    do j = 2, ny
       do i = 1, nx
          data3m(:, i , j-1) = 0.5 * ( data3d(:, i , j  ) &
                                     + data3d(:, i , j-1) )
       enddo
    enddo

    ! north border extrapolation
    j = ny
    do i = 1, nx
          data3m(:, i , j  ) = 1.5 *   data3d(:, i , j  ) &
                             - 0.5 *   data3d(:, i , j-1)
    enddo
    deallocate(data3d)

    ntot = nx2 * ny2
    allocate(data3d(3,0:nx,0:ny))

    ! west border extrapolation
    i = 1
    do j = 0, ny
          data3d(:, i-1, j) = 1.5 *   data3m(:, i   , j) &
                            - 0.5 *   data3m(:, i+1 , j)
    enddo

    ! central averaging
    do j = 0, ny
       do i = 2, nx
          data3d(:, i-1, j) = 0.5 * ( data3m(:, i   , j) &
                                    + data3m(:, i-1 , j) )
       enddo
    enddo

    ! east border extrapolation
    i = nx
    do j = 0, ny
          data3d(:, i  , j) = 1.5 *   data3m(:, i   , j) &
                            - 0.5 *   data3m(:, i-1 , j)
    enddo
    deallocate(data3m)

    ! allocate corner mesh points
    allocate( grid % lld % lonc(0:nx,0:ny), &
              grid % lld % latc(0:nx,0:ny) )
    ! convert back
    call ez_cal( grid % lld % lonc, grid % lld % latc, data3d, ntot)
    deallocate(data3d)

  end subroutine create_corner_mesh

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  subroutine get_delta_i_xyz( grid, di3d )
    implicit none
    ! arguments
    type(type_grid), target :: grid
    real(kind=4), dimension(:,:,:) :: di3d
    ! locals
    real(kind=4), dimension(:,:,:), allocatable :: di3dm
    integer ntot, ntot2, nx, ny, nx2, ny2
    integer i,j,k

    nx = grid % nx; ny = grid % ny
    if (.not.associated(grid % lld % t3d )) then
       allocate( grid % lld % t3d(3,nx,ny) )
       call ez_lac(grid % lld % t3d, grid % lld % lont, grid % lld % latt, nx*ny)
    endif

    ! coompute delta of the cartesian coordinate along the i-direction, centered on U-point
    allocate(di3dm(3,0:nx,ny))
    do j = 1, ny
       do i = 1, nx - 1
          do k = 1, 3
             di3dm(k, i, j) = 0.5 * ( grid % lld % t3d(k, i+1 , j) &
                                    - grid % lld % t3d(k, i   , j) )
          enddo
       enddo
       ! extrapolation for last row
       i = 0
          do k = 1, 3
             di3dm(k, i , j) = 2.0 * ( di3dm(k, i+1 , j) &
                                     - di3dm(k, i+2 , j) )
          enddo
       i = nx
          do k = 1, 3
             di3dm(k, i , j) = 2.0 * ( di3dm(k, i-1 , j) &
                                     - di3dm(k, i-2 , j) )
          enddo
    enddo
    ! central averaging from U-point to T-pointn following NEMO convention
    do j = 1, ny
       do i = 1, nx
          do k = 1, 3
             di3d(k, i , j) = 0.5 * ( di3dm(k, i   , j) &
                                    + di3dm(k, i-1 , j) )
          enddo
       enddo
    enddo
    deallocate(di3dm)

    end subroutine get_delta_i_xyz

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
  subroutine get_dlatlon_xyz( grid, di3d )
    implicit none
    ! arguments
    type(type_grid), target :: grid
    real(kind=4), dimension(:,:,:), intent(out) :: di3d
    ! locals
    real(kind=4), dimension(:,:,:), allocatable :: di3dm
    integer ntot, ntot2, nx, ny, nx2, ny2
    integer i,j,k

    nx = grid % nx; ny = grid % ny
    if (.not.associated(grid % lld % t3d )) then
       allocate( grid % lld % t3d(3,nx,ny) )
       call ez_lac(grid % lld % t3d, grid % lld % lont, grid % lld % latt, nx*ny)
    endif

    ! force vector to follow tangent to longitude circle
    do j = 1, ny
       do i = 1, nx
          di3d(1, i, j) = - grid % lld % t3d(2, i, j)
          di3d(2, i, j) =   grid % lld % t3d(1, i, j)
          di3d(3, i, j) =  zerd
       enddo
    enddo

    end subroutine get_dlatlon_xyz

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------

  !! --------------------------------------------------------------------------------------
  !! --------------------------------------------------------------------------------------
END MODULE grids
