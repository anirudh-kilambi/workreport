MODULE triangle_proj

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE nemo_proj ***
  !!
  !!  **  Purpose: 1) Read nemo coordinates and reperes
  !!               (getll_nemo_ext,getreps_nemo)
  !!               2) Convert ll to xy and xy to ll
  !!               (ll2xy_nemo,...)
  !!
  !!  **  Method: 1) Read lat-lons
  !!              2) Manage localisation into nemo grid
  !!              using local PS projections
  !!
  !! history :
  !!     Original code : F. Roy (October 2012)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

  USE grids
  USE mesh
  USE kdtree2_module
  USE mympp

  IMPLICIT NONE

  PUBLIC ll2xy_triangle

  Type (kdtree2), pointer :: my_kdtree
  Type (type_mesh), pointer :: meshp

  CONTAINS

    !! --------------------------------------------------------------------------------------
    SUBROUTINE ll2xy_triangle(grid,lat,lon,AA,BB,ii,jj,M,nx,ny)

    !! Compute X,Y components from nemo grid
    !! Return values of -9. for masked points, M=0

    !! I/O Variables
    IMPLICIT NONE
    TYPE(type_grid), TARGET :: grid
    INTEGER , INTENT(in) :: nx,ny
    REAL(KIND=4), INTENT(in), DIMENSION (nx,ny) :: lat,lon
    REAL(KIND=8), INTENT(inout), DIMENSION (nx,ny) :: AA,BB
    INTEGER,      INTENT(inout), DIMENSION (nx,ny) :: ii,jj
    LOGICAL,      INTENT(inout), DIMENSION (nx,ny) :: M
    !! Local variables
    LOGICAL :: maskr
    INTEGER :: ie,je    !Grille entree
    REAL(KIND=4), PARAMETER :: r=1.0                    !Projection scale
    INTEGER n_series, ipt, icell
    REAL(KIND=4) :: start, finish, my_xyz(3)
    INTEGER,      allocatable, dimension(:) :: list_src
    TYPE(kdtree2_result), allocatable, dimension(:), target :: my_res

    my_kdtree => grid % kdtree
    meshp => grid % meshd

   !--------
   ! construct source mesh
   !--------
    nn => meshp % nn
    ne => meshp % ne
    nconnect  => meshp % nconnect
    longr => meshp % longr
    latgr => meshp % latgr
    in => meshp % in
    call pointer_mesh( meshp, .true. )

   !--------
    ! define and allocate search radius
   !--------
    n_series = 1 ! just give the closest!
    allocate (my_res(n_series))
    allocate (list_src(n_series))

    write(*,*) 'my_kdtree',my_kdtree%n, my_kdtree%dimen
    call cpu_time(start)

    !! Boucle sur la grille d'entree
    DO je=1,ny
    DO ie=1,nx
      
      !! ... find closest index from regular grid
      call ez_lac(my_xyz, lon(ie,je), lat(ie,je), 1)
      n_series = 1
      call kdtree2_n_nearest(my_kdtree,my_xyz,n_series,my_res)
      do ipt = 1, n_series
         list_src(ipt) = my_res(ipt) % idx
      enddo
      call refine_search_tri(list_src, icell, n_series, lon(ie,je), lat(ie,je))
      call calculate_tri_on_src(AA(ie, je), BB(ie, je), &
                         ii(ie, je), jj(ie, je), M(ie, je), lon(ie,je), lat(ie,je), &
                         icell, r)
    ENDDO
    ENDDO

    deallocate(my_res, list_src)
!    call dealloc_llm_nemo(maskr)
    call cpu_time(finish)
    print '("ll2xy Time = ",f6.3," seconds.")',finish-start

    END SUBROUTINE ll2xy_triangle
    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
    SUBROUTINE refine_search_tri(list_node, icell, n_series, lonw, latw)
      !! Modules
      USE intrp_oce

    IMPLICIT NONE
    ! arguments
    INTEGER, INTENT(in) :: n_series
    INTEGER, DIMENSION(n_series), INTENT(IN) :: list_node
    REAL(KIND=4), INTENT(in) :: lonw, latw
    INTEGER, INTENT(out) :: icell
    ! locals
    REAL(KIND=8) :: srchpt_x, srchpt_y, latwd, lonwd
    INTEGER ipt, j, i_loc, t_loc
    real(kind=8), dimension(nconnect) :: lat_tri, lon_tri
    real(kind=8), dimension(nconnect) :: xpol_src, ypol_src
    logical inpoly

    ! the point on the destination grid is at the center of the projection
    srchpt_x = zerd; srchpt_y = zerd
    lonwd = lonw; latwd = latw

    do ipt = 1, n_series
       i_loc = list_node(ipt)
       do t_loc = meshp % nmit(i_loc) , meshp % nmit(i_loc+1) - 1
          icell = meshp % tmit(t_loc)
          lon_tri(:) = longr( in (:,icell) )
          lat_tri(:) = latgr( in (:,icell) )
          do j = 1, nconnect
             call ll_to_polars(oned, latwd, lonwd, lat_tri(j), lon_tri(j), xpol_src(j), ypol_src(j))
          enddo
          inpoly = inpoly_complex(xpol_src, ypol_src, nconnect, srchpt_x, srchpt_y)
          if( inpoly ) return
       enddo
    enddo

    ! if it fails finding a close-by element, set the cell id to zero
    icell = 0

    END SUBROUTINE refine_search_tri
    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
    SUBROUTINE calculate_tri_on_src(AA, BB, ii, jj, msk, lonw, latw, icell, r)
      ! Object: calclualte X,Y and msk of the target grid with respect to the source grid

      !! Modules
      USE intrp_oce


      implicit none

      ! output coordinates and mask of the target point wrt source grid
      real (kind=8), intent(inout) :: AA, BB
      integer,       intent(inout) :: ii,jj
      logical,       intent(inout) :: msk

      real (kind=4), intent(in)  :: lonw, latw ! coordinates of the target grid point
      !ranges of indices of the source gridpoints to be considered for the given target
      integer, intent (in) :: icell
      REAL(KIND=4), intent (in) :: r ! projection scale
      ! locals
      REAL(KIND=8) :: srchpt_x, srchpt_y, area, vecp, dxa, dya, dxb, dyb, vect, lonwd, latwd, rd
      real(kind=8), dimension(nconnect) :: lat_tri, lon_tri
      real(kind=8), dimension(nconnect) :: xpol_src, ypol_src
      integer j, i1, i2, i3

      rd = real(r,kind=8) ! conversion to double precision

      AA = zerd; BB = zerd; ii = -9; jj = -9; msk = .false.
      if( icell == 0 ) return

      msk = .true. ! if icell > 0 then we can find a value

      ! the point on the destination grid is at the center of the projection
      srchpt_x = zerd; srchpt_y = zerd
      lonwd = lonw; latwd = latw

      lon_tri(:) = longr( in (:,icell) )
      lat_tri(:) = latgr( in (:,icell) )
      do j = 1, nconnect
         call ll_to_polars(oned, latwd, lonwd, lat_tri(j), lon_tri(j), xpol_src(j), ypol_src(j))
      enddo

      do i1 = 1, nconnect - 1
        i2=mod(i1,nconnect)+1
        i3=mod(i2,nconnect)+1
        dxa = xpol_src(i2) - xpol_src(i1)
        dya = ypol_src(i2) - ypol_src(i1)
        dxb = srchpt_x     - xpol_src(i1)
        dyb = srchpt_y     - ypol_src(i1)
        vecp = dyb*dxa - dya*dxb
        dxb = xpol_src(i3) - xpol_src(i1)
        dyb = ypol_src(i3) - ypol_src(i1)
        vect = dyb*dxa - dya*dxb
        if (vect < epsilon ) then
           write(*,*) 'bad triangle, bailing out'
        endif
        if (i1==1) then 
           AA = MIN( MAX( zerd, vecp / vect), oned)
           ii = icell
        endif
        if (i1==2) then 
           BB = MIN( MAX( zerd, vecp / vect), oned)
           jj = 1
        endif
      enddo

    END SUBROUTINE calculate_tri_on_src
    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
END MODULE triangle_proj
 
