MODULE llxy

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE llxy ***
  !!
  !!  **  Purpose: General grid conversion program lat-lon to x,y
  !!      the inverse, F. Roy
  !!
  !!  **  Method: Manage custom structured grids with masks and angles
  !!  **          assuming angle rotation at center of the grid
  !!
  !! history :
  !!     Original code  :   F. Roy (November 2009)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

  !! Modules
  USE nemo_proj
  USE triangle_proj

  IMPLICIT NONE

  PUBLIC ll2xy

  CONTAINS

    SUBROUTINE ll2xy(grds,grdr,loc)

    !! Compute X,Y components and angle A from grid destination
    !! relative to source grid
    !! A is relative to X direction (i==>nx) 
    !! and positive anticlockwise (radians)
    !! Return values of -9. for masked points, M=0

    USE intrp_oce

    IMPLICIT NONE
    ! arguments
    TYPE(type_grid), TARGET, INTENT(inout) :: grds, grdr
    TYPE(type_location), TARGET, INTENT(inout) :: loc
    ! locals
    REAL(KIND=4), DIMENSION(:,:), POINTER :: lat,lon,A
    REAL(KIND=8), DIMENSION(:,:), POINTER :: AA,BB
    INTEGER,      DIMENSION(:,:), POINTER :: II,JJ
    LOGICAL,      DIMENSION(:,:), POINTER :: M
    INTEGER nx, ny, nxs, nys
    REAL(KIND=4), DIMENSION (:,:), ALLOCATABLE :: x,y
    integer ier, gid
    REAL(KIND=8) :: w1,w2,w3,w4
    REAL(KIND=8) :: as1,as2,as3,as4,asm
    REAL(KIND=8) :: csm,snm,cs1,cs2,cs3,cs4,sn1,sn2,sn3,sn4
    integer i,j,is,js, i1, i2, i3
    integer, external :: gdxyfll ! rmn functions
    REAL(KIND=4), allocatable, dimension(:,:,:) :: di3ds, di3dr
    REAL(Kind=4), dimension(3) :: di1, di2, di3, di4, dir, dis
    

    !-----------------------------------------------------------------------------------------
    ! find the fractional position on the source grid of the lat/lon of the destination grid
    !-----------------------------------------------------------------------------------------

    lon => grdr % lld % lont
    lat => grdr % lld % latt
    nx  =  grdr % nx
    ny  =  grdr % ny
    gid  = grds % std_grd % grid_id
    allocate( loc % ag(nx,ny) &
            , loc % bg(nx,ny) &
            , loc % ig(nx,ny) &
            , loc % jg(nx,ny) &
            , loc % pg(nx,ny) &
            , loc % mg(nx,ny) )
    if (grds % std_grd % grtyp .eq. 'Y') then
      write(*,*) 'CSTINTRP: This call of ll2xy should not have happened, gdefs="cloud"'
      call my_abort
    endif
    AA => loc % ag
    BB => loc % bg
    II => loc % ig
    JJ => loc % jg
    A  => loc % pg; A(:,:) = 0.
    M  => loc % mg; M(:,:) = .true.

    select case( grds % std_grd % grtyp )
    case('Z','B','U','A','N','S','L')
      allocate( x(nx,ny), y(nx,ny) )
      ier = gdxyfll(gid,x,y,lat,lon,nx*ny)
      II(:,:) = int(x(:,:))
      JJ(:,:) = int(y(:,:))
      AA(:,:) = x(:,:) - II(:,:)
      BB(:,:) = y(:,:) - JJ(:,:)
      deallocate(x,y)
    case('M')
      call ll2xy_triangle(grds,lat,lon,AA,BB,II,JJ,M,nx,ny)
    case default
      call ll2xy_nemo(grds,lat,lon,AA,BB,II,JJ,M,nx,ny)
    end select

    ! set the mask to zero if along or outside the border of source grid
    ! is this required for all grids/interpolator (i.e. nearst) ??
    nxs = grds % std_grd % ni
    nys = grds % std_grd % nj
    if ( grds % std_grd % grtyp /= 'M' ) &
       where(ii < 1 .or. ii > nxs-1 .or. jj < 1 .or. jj > nys-1) M = .false.

    !-----------------------------------------------------------------------
    ! Compute the difference angle between the destination and source grids
    !-----------------------------------------------------------------------

    allocate( di3ds(3, nxs, nys ))
    allocate( di3dr(3, nx , ny  ))

    select case( grds % std_grd % grtyp )
    case('M','Y')
      call get_dlatlon_xyz(  grds, di3ds )
    case default
      call get_delta_i_xyz( grds, di3ds )
    end select

    select case( grdr % std_grd % grtyp )
    case('M','Y')
      call get_dlatlon_xyz(  grdr, di3dr )
    case default
      call get_delta_i_xyz( grdr, di3dr )
    end select

    do j=1,ny
    do i=1,nx
      if (M(i,j)) then
        dir = di3dr(:,i,j)
        is = ii(i,j)
        js = jj(i,j)
        if ( grds % std_grd % grtyp == 'M') then
          w1 = AA(i,j)
          w2 = BB(i,j)
          w3 = MIN( MAX(zerd, oned - w1 - w2), oned)
          i1 = grds % meshd % in(1,is)
          i2 = grds % meshd % in(2,is)
          i3 = grds % meshd % in(3,is)
          di1 = di3ds(:,i1,js)
          di2 = di3ds(:,i2,js)
          di3 = di3ds(:,i3,js)
          dis = w1 * di1 + w2 * di2 + w3 * di3 ! interpolated vector onto destination point
        else
          di1 = di3ds( :, is  , js  )
          di2 = di3ds( :, is+1, js  )
          di3 = di3ds( :, is+1, js+1)
          di4 = di3ds( :, is  , js+1)
          w1 = (1._8-AA(i,j)) * ( 1._8-BB(i,j))
          w2 =       AA(i,j)  * ( 1._8-BB(i,j))
          w3 =       AA(i,j)  *        BB(i,j)
          w4 = (1._8-AA(i,j)) *        BB(i,j)
          dis = w1 * di1 + w2 * di2 + w3 * di3 + w4 * di4 ! interpolated vector onto destination point
        endif
        csm = dot_product( dis, dir )
        snm = vector_pro( dis, dir, grdr % lld % t3d(:,i,j) ) ! the cross product is reduced to the normal component to the sphere
        asm = atan2(snm, csm)
        A(i,j) = asm
      endif
    enddo
    enddo
    deallocate(di3ds, di3dr)

    END SUBROUTINE ll2xy

    real(kind=4) function vector_pro( d1, d2, xyz )
    implicit none
    ! arguments
    real(kind=4), dimension(3) :: d1, d2, xyz
    ! local
    real(kind=4), dimension(3) :: d3
    real(kind=4) :: nd

    ! compute cross product between d1 and d2
    d3(1) = d1(2) * d2(3) - d1(3) * d2(2)
    d3(2) = d1(3) * d2(1) - d1(1) * d2(3)
    d3(3) = d1(1) * d2(2) - d1(2) * d2(1)

    ! compute the vector product as the component along the normal to the sphere, which xyz here (unit vector)
    vector_pro = dot_product( d3, xyz )

    end function vector_pro

END MODULE llxy
