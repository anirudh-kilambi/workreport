MODULE nemo_proj

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
  USE kdtree2_module
  USE mympp

  IMPLICIT NONE

  PUBLIC ll2xy_nemo

  Type (kdtree2), pointer :: my_kdtree
  INTEGER werepeat, nrepeat, dxtree, dytree
  LOGICAL webdy
  CHARACTER(LEN=1) nbdy
  REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: lonred, latred
  INTEGER nxred, nyred

  INTEGER, DIMENSION(:), POINTER :: iil, jjl
  REAL(KIND=4), DIMENSION(:,:), POINTER :: longrd, latgrd
  INTEGER nxgrd, nygrd
  integer ipl, jpl, ipr, jpr, ipt, jpt

  CONTAINS

    !! --------------------------------------------------------------------------------------
    SUBROUTINE ll2xy_nemo(grid,lat,lon,AA,BB,ii,jj,M,nx,ny)

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
    REAL(KIND=8) Xw,Yw
    INTEGER ntot, n_series, ipt, ii0, jj0, ilast, jlast, idx, ilast0, jlast0
    REAL(KIND=4) :: start, finish, my_xyz(3)
    INTEGER,      allocatable, dimension(:) :: list_isrc, list_jsrc
    TYPE(kdtree2_result), allocatable, dimension(:), target :: my_res
    INTEGER dxtree, dytree, nxred, nyred

    my_kdtree => grid % kdtree
    dxtree    =  grid % dxtree
    dytree    =  grid % dytree
    nxred     =  grid % nxred
    nyred     =  grid % nyred
    longrd    => grid % lld % lont
    latgrd    => grid % lld % latt
    nxgrd     =  grid % nx
    nygrd     =  grid % ny
    iil       => grid % ired
    jjl       => grid % jred
    nbdy      =  grid % nbdy
    nrepeat   =  grid % nrepeat
    webdy     =  grid % webdy
    werepeat  =  grid % werepeat
    select case(werepeat)
    case(0) ! CICE case
      ipl = 1
      ipr = nxgrd
    case(1) ! GEM global grid
      ipl = 1
      ipr = nxgrd - 1
    case(2) ! ORCA grid
      ipl = 2
      ipr = nxgrd - 1
    end select
    jpt = nygrd - nrepeat

    ! define and allocate search radius
    n_series = 1 ! just give the closest!
    allocate (my_res(n_series))
    n_series = 4
    allocate (list_isrc(n_series), &
              list_jsrc(n_series))

    write(*,*) 'my_kdtree',my_kdtree%n, my_kdtree%dimen
    call cpu_time(start)

    !! Boucle sur la grille d'entree
    DO je=1,ny
    DO ie=1,nx
      ii(ie,je) =-9
      jj(ie,je) =-9
      M (ie,je) = .false.
      AA(ie,je) = 0.
      BB(ie,je) = 0.
      !! ... find closest index from regular grid
      if (mod(ie-1,dxtree)==0) then
         call ez_lac(my_xyz, lon(ie,je), lat(ie,je), 1)
         n_series = 1
         call kdtree2_n_nearest(my_kdtree,my_xyz,n_series,my_res)
         ipt = 1
         idx = my_res(ipt) % idx
         jj0 = (idx-1) / nxred + 1
         ii0 = idx - nxred * (jj0 - 1)
         ilast = iil(ii0)
         jlast = jjl(jj0)
         ilast0 = ilast; jlast0 = jlast
      endif
      n_series = 4
      call refine_search(ilast, jlast, n_series, list_isrc, list_jsrc, lon(ie,je), lat(ie,je))
      if( n_series == 0 ) then
         ! reset the search indices to the last successful ones
         ilast = ilast0; jlast = jlast0
         cycle
      endif
      call calculate_XYM_dst_on_src(AA(ie, je), BB(ie, je), &
                         ii(ie, je), jj(ie, je), M(ie, je), lon(ie,je), lat(ie,je), &
                         n_series, list_isrc, list_jsrc, r)
      if (M(ie,je)) then
         ilast = ii(ie,je)
         jlast = jj(ie,je)
         ilast0 = ilast; jlast0 = jlast
      endif
    ENDDO
    ENDDO

    deallocate(my_res, list_isrc, list_jsrc)
!    call dealloc_llm_nemo(maskr)
    call cpu_time(finish)
    print '("ll2xy Time = ",f6.3," seconds.")',finish-start

    END SUBROUTINE ll2xy_nemo
    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
    SUBROUTINE calculate_XYM_dst_on_src(AA, BB, ii, jj, msk, lonw, latw, npt, ipt, jpt, r)
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
      integer, intent (in) :: npt, ipt(npt), jpt(npt)
      REAL(KIND=4), intent (in) :: r ! projection scale
      ! locals
      integer :: iii
      real (kind=8) :: Xw, Yw, w1, w2, w3, w4, rd
      REAL(KIND=8), DIMENSION (4) :: lonin, latin, xpo, ypo  !nemo projection on PS
      logical :: inp ! =true if inpoly
      integer iip

      rd = real(r,kind=8) ! conversion to double precision

      DO iip = 1, npt
        ii = ipt(iip)
        jj = jpt(iip)
        ! take the points anticlockwise and convert to double precision
        lonin(1)=longrd(ii,jj)
        latin(1)=latgrd(ii,jj)
        lonin(2)=longrd(ii+1,jj)
        latin(2)=latgrd(ii+1,jj)
        lonin(3)=longrd(ii+1,jj+1)
        latin(3)=latgrd(ii+1,jj+1)
        lonin(4)=longrd(ii,jj+1)
        latin(4)=latgrd(ii,jj+1)

        !! ... Project four corners on PS grid
        Xw = lonw; Yw = latw ! convert to double precision
        do iii=1,4
        ! all double precision from here on
            call ll_to_polars(rd,Yw,Xw,    &
   &                          latin(iii),lonin(iii), &
   &                          xpo(iii),ypo(iii))
        enddo
        Xw = zerd; Yw = zerd ! the projection is centered on the interpolated point
        !! ... Project entry grid lat-lon
        !! ... Check if point is inside nemo polygon
        inp=inpoly_complex(xpo,ypo,4,Xw,Yw)

        if (inp) then
          !! ... Project four corners on PS grid
          !! Interpolate x,y coordinates according to grid index ii,jj
          call interp_gen_2D(xpo(1),ypo(1),xpo(2),ypo(2),    &
   &                         xpo(3),ypo(3),xpo(4),ypo(4),    &
   &                         Xw,Yw,AA,BB)
          msk = .true.
          return ! bail out, no need to continue, the first accceptable polygon is always the best!
        endif ! end condition on inside polygon
      ENDDO

    END SUBROUTINE calculate_XYM_dst_on_src

    !! --------------------------------------------------------------------------------------

    !! --------------------------------------------------------------------------------------
    SUBROUTINE refine_search(ilast, jlast, n_series, list_isrc, list_jsrc, lonw, latw)
      !! Modules
      USE intrp_oce

    IMPLICIT NONE
    ! arguments
    INTEGER ilast, jlast, n_series, list_isrc(n_series), list_jsrc(n_series)
    REAL(KIND=4) :: lonw, latw
    ! locals
    REAL(KIND=8) :: lonwd, latwd, x1, y1, x2, y2, xg, yg
    REAL(KIND=8) :: lonin, latin, m012, x0g2, x01g, a, b
    INTEGER inext, jnext, di, dj, d0, d, counter, ilast0, jlast0

    a = 1; b = 1
    inext = -nxgrd ; jnext = -nygrd
    ilast0 = inext; jlast0 = jnext
    d = abs(ilast - inext) + abs(jlast - jnext)
    d0 = d + 1
    counter = 0
!    ! apply periodic bc on ref point
    call apply_periodic_condition_search(ilast, jlast, 0, 0, inext, jnext)
    ilast = inext; jlast = jnext

    do while ( d < d0 ) ! the next iteration should improve on the previous. if not bail out

      ! all double precision from here on
      ! projection done relative to point ilast, jlast => (0,0)
      lonwd = longrd(ilast,jlast); latwd = latgrd(ilast,jlast)

      ! find projected point +1,0 from ilast,jlast
      call apply_periodic_condition_search(ilast, jlast, 1, 0, inext, jnext)
      lonin = longrd(inext,jnext); latin = latgrd(inext,jnext)
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, x1, y1)

      ! find projected point 0,+1 from ilast,jlast
      call apply_periodic_condition_search(ilast, jlast, 0, 1, inext, jnext)
      lonin = longrd(inext,jnext); latin = latgrd(inext,jnext)
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, x2, y2)

      ! the projected target point relative to ilast, jlast
      lonin = lonw; latin = latw
      call ll_to_polars(oned, latwd, lonwd, latin, lonin, xg, yg)

      ! compute strides in x
!--------------------------------------------------
! project point onto local axes (x0,x1,x2)
!--- vectorial equation
! X0G = a X0X1 + b X0X2
!-apply the right cross product with x X0X2 and divide by X0X1 x X0X2
!--- leads to
! x' = a
!--------------------------------------------------

      m012 = y2 * x1 - y1 * x2
      if ( ABS(m012) < epsilon ) then
         n_series = 0
         exit ! likely a grid singularity, we do not want to interpolate close to it, exiting
      endif
      x0g2 = y2 * xg - yg * x2
      a = x0g2 / m012

      ! compute strides in y
!--------------------------------------------------
!-apply the left cross product with X0X1 x and divide by X0X1 x X0X2
!--- leads to
! x" = b
!--------------------------------------------------

      x01g = yg * x1 - y1 * xg
      b = x01g / m012

      di = nint(a); dj = nint(b)

      call apply_periodic_condition_search(ilast,jlast,di,dj,inext,jnext)

      ! preparation for next iteration
      d0 = d
      ilast0 = ilast
      jlast0 = jlast
      d = abs(di) + abs(dj)
      ilast = inext
      jlast = jnext

      counter = counter + 1
    enddo

    ! prepare list of closest cells
    !!! need to rethink this in terms of periodicity !!!
    ilast = ilast0; jlast = jlast0

    ! first search is closest point
    list_isrc(1) = ilast
    list_jsrc(1) = jlast

    ! second search is closest point to the right or left
    if (a<0.5_8) then 
      di = -1
    else
      di =  1
    endif
    call apply_periodic_condition_search(ilast,jlast,di,0,inext,jnext)
    list_isrc(2) = inext
    list_jsrc(2) = jnext

    ! third search is closest point on top or bottom
    if (b<0.5_8) then 
      dj = -1
    else
      dj =  1
    endif
    call apply_periodic_condition_search(ilast,jlast,0,dj,inext,jnext)
    list_isrc(3) = inext
    list_jsrc(3) = jnext

    ! fourth search is closest point in both directions
    call apply_periodic_condition_search(ilast,jlast,di,dj,inext,jnext)
    list_isrc(4) = inext
    list_jsrc(4) = jnext

    END SUBROUTINE refine_search

 
    SUBROUTINE apply_periodic_condition_search(ilast,jlast,di,dj,inext,jnext)
    ! -----------------------------------------------------------------------------
    ! apply periodic conditions on global grids when searching for closest points
    ! -----------------------------------------------------------------------------
    ! arguments
    integer, intent(in) :: ilast, & ! i-index of the starting point
                           jlast, & ! j-index of the starting point
                           di, dj   ! jump stride in i/j
    integer, intent(out) :: inext,jnext ! final position after jump
    ! locals
    integer stri, strj

    ! start with west-east
    inext = ilast + di ! apply stride
    if (webdy) then ! w-e periodicity active
      if (inext < ipl ) then
         stri = di - ipl + ilast + 1  ! remaining stride from right border
         inext = ipr + stri
      endif
      if (inext > ipr ) then
         stri = di - ipr + ilast - 1  ! remaining stride from left border
         inext = ipl + stri
      endif
    else ! simple case, no periodicity, force index to stay within grid
      if (inext >= nxgrd ) inext = nxgrd - 1
      if (inext <      1 ) inext = 1
    endif

    ! y-axis north fold
    jnext = jlast + dj
    select case(nbdy)
    case('T') ! T-folding
      if (jnext > jpt) then
         strj = dj - jpt + jlast  ! remaining stride from top border
         jnext = jpt - strj
         inext = nxgrd - inext + 2
      endif
    case('F') ! F-folding
      if (jnext > jpt) then
         strj = dj - jpt + jlast - 1  ! remaining stride from top border
         jnext = jpt - strj
         inext = nxgrd - inext + 1
      endif
    case default
      if (jnext >= nygrd ) jnext = nygrd - 1
    end select

    if (jnext <      1 ) jnext = 1 ! south boundary

    END SUBROUTINE apply_periodic_condition_search

END MODULE nemo_proj
