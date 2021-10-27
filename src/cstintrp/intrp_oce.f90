MODULE intrp_oce
   USE mympp
   use grids

  !!-----------------------------------------------------------------------
  !!                 ***  MODULE intrp_oce  ***
  !!
  !!  **  Purpose: interpolation tools, F. Roy
  !!
  !!  **  Method: To have fun ...
  !!
  !! history :
  !!     Original code  :   F. Roy (March 2009-November 2009)
  !!-----------------------------------------------------------------------
  !!--------------------------------------------------------------
  !!

    IMPLICIT NONE
    ! * Arguments

    PUBLIC interp_gen_2D,ll_to_polars,inpoly_complex,mapospread,cintrp,fillm

    CONTAINS

!---------------------------------------------------------------------------------------------
      SUBROUTINE interp_gen_2D(x1,y1,x2,y2,x3,y3,x4,y4,x,y,a,b)
!---------------------------------------------------------------------------------------------

! General 2D bilinear interpolation for arbitrary 4 sides polygon
!
!                ------------X3
!      X4---------
!                |
!               |
!              |
!             |
!            |
!           G (point to interpolate to)
!          |
!    X1----------------X2
!
!
      implicit none
! Inputs
      real(kind=8), intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4,x,y ! orientation 1-2-3-4 anticlockwise forms the quadrangle
      real(kind=8), intent(out) :: a, b
! Work
      real*8 m124,x1g4,x134,x12g,x123,den,d,e,ax,bx,cx,dx,u1,u2,u
      logical :: debug, calc

      debug=.false.
      calc = .false.

!--------------------------------------------------
! project point and quadrangle onto local axes
!--- vectorial equation
! X1G = a(1-b) X1X2 + (1-a)b X1X4 + ab X1X3
!-apply the right cross product with x X1X4 and divide by X1X2 x X1X4
!--- leads to
! x' = a(1-b) + ab x3'
!--------------------------------------------------

      m124 = (y4-y1) * (x2-x1) - (y2-y1) * (x4-x1)
      x1g4 = (y4-y1) * (x -x1) - (y -y1) * (x4-x1)
      x134 = (y4-y1) * (x3-x1) - (y3-y1) * (x4-x1)
      if ( ABS(m124) < epsilon ) then
         a = -9.; b = -9.
         write(*,*) 'degenerate quadrangle. exiting'
         return
      endif
      x1g4 = x1g4 / m124
      x134 = x134 / m124

!--------------------------------------------------
!-apply the left cross product with X1X2 x and divide by X1X2 x X1X4
!--- leads to
! x" = (1-a)b + ab x3"
!--------------------------------------------------

      x12g = (y -y1) * (x2-x1) - (y2-y1) * (x -x1)
      x123 = (y3-y1) * (x2-x1) - (y2-y1) * (x3-x1)
      x12g = x12g / m124
      x123 = x123 / m124
     
!--------------------------------------------------
! consider u=ab
! x' = a + u (x3'-1) = a + du
! x" = b + u (x3"-1) = b + eu
!--------------------------------------------------

!--------------------------------------------------
! test for degenerate position of X3 (parallel to X1X2 or X1X4)
!--------------------------------------------------

     d = x134 - oned; e = x123 - oned

     if ( ABS(d) < epsilon ) then
! a=x' ; u = x'b ; x" = b + x' be
        a = x1g4
        den = oned + e * x1g4
        b = x12g / den
        calc = .true.
     endif

     if ( ABS(e) < epsilon .and. .not.calc) then
! b=x" ; u = ax" ; x' = a + adx"
        b = x12g
        den = oned + d * x12g
        a = x1g4 / den
        calc = .true.
     endif

!--------------------------------------------------
! general case
!--------------------------------------------------
! [x'-ud][x"-eu] = ab = u
!--- solve for
! ed u^2 - (1+x'e+x"d) u +x'x" = 0
! with u in [0:1]
!--------------------------------------------------

     if ( .not.calc ) then
        ax = e*d; bx = oned + x1g4 * e + x12g * d; cx = x1g4*x12g
        dx = sqrt( bx**2 - 4.0_8 * ax * cx )
        u1 = ( bx - dx ) / ax * 0.5_8
        u2 = ( bx + dx ) / ax * 0.5_8
        if ( u1 < zerd ) u = u2
        if ( u1 > oned ) u = u2
        if ( u2 < zerd ) u = u1
        if ( u2 > oned ) u = u1
!--------------------------------------------------
! and finally
! a = x' - du
! b = x" - eu
!--------------------------------------------------
        a = x1g4 - d * u 
        b = x12g - e * u
     endif

!--------------------------------------------------
! add bound clipping in case of truncation error
! (should be very small with double precision
! and calculates weights for X1,2,3,4
!--------------------------------------------------
      a = max(min(a, oned), zerd)
      b = max(min(b, oned), zerd)

      return
      END SUBROUTINE interp_gen_2D    


!---------------------------------------------------------------------------------------------
      SUBROUTINE ll_to_polars(r,latc,lonc,lat,lon,xp,yp)
!---------------------------------------------------------------------------------------------
      implicit none

! Lat-Lon to polar orthographic projection with
! center at latc,lonc

! Input/Outputs
      real*8 r,latc,lonc,lat,lon,xp,yp

! ... here latc and lonc and lat and lon are assumed to be degrees
! ... negative west of greenwich and south of equator

! ... latc,lonc are central coordinates
! ... r is the sphere radius
! ... k is the scale factor in both xp,yp directions

! Local variables
      real*8 lambda,lambda0,phi,phi1,k
      real(kind=8), parameter :: pi = atan(1.0_8) * 4.0_8, a=pi/180.0_8

      phi1=latc*a
      phi =lat*a
      lambda0=lonc*a
      lambda =lon*a

!     if stereographic:
!      k=2.d0/( 1.d0+sin(phi1)*sin(phi)    &
!     &             +cos(phi1)*cos(phi)*cos(lambda-lambda0) )

! the orthographic projection is actually easier because k=1
! the equations can be derived by first projecting to cartesian (x,y,z)
! and then applying a rotation that places at the pole the central location
      k = 1.0_8

      xp=r*k*cos(phi)*sin(lambda-lambda0)
      yp=r*k*( cos(phi1)*sin(phi)    &
     &        -sin(phi1)*cos(phi)*cos(lambda-lambda0) )

      return
      END SUBROUTINE ll_to_polars

!-------------------------------------------------------------------------------
      FUNCTION inpoly_complex(x,y,n,xx,yy)
!-------------------------------------------------------------------------------
! Determines if a point xx,yy is inside a complex polygon
! defined by x(n),y(n), the points must be correctly ordered anticloskwise
!-------------------------------------------------------------------------------
      implicit none
! arguments
      integer n
      real*8 x(n),y(n),xx,yy 
      logical inpoly_complex
! Local variables
      integer ncross,i1,i2
      real*8  dxa, dya, dxb, dyb, vecp

      inpoly_complex=.false.

! Travel around segments
      ncross=0
      do i1=1,n        
        i2=mod(i1,n)+1
        dxa = x(i2)-x(i1)
        dya = y(i2)-y(i1)
        dxb = xx   -x(i1)
        dyb = yy   -y(i1)
        vecp = dyb*dxa - dya*dxb
        if (vecp > - epsilon ) then
          ncross = ncross + 1
        endif
      enddo

! Now count nodes
      if (ncross == n) then
        inpoly_complex=.true.
      endif
   
      END FUNCTION inpoly_complex

      SUBROUTINE mapospread( grd, spr, M, nspread, use_local_geo, nearst )
      implicit none
      type(type_grid),    intent(in) :: grd
      type(type_weights), intent(inout) :: spr
      real(kind=8), dimension(:,:), intent(in) :: M
      integer, intent(in) :: nspread
      logical, intent(in) :: use_local_geo, nearst

! Spread data d around unmasked area
! If use_local_geo=.true. apply a weight average according to surrounding lat-lon's
! Assumes the mask is binary (0 or 1)
! Loop around nspread times

! nwgts should be set to 4*nspread, safest choice
! Theoretical max for island shape like this around 'x' for spread:
!    1    2      3        4
!  4*1   4*2    4*3      4x4
                          !
                 !       ! !
          !     ! !     !   !
     !   ! !   !   !   !     !
    !x! ! x ! !  x  ! !   x   !
     !   ! !   !   !   !     !
          !     ! !     !   !
                 !       ! !
                          !

! Local variables
      integer nwgts
      real(kind=4), dimension(:,:),   pointer :: lat,lon
      real(kind=8), dimension(:,:,:), pointer :: wgts
      integer,      dimension(:,:,:), pointer :: iwgts,jwgts
      real(kind=8) :: d(4),dtot12,dtot34,dtot   !Distances
      real(kind=8) :: mw
      real(kind=8) :: xp(4),yp(4)   !Projected lat-lon
      real(kind=8), parameter :: r = oned

      real(kind=8), allocatable, dimension(:,:) :: Ms,Mtmp  ! Spread work mask
      integer i,j,is,iw,iiw,nw,nwadd,iwgts_save,jwgts_save,inrs,jnrs
      real(kind=8) wgts_save
      logical first
      real(kind=8) lond,latd,lond1,latd1
      integer isp
      integer, dimension(4) :: ispp = (/-1, 1, 0, 0 /) ! i-array of adjacent sides
      integer, dimension(4) :: jspp = (/ 0, 0,-1, 1 /) ! j-array of adjacent sides
      real(kind=8), allocatable, dimension(:) :: wa
      integer, allocatable, dimension(:) :: iwa, jwa
      integer, dimension(:,:), pointer :: nws
      integer inrs2, jnrs2
      real(kind=8) wgts_save2
      integer nx,ny

      ! allocate the main arrays and some auxiliaries
      nx = grd % nx
      ny = grd % ny
      nwgts = spr % numwgt

      if (associated(spr % wgt)) deallocate( spr % wgt, spr % iwgt, spr % jwgt, spr % nwgt, spr % mwgt)
      allocate(spr %  wgt(nwgts,nx,ny), &
     &         spr % iwgt(nwgts,nx,ny), &
     &         spr % jwgt(nwgts,nx,ny), &
     &         spr % nwgt(nx,ny), &
     &         spr % mwgt(nx,ny))
      allocate(Ms(nx,ny),Mtmp(nx,ny))

      write(*,*) 'CSTINTRP: ETALEMENT DES DONNEES SUR LES REGIONS MASQUEES nsprd:',nspread
      write(*,*) 'MAPOSPREAD: NUMBER OF WEIGHTS FOR SPREADING=',nwgts
      nw = 4 * nspread * 4
      allocate(wa(nw),iwa(nw),jwa(nw))

      ! default definitions
      Ms=M
      Mtmp=M
       wgts => spr %  wgt;  wgts = 0.
      iwgts => spr % iwgt; iwgts = 1
      jwgts => spr % jwgt; jwgts = 1
        nws => spr % nwgt;   nws = 0
      lon => grd % lld % lont
      lat => grd % lld % latt

      wgts (1,:,:)=1.
      DO j=1,ny
      DO i=1,nx
        iwgts(1,i,j)=i
        jwgts(1,i,j)=j
      ENDDO
      ENDDO
      !
      ! start spreading iteratively
      !
      DO is=1,nspread
        !... defines the weights for the closest points around i,j using inverse distance
        DO j=2,ny-1
        DO i=2,nx-1
          mw=Ms(i-1,j)+Ms(i+1,j)+Ms(i,j-1)+Ms(i,j+1)
          if (Ms(i,j).lt.0.5_8.and.mw.gt.0.5_8) then
            if (nearst) then ! why nearst would use a different spreading technique? Moreover the technique is biased towards the last wet point found
              do isp = 1,4
                inrs = ispp(isp) + i
                jnrs = jspp(isp) + j
                if (Ms(inrs,jnrs).gt.0.5_8) then
                  wgts (1,i,j) = 1.
                  iwgts(1,i,j) = inrs
                  jwgts(1,i,j) = jnrs
                  nws  (  i,j) = 1
                endif
              enddo
            else
              if (use_local_geo) then
                ! conversion to double precision
                lond = lon(i,j) ; latd = lat(i,j)
                do isp = 1,4
                  inrs = ispp(isp) + i ! index of adjacent point
                  jnrs = jspp(isp) + j
                ! ... the center lat-lon is at x,y=0,0
                  lond1 = lon(inrs,jnrs); latd1 = lat(inrs,jnrs)
                  call ll_to_polars(r,latd,lond,latd1,lond1,xp(isp),yp(isp))
                  d(isp) = sqrt((xp(isp)**2+yp(isp)**2))
                  d(isp) = max( d(isp), 1e-7_8 ) ! because some points are actually co-located on the ORCA grid
                  d(isp) = r / d(isp) * Ms(inrs,jnrs) ! use the inverse distance metric for weight => the closer, the higher the weight
                enddo
              else
                do isp = 1,4
                  inrs = ispp(isp) + i ! index of adjacent point
                  jnrs = jspp(isp) + j
                  d(isp) = Ms(inrs,jnrs)
                enddo
              endif
              dtot = SUM( d )
              nw = 0
              do isp = 1,4
                if (d(isp) > 0.0_8) then
                  nw = nw + 1
                  inrs = ispp(isp) + i ! index of adjacent point
                  jnrs = jspp(isp) + j
                   wgts(nw,i,j) = d(isp) / dtot
                  iwgts(nw,i,j) = inrs
                  jwgts(nw,i,j) = jnrs
                endif
              enddo
              nws(i,j) = nw
            endif
            Mtmp(i,j)=1.
          endif
        ENDDO
        ENDDO
! ... Fix weights if is > 1
! ... Goal here is to trace back points that were precedingly used for spread
! ... so that all weights correspond to wet points
        if (is.gt.1) then
          DO j=2,ny-1
          DO i=2,nx-1
            mw=Ms(i-1,j)+Ms(i+1,j)+Ms(i,j-1)+Ms(i,j+1)
            if (Ms(i,j).lt.0.5_8.and.mw.gt.0.5_8) then
              nw = nws(i,j)
              iwa(1:nw) = iwgts(1:nw,i,j); iwgts(1:nw,i,j) = 1
              jwa(1:nw) = jwgts(1:nw,i,j); jwgts(1:nw,i,j) = 1
               wa(1:nw) =  wgts(1:nw,i,j);  wgts(1:nw,i,j) = 0.0_8
              DO iw = 1, nws(i,j)
                inrs = iwa(iw)
                jnrs = jwa(iw)
                if (M(inrs,jnrs).lt.0.5_8) then !Dry point weight, need to find a substitute
                  wgts_save = wa(iw)
                  iwa(iw) = 1
                  jwa(iw) = 1
                   wa(iw) = 0.0_8
                  DO iiw = 1, nws(inrs,jnrs)
                    inrs2 = iwgts(iiw,inrs,jnrs)
                    jnrs2 = jwgts(iiw,inrs,jnrs)
                    wgts_save2 = wgts(iiw,inrs,jnrs)
                    nw = nw + 1
                    iwa(nw) = inrs2
                    jwa(nw) = jnrs2
                     wa(nw) = wgts_save2 * wgts_save
                  ENDDO
                endif
              ENDDO
              ! ... eliminate redundancy
              DO iw= 1, nw-1
                inrs = iwa(iw)
                jnrs = jwa(iw)
                DO iiw = iw+1, nw
                  ! find the point if already stored in wa
                  inrs2 = iwa(iiw)
                  jnrs2 = jwa(iiw)
                  if (inrs2 == inrs .and. jnrs2 == jnrs .and. wa(iiw) > 0.0_8) then
                        wa(iw) = wa(iw) + wa(iiw)
                        wa(iiw) = 0.0_8
                        iwa(iiw) = 1
                        jwa(iiw) = 1
                   endif
                ENDDO
              ENDDO
              ! ... compaction of the resulting array
              iw = 1
              DO WHILE (iw <= nw )
                if (wa(iw) == 0.0_8) then
                  iwa(iw:nw-1) = iwa(iw+1:nw)
                  jwa(iw:nw-1) = jwa(iw+1:nw)
                   wa(iw:nw-1) =  wa(iw+1:nw)
                  nw = nw - 1
                else
                  iw = iw + 1
                endif
              ENDDO
              if (nw.gt.nwgts) then
                write(*,*) 'i,j,nw=',i,j,nw,nwgts
                write(*,*) 'nw too large'
                call my_abort
              endif
              ! copy back to main array
              iwgts(1:nw,i,j) = iwa(1:nw)
              jwgts(1:nw,i,j) = jwa(1:nw)
               wgts(1:nw,i,j) =  wa(1:nw)
               nws(i,j) = nw
            endif
          ENDDO
          ENDDO
        endif
        Ms=Mtmp
      ENDDO   !is=1,nspread

! ... Apply zero gradient at corners and extremities if needed
      i=1
      DO j=1,ny
        if (Ms(i,j).lt.0.5) then
          Mtmp(i,j)=Ms(i+1,j)
            DO iw=1,nwgts
              if (Ms(i+1,j).gt.0.5) then
                iwgts(iw,i,j)=iwgts(iw,i+1,j)
                jwgts(iw,i,j)=jwgts(iw,i+1,j)
                wgts (iw,i,j)=wgts (iw,i+1,j)
              endif
            ENDDO
        endif
      ENDDO
      Ms=Mtmp
      i=nx
      DO j=1,ny
        if (Ms(i,j).lt.0.5) then
          Mtmp(i,j)=Ms(i-1,j)
            DO iw=1,nwgts
              if (Ms(i-1,j).gt.0.5) then
                iwgts(iw,i,j)=iwgts(iw,i-1,j)
                jwgts(iw,i,j)=jwgts(iw,i-1,j)
                wgts (iw,i,j)=wgts (iw,i-1,j)
              endif
            ENDDO
        endif
      ENDDO
      Ms=Mtmp
      j=1
      DO i=1,nx
        if (Ms(i,j).lt.0.5) then
          Mtmp(i,j)=Ms(i,j+1)
            DO iw=1,nwgts
              if (Ms(i,j+1).gt.0.5) then
                iwgts(iw,i,j)=iwgts(iw,i,j+1)
                jwgts(iw,i,j)=jwgts(iw,i,j+1)
                wgts (iw,i,j)=wgts (iw,i,j+1)
              endif
            ENDDO
        endif
      ENDDO
      Ms=Mtmp
      j=ny
      DO i=1,nx
        if (Ms(i,j).lt.0.5) then
            DO iw=1,nwgts
              if (Ms(i,j-1).gt.0.5) then
                iwgts(iw,i,j)=iwgts(iw,i,j-1)
                jwgts(iw,i,j)=jwgts(iw,i,j-1)
                wgts (iw,i,j)=wgts (iw,i,j-1)
              endif
            ENDDO
        endif
      ENDDO

      ! reset nwgts based on the compacted array
      nwgts = maxval(nws)
      spr % numwgt = nwgts

      deallocate(Ms,Mtmp)
      deallocate(iwa,jwa,wa)
      write(*,*) 'CSTINTRP: MAX NUMBER OF EXTRA WEIGHTS AFTER SPREADING = ', nwgts

      END SUBROUTINE mapospread

      SUBROUTINE cintrp( loc, intr, bicub, nxsrc, nysrc )
      implicit none
      ! arguments
      type(type_location), intent(in), target :: loc
      type(type_weights), intent(inout), target :: intr
      logical,         intent(in) :: bicub
      integer,         intent(in) :: nxsrc, nysrc
      ! locals
      integer nxdst,nydst
      real(kind=8),    dimension(:,:,:), pointer :: wgt
      integer(kind=4), dimension(:,:,:), pointer :: iwgt,jwgt
      integer,         dimension(:,:)  , pointer :: II,JJ
      real(kind=8),    dimension(:,:)  , pointer :: AA, BB
      logical,         dimension(:,:)  , pointer :: ML
      integer,         dimension(:,:)  , pointer :: nwgt

! knowing that Idst,Jdst correspond to localisation indexes 
! (linking source grid to destination grid)
! If defrm=.true. consider lat-lon grid deformation
! by projecting on a locally centered PS grid
! Exclude computation over masked data (Mdst).
! Bi-Linear, default if bicub=.false.
! Bi-cubique, if bicub=.true.

      integer i,j,k,is,js
      real(kind=8) :: x1,x2,x3,x4,y1,y2,y3,y4,x,y  !Projected lat-lon
      real(kind=8), parameter :: r=1.0_8
      real(kind=8) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4
      real(kind=8) :: one=1.,two=2.,ov6=1./6.
      real(kind=8) lond,latd,lond1,latd1
     
      II => loc % ig; JJ => loc % jg; ML => loc % mg
      AA => loc % ag; BB => loc % bg
       wgt => intr %  wgt;  wgt=0.
      iwgt => intr % iwgt; iwgt=1  !To avoid index overflow later
      jwgt => intr % jwgt; jwgt=1
      nwgt => intr % nwgt; nwgt=0
      nxdst = intr % nx
      nydst = intr % ny

      DO j=1,nydst
      DO i=1,nxdst
        if( ML(i,j))  then
          is=II(i,j)
          js=JJ(i,j)
          if (is.lt.1) then
            write(*,*) 'Idst too small',II(i,j),i,j
            call my_abort
          endif
          if (js.lt.1) then
            write(*,*) 'Jdst too small',JJ(i,j),i,j
            call my_abort
          endif
          if (is.ge.nxsrc) then
            write(*,*) 'Idst too large',II(i,j),i,j
            call my_abort
          endif
          if (js.ge.nysrc) then
            write(*,*) 'Jdst too large',JJ(i,j),i,j
            call my_abort
          endif
          if (bicub) then
            c3 = AA(i,j) ! a, or the fractional position in i
            d3 = BB(i,j) ! b, or the fractional position in j
            c1 = one - c3
            c0 = - ov6 * c1 * c3
            c2 = c0 * ( two - c3 )
            c4 = c0 * ( one + c3 )
           
            d1 = one - d3
            d0 = - ov6 * d1 * d3
            d2 = d0 * ( two - d3  )
            d4 = d0 * ( one + d3 )
           
            wgt ( 1,i,j) = d1*(c1-2.*c2+c4) - 2.*d2*c1 + d4*c1 ! f_i,j
            iwgt( 1,i,j) = is
            jwgt( 1,i,j) = js
            wgt ( 2,i,j) = d3*(c1-2.*c2+c4) - 2.*d4*c1 + d2*c1 ! f_i,j+1
            iwgt( 2,i,j) = is
            jwgt( 2,i,j) = js+1
            wgt ( 3,i,j) = d3*(c3+c2-2.*c4) + d2*c3 - 2.*d4*c3 ! f_i+1,j+1
            iwgt( 3,i,j) = is+1
            jwgt( 3,i,j) = js+1
            wgt ( 4,i,j) = d1*(c3+c2-2.*c4) - 2.*d2*c3 + d4*c3 ! f_i+1,j
            iwgt( 4,i,j) = is+1
            jwgt( 4,i,j) = js
            wgt ( 5,i,j) = d1*c2 ! f_i-1,j
            iwgt( 5,i,j) = MAX(is-1,1)
            jwgt( 5,i,j) = js
            wgt ( 6,i,j) = d2*c1 !f_i,j-1
            iwgt( 6,i,j) = is
            jwgt( 6,i,j) = MAX(js-1,1)
            wgt ( 7,i,j) = c3*d2 ! f_i+1,j-1
            iwgt( 7,i,j) = is+1
            jwgt( 7,i,j) = MAX(js-1,1)
            wgt ( 8,i,j) = d3*c2 ! f_i-1,j+1
            iwgt( 8,i,j) = MAX(is-1,1)
            jwgt( 8,i,j) = js+1
            wgt ( 9,i,j) = d1*c4 ! f_i+2,j
            iwgt( 9,i,j) = MIN(is+2,nxsrc)
            jwgt( 9,i,j) = js
            wgt (10,i,j) = d3*c4 ! f_i+2,j+1
            iwgt(10,i,j) = MIN(is+2,nxsrc)
            jwgt(10,i,j) = js+1
            wgt (11,i,j) = d4*c1 ! f_i,j+2
            iwgt(11,i,j) = is
            jwgt(11,i,j) = MIN(js+2,nysrc)
            wgt (12,i,j) = d4*c3 ! f_i+1,j+2
            iwgt(12,i,j) = is+1
            jwgt(12,i,j) = MIN(js+2,nysrc)
            nwgt(   i,j) = 12
            ! WARNING
          else
            !linear interpolation default (bicub=.false.)
            wgt (1,i,j)= (1._8-AA(i,j)) * (1._8-BB(i,j))
            iwgt(1,i,j)=is
            jwgt(1,i,j)=js
            wgt (2,i,j)=       AA(i,j)  * (1._8-BB(i,j))
            iwgt(2,i,j)=is+1
            jwgt(2,i,j)=js
            wgt (3,i,j)=       AA(i,j)  *       BB(i,j)
            iwgt(3,i,j)=is+1
            jwgt(3,i,j)=js+1
            wgt (4,i,j)= (1._8-AA(i,j)) *       BB(i,j)
            iwgt(4,i,j)=is
            jwgt(4,i,j)=js+1
            nwgt(  i,j)=4
            ! interpolate
          endif

        endif
      ENDDO
      ENDDO

      END SUBROUTINE cintrp

      SUBROUTINE caggrt( coli, grds, grdr, intr, defrm, nadis, Ms)
      implicit none
      ! arguments
      type(type_location), intent(in), target :: coli ! reverse location
      type(type_grid)    , intent(in), target :: grds, grdr
      type(type_weights), intent(inout), target :: intr
      real(kind=8),dimension(:,:), intent(in) :: Ms
      logical, intent(in) :: defrm
      integer, intent(in) :: nadis
      ! locals
      integer nxdst,nydst, nxsrc, nysrc
      real(kind=8),    dimension(:,:,:), pointer :: wgt
      integer(kind=4), dimension(:,:,:), pointer :: iwgt,jwgt
      integer,         dimension(:,:)  , pointer :: II,JJ
      real(kind=8),    dimension(:,:)  , pointer :: AA, BB
      logical,         dimension(:,:)  , pointer :: ML, Mdsta
      integer, dimension(:,:), pointer :: Nadst
      real(kind=4),    dimension(:,:), pointer :: latsrc,lonsrc
      real(kind=4),    dimension(:,:), pointer :: latdst,londst

! Return mask where no data on destination grid (Mdsta)
! No computation on source mask Msrc

      integer :: naggrmax, naggrtot
      integer i,j,k,is,js,icl,jcl,nwgtw,icldm,jcldm
      logical, dimension(:,:)  , allocatable :: Msrc
      real(kind=8), allocatable, dimension(:,:)   :: surtot 
      real(kind=8) :: surwrk
      real(kind=8), parameter :: r=6371.0_8 !Earth avg radius
      real(kind=8) :: dxi,dyi,dxj,dyj,dlati,dlatj,dloni,dlonj, &
     &                lati,latj,loni,lonj
      integer nmixt, naggr, nwarn
      logical :: print_warning, printw

      II => coli % ig; JJ => coli % jg; ML => coli % mg
      AA => coli % ag; BB => coli % bg
       wgt => intr %  wgt; wgt=0.
      iwgt => intr % iwgt; iwgt=1  !To avoid index overflow later
      jwgt => intr % jwgt; jwgt=1
      Mdsta=> intr % mwgt; Mdsta=.false.
      Nadst=> intr % nwgt; Nadst=0
      nxdst = intr % nx
      nydst = intr % ny
      nxsrc = SIZE( Ms, 1)
      nysrc = SIZE( Ms, 2)
      naggrmax = intr % numwgt
      lonsrc => grds % lld % lont
      latsrc => grds % lld % latt
      londst => grdr % lld % lont
      latdst => grdr % lld % latt

      ! Action

      allocate(surtot(nxdst,nydst)); surtot=0.
      allocate( Msrc(nxsrc,nysrc) ); Msrc  = .false.; where( Ms > 0.5_8) Msrc = .true.

      call remove_periodic_points(grds, Msrc)
      where( .not.ML ) Msrc = .false. ! where the destination grid does not match the source grid

      ! ... First compile closest point in 3D table
      naggrtot=0
      nwgtw=max(nadis-1,0)
      print_warning=.false.
      nwarn=0

! Compute weighting surfaces
! from closest point in defined neighbour
      do js=1,nysrc
      do is=1,nxsrc
        if ( Msrc(is,js) ) then
         icl=II(is,js) + nint(AA(is,js))
         jcl=JJ(is,js) + nint(BB(is,js))
         do j=jcl-nwgtw,jcl+nwgtw
         do i=icl-nwgtw,icl+nwgtw
           if (i.ge.1.and.i.le.nxdst.and.j.ge.1.and.j.le.nydst) then
             if (defrm) then ! Compute weight surface for each point
                             ! A relatively precise estimate
                             ! Neglecting surface overflows from
                             ! destination grid point
               if (is.gt.1.and.is.lt.nxsrc) then
                 loni= 0.5*(lonsrc(is,js)  +lonsrc(is-1,js))
                 dloni=0.5*(lonsrc(is+1,js)-lonsrc(is-1,js))
                 lati= 0.5*(latsrc(is,js)  +latsrc(is-1,js))
                 dlati=0.5*(latsrc(is+1,js)-latsrc(is-1,js))
                elseif (is.eq.1) then
                 loni= lonsrc(is,js)-0.5*(lonsrc(is+1,js)-lonsrc(is,js))
                 dloni=lonsrc(is+1,js)-lonsrc(is,js)
                 lati= latsrc(is,js)-0.5*(latsrc(is+1,js)-latsrc(is,js))
                 dlati=latsrc(is+1,js)-latsrc(is,js)
                elseif (is.eq.nxsrc) then
                 loni= lonsrc(is,js)-0.5*(lonsrc(is,js)-lonsrc(is-1,js))
                 dloni=lonsrc(is,js)-lonsrc(is-1,js)
                 lati= latsrc(is,js)-0.5*(latsrc(is,js)-latsrc(is-1,js))
                 dlati=latsrc(is,js)-latsrc(is-1,js)
               endif
               if (js.gt.1.and.js.lt.nysrc) then
                 lonj= 0.5*(lonsrc(is,js)  +lonsrc(is,js-1))
                 dlonj=0.5*(lonsrc(is,js+1)-lonsrc(is,js-1))
                 latj= 0.5*(latsrc(is,js)  +latsrc(is,js-1))
                 dlatj=0.5*(latsrc(is,js+1)-latsrc(is,js-1))
                elseif (js.eq.1) then
                 lonj= lonsrc(is,js)-0.5*(lonsrc(is,js+1)-lonsrc(is,js))
                 dlonj=lonsrc(is,js+1)-lonsrc(is,js)
                 latj= latsrc(is,js)-0.5*(latsrc(is,js+1)-latsrc(is,js))
                 dlatj=latsrc(is,js+1)-latsrc(is,js)
                elseif (js.eq.nysrc) then
                 lonj= lonsrc(is,js)-0.5*(lonsrc(is,js)-lonsrc(is,js-1))
                 dlonj=lonsrc(is,js)-lonsrc(is,js-1)
                 latj= latsrc(is,js)-0.5*(latsrc(is,js)-latsrc(is,js-1))
                 dlatj=latsrc(is,js)-latsrc(is,js-1)
               endif
               call ll_to_polars(r,lati,loni,    &
     &           lati+dlati,loni+dloni,dxi,dyi)
               call ll_to_polars(r,latj,lonj,    &
     &           latj+dlatj,lonj+dlonj,dxj,dyj)
               surwrk=sqrt( dxi**2 + dyi**2 ) * sqrt( dxj**2 + dyj**2 )
              else
               surwrk=1.  !Equal weighting (simple average)
             endif
             k=Nadst(i,j)+1
             if (k.le.naggrmax) then
               Nadst(i,j)=Nadst(i,j)+1
               if (k.gt.naggrtot) then
	         icldm=i
	         jcldm=j
	       endif
               naggrtot=max(naggrtot,k)
               if (k.eq.naggrmax) then
                 print_warning=.true. 
                 nwarn=nwarn+1
	       endif
               surtot(i,j)=surtot(i,j)+surwrk
               wgt (k,i,j)=surwrk
               iwgt(k,i,j)=is
               jwgt(k,i,j)=js
             endif
           endif
         enddo
         enddo
        endif
      enddo
      enddo
      if (print_warning) then
        write(*,*) 'CSTINTRP: WARNING not enough storage to aggregate all points'
        write(*,*) 'CSTINTRP: FOR THIS NUMBER OF DESTINATION POINTS ==>', nwarn
        write(*,*) 'CSTINTRP: naggrmax=',naggrmax
      endif
   
! Now aggregate over destination grid
      naggr=0
      do j=1,nydst
      do i=1,nxdst
        if (Nadst(i,j).gt.0) then
          naggr = naggr + 1
          Mdsta(i,j) = .true.
          if (surtot(i,j).gt.0.) then
            do k=1,Nadst(i,j)
              wgt(k,i,j)=wgt(k,i,j)/surtot(i,j)
            enddo
          else !case where all weights are zero (equal weights, no surfaces)
            write(*,*) 'SOMETHING WRONG WITH WEIGHTING AT DESTINATION',i,j
            write(*,*) 'WEIGHTS WERE ALL FIXED TO ZERO, ... HOPING IT IS A DRY POINT'
          endif
        endif
      enddo
      enddo

      write(*,*) 'CSTINTRP: Local max number of points in aggregation, naggrtot=',naggrtot
      write(*,*) 'CSTINTRP: Local max number of points at i,j=',icldm,jcldm
      write(*,*) 'CSTINTRP: Total number of aggregated points=',naggr

      if (.not.defrm) then
        write(*,*) 'CSTINTRP: Local lat-lon deformation neglected'
      endif

      deallocate(surtot, Msrc)
      intr % numwgt = naggrtot

      END SUBROUTINE caggrt


  SUBROUTINE adjust_wgtb( pintr, psprd, Ms )

    implicit none
    ! arguments
    type(type_weights), intent(inout), target :: pintr
    type(type_weights), intent(in)   , target :: psprd
    real(kind=8), intent(in)     :: Ms(:,:)
    ! locals
    integer nxdst,nydst, nxsrc, nysrc
    integer nwgtb, nwgts, nwgtbmax
    real(kind=8), dimension(:,:,:), pointer  :: wgtb, wgts
    integer,      dimension(:,:,:), pointer  :: iwgtb, jwgtb, iwgts, jwgts
    integer, dimension(:,:), pointer :: nws
    integer i,j,nwbase, nw, iw, iiw, nnw
    integer inrs, jnrs, inrs2, jnrs2
    real(kind=8) wgt_save, wgt_save2
    integer ib,jb
    real(kind=8), allocatable, dimension(:) :: wa
    integer, allocatable, dimension(:) :: iwa, jwa

    write(*,*) 'CSTINTRP: ADJUSTING WEIGHTS ACCORDING TO SPREAD'

     wgtb => pintr %  wgt
    iwgtb => pintr % iwgt
    jwgtb => pintr % jwgt
      nws => pintr % nwgt
     wgts => psprd %  wgt
    iwgts => psprd % iwgt
    jwgts => psprd % jwgt

    nwbase = SIZE( pintr % wgt, 1 )
    nwgtb  = pintr % numwgt
    nwgtbmax=0
    nxdst = pintr % nx
    nydst = pintr % ny
    nxsrc = SIZE( Ms, 1)
    nysrc = SIZE( Ms, 2)
    nw = 4 * nwbase
    allocate(wa(nw),iwa(nw),jwa(nw))

    DO j=1,nydst
    DO i=1,nxdst
      nw = nwgtb
      iwa(1:nw) = iwgtb(1:nw,i,j); iwgtb(1:nw,i,j) = 1
      jwa(1:nw) = jwgtb(1:nw,i,j); jwgtb(1:nw,i,j) = 1
       wa(1:nw) =  wgtb(1:nw,i,j);  wgtb(1:nw,i,j) = 0.0_8
      DO iw = 1, pintr % numwgt
        inrs = iwa(iw)
        jnrs = jwa(iw)
        nnw = psprd % nwgt(inrs,jnrs)
        if (wa(iw).ne.0.0_8.and.Ms(inrs,jnrs).lt.0.5_8.and.nnw>0) then !Dry point weight, need to find a substitute
          wgt_save = wa(iw)
          iwa(iw) = 1
          jwa(iw) = 1
           wa(iw) = 0.0_8
          DO iiw =1 , nnw
            inrs2 = iwgts(iiw,inrs,jnrs)
            jnrs2 = jwgts(iiw,inrs,jnrs)
            wgt_save2 = wgts(iiw,inrs,jnrs)
            nw = nw + 1
            iwa(nw) = inrs2
            jwa(nw) = jnrs2
             wa(nw) = wgt_save2 * wgt_save
          ENDDO
        endif
      ENDDO
      ! ... eliminate redundancy
      DO iw = 1, nw-1
        inrs = iwa(iw)
        jnrs = jwa(iw)
        DO iiw = iw+1, nw
          ! find the point if already stored in wa
          inrs2 = iwa(iiw)
          jnrs2 = jwa(iiw)
          if (inrs2 == inrs .and. jnrs2 == jnrs .and. wa(iiw) .ne. 0.0_8) then
                wa(iw) = wa(iw) + wa(iiw)
                wa(iiw) = 0.0_8
                iwa(iiw) = 1
                jwa(iiw) = 1
           endif
        ENDDO
      ENDDO
      ! ... compaction of the resulting array
      iw = 1
      DO WHILE (iw <= nw )
        if (wa(iw) == 0.0_8) then
          iwa(iw:nw-1) = iwa(iw+1:nw)
          jwa(iw:nw-1) = jwa(iw+1:nw)
           wa(iw:nw-1) =  wa(iw+1:nw)
          nw = nw - 1
        else
          iw = iw + 1
        endif
      ENDDO
      if (nw.gt.nwbase) then
        write(*,*) 'i,j,nw=',i,j,nw,nwbase
        write(*,*) 'nw too large'
        call my_abort
      endif
      ! copy back to main array
      iwgtb(1:nw,i,j) = iwa(1:nw)
      jwgtb(1:nw,i,j) = jwa(1:nw)
       wgtb(1:nw,i,j) =  wa(1:nw)
       nws(i,j) = nw
      nwgtbmax = max(nwgtbmax,nw)
    ENDDO
    ENDDO

    write(*,*) 'CSTINTRP: MAX NUMBER OF WEIGHTS AFTER SPREAD = ', nwgtbmax
    pintr % numwgt = nwgtbmax
    deallocate(wa,iwa,jwa)

  END SUBROUTINE adjust_wgtb


      SUBROUTINE rotate_vector_ninj(u,v,theta,ni,nj)
!--------------------------------------------------------------------------
! Sous-Routine rotate_vector
! Tourne un vecteur (u,v) dans une base orthonormee ayant un angle
! de theta (dans le sens horaire) par rapport a l'ancienne base
!--------------------------------------------------------------------------

! Dans le modele ROM theta=49 degres (il faut transformer en radians)
! si on va de (x,y) vers (est,nord) theta= 49 deg
! si on va de (est,nord) vers (x,y) theta=-49 deg

      implicit none
      integer ni,nj
      real*8 u(ni,nj),v(ni,nj)
      real theta(ni,nj)
      real, allocatable :: utmp(:,:),vtmp(:,:)

      allocate(utmp(ni,nj),vtmp(ni,nj))
      utmp(:,:) = cos(theta(:,:))*u(:,:) - sin(theta(:,:))*v(:,:)
      vtmp(:,:) = sin(theta(:,:))*u(:,:) + cos(theta(:,:))*v(:,:)

      u = utmp
      v = vtmp

      deallocate(utmp,vtmp)
      return
      END SUBROUTINE rotate_vector_ninj

      SUBROUTINE fillm(var,mskd,mskr,nx,ny,npts,npass)
      ! Fill holes in a field using mask (msk)
      ! based on a pseudo-inhouse kriging method (linear)
      implicit none
      integer, intent(in)                      :: nx,ny
      real, intent(inout),    dimension(nx,ny) :: var
      integer, intent(inout), dimension(nx,ny) :: mskd
      integer, intent(in),    dimension(nx,ny) :: mskr
      integer, intent(in), optional            :: npts,npass

      real,     dimension(nx,ny) :: var_tmp
      integer,  dimension(nx,ny) :: msk_tmp

      integer i,j, is, js, ipass, npoints, nmid, nwgt
      integer ileft,iright,jup,jdown,idleft,idright
      logical left,right,up,down,dleft,dright
      real    val,d1,d2

      real,    allocatable, dimension(:,:) :: subxx
      integer, allocatable, dimension(:,:) :: submr,submd

      integer npasse
      logical debug

      npoints=11
      npasse=5
      if (present(npts))  npoints=npts
      if (present(npass)) npasse=npass

      if (mod(npoints,2) == 0) then
        write(*,*) 'fillm: npoints must not be even'
        call my_abort
      endif

      var_tmp(:,:) = var(:,:)
      msk_tmp(:,:) = mskd(:,:)

      !Add 3 extra passes with reduced number of points to deal with boundaries
      do ipass=1,npasse+3
        nmid=int(real(npoints)/2.)
        allocate(subxx(npoints,npoints),submr(npoints,npoints),submd(npoints,npoints))
        do js=1,ny-npoints+1
        do is=1,nx-npoints+1
          subxx(:,:)=var( is:is+npoints-1,js:js+npoints-1)
          submr(:,:)=mskr(is:is+npoints-1,js:js+npoints-1)
          submd(:,:)=mskd(is:is+npoints-1,js:js+npoints-1)
          if ( submr(nmid,nmid) == 1 .and. submd(nmid,nmid) == 0 ) then
            ! Check horizontal 
            left=.false.
            right=.false.
            j=nmid
            do i=1,nmid-1
              if (submd(i,j)==1) then
                left=.true.
                ileft=i
              endif
            enddo
            do i=npoints,nmid+1,-1
              if (submd(i,j)==1) then
                right=.true.
                iright=i
              endif
            enddo
            ! Check vertical
            up=.false.
            down=.false.
            i=nmid
            do j=1,nmid-1
              if (submd(i,j)==1) then
                down=.true.
                jdown=j
              endif
            enddo
            do j=npoints,nmid+1,-1
              if (submd(i,j)==1) then
                up=.true.
                jup=j
              endif
            enddo
            ! Check diagonal
            dleft=.false.
            dright=.false.
            do i=1,nmid-1
              j=i
              if (submd(i,j)==1) then
                dleft=.true.
                idleft=i
              endif
            enddo  
            do i=npoints,nmid+1,-1
              j=i
              if (submd(i,j)==1) then
                dright=.true.
                idright=i
              endif
            enddo
            nwgt=0
            val=0.
            if ( left .and. right ) then
              nwgt=nwgt+1
              d1=nmid-ileft
              d2=iright-nmid
              val=val+( d1 * subxx(iright,nmid) + d2 * subxx(ileft,nmid) ) / (d1+d2)
            endif
            if ( up .and. down ) then
              nwgt=nwgt+1
              d1=nmid-jdown
              d2=jup-nmid
              val=val+( d1 * subxx(nmid,jup)    + d2 * subxx(nmid,jdown)  ) / (d1+d2)
            endif
            if ( dleft .and. dright ) then
              nwgt=nwgt+1
              d1=nmid-idleft
              d2=idright-nmid
              val=val+( d1 * subxx(idright,idright) + d2 * subxx(idleft,idleft) ) / (d1+d2)
            endif
            if ( nwgt > 0 ) then
              val=val/nwgt
              var_tmp(is+nmid-1,js+nmid-1) =  val
              msk_tmp(is+nmid-1,js+nmid-1) =  1
            endif
          endif
        enddo
        enddo
        var(:,:)  = var_tmp(:,:)
        mskd(:,:) = msk_tmp(:,:)
        deallocate(subxx,submr,submd)
        if ( ipass >= npasse ) then
          npoints = npoints / 2
          if (mod(npoints,2) == 0) npoints = npoints + 1
        endif
      enddo
 
      return
      END SUBROUTINE fillm

END MODULE intrp_oce
