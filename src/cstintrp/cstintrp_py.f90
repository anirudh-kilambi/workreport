module cstintrp_py
    
  ! module declaration
    USE iso_c_binding
    USE cstintrp_mod
    USE intrp_itfc
    USE utils
    USE mympp
    USE std
    USE grids
    USE llxy

    IMPLICIT NONE

    PUBLIC

!! Argument variables:

    CONTAINS

    subroutine cstintrp_from_call(xsrc,ysrc,fsrc,xdst,ydst,fdst) bind(C, name="cstintrp_from_call_c")  
    !------------------------------------------------------------------------------
    ! main call for treating interpolation between src and dst grids`
    !------------------------------------------------------------------------------
    ! arguments
    real(c_double), dimension(:,:) :: xsrc,ysrc,fsrc,xdst,ydst,fdst
    type(type_grid), target &
                        :: F_grds, F_grdr
    ! locals
    type(type_data), target  :: psrc, pdst
    type(std_grid),  pointer :: grdp
    type(type_grid), pointer :: subg, subs, subr 
    type(type_location), target :: plocr, ploci
    type(type_weights), target  :: pintr, psprd
    integer nx, ny
    logical vector
    type(variable)      :: var_grds

    !---------------------------------------------------------------------------------------
    call read_vector_table

    ! Initialize parameters and validate options with source and destination grid

    nx = size(xsrc,1)
    ny = size(xsrc,2)
    grdp => F_grds % std_grd
    grdp % ni = nx
    grdp % nj = ny
    grdp % grtyp = 'O' ! by default assume the grid of being irregular in type
    F_grds % dxtree  = dxred
    F_grds % dytree  = dyred
    F_grds % nbdy    = snbdy
    F_grds % nrepeat = snrep
    F_grds % webdy   = sebdy
    F_grds % werepeat= serep
    F_grds % nx      = nx
    F_grds % ny      = ny
    F_grds % ng      = 1
    allocate( F_grds % subgrd(1) )
    subg => F_grds % subgrd(1)
    subg % std_grd = grdp
    subg % nx      = F_grds % nx
    subg % ny      = F_grds % ny
    subg % dxtree  = F_grds % dxtree
    subg % dytree  = F_grds % dytree
    subg % nbdy    = F_grds % nbdy
    subg % nrepeat = F_grds % nrepeat
    subg % webdy   = F_grds % webdy
    subg % werepeat= F_grds % werepeat
    allocate( subg % lld % latt(nx,ny), &
              subg % lld % lont(nx,ny) )
    subg % lld % lont(:,:) = xsrc(:,:)
    subg % lld % latt(:,:) = ysrc(:,:)


    nx = size(xdst,1)
    ny = size(xdst,2)
    grdp => F_grdr % std_grd
    grdp % ni = nx
    grdp % nj = ny
    grdp % grtyp = 'O' ! by default assume the grid of being irregular in type
    F_grdr % dxtree  = dxred
    F_grdr % dytree  = dyred
    F_grdr % nbdy    = dnbdy
    F_grdr % nrepeat = dnrep
    F_grdr % webdy   = debdy
    F_grdr % werepeat= derep
    F_grdr % nx      = grdp % ni
    F_grdr % ny      = grdp % nj
    F_grdr % ng      = 1
    allocate( F_grdr % subgrd(1) )
    subg => F_grdr % subgrd(1)
    subg % std_grd = grdp
    subg % nx      = F_grdr % nx
    subg % ny      = F_grdr % ny
    subg % dxtree  = F_grdr % dxtree
    subg % dytree  = F_grdr % dytree
    subg % nbdy    = F_grdr % nbdy
    subg % nrepeat = F_grdr % nrepeat
    subg % webdy   = F_grdr % webdy
    subg % werepeat= F_grdr % werepeat
    allocate( subg % lld % latt(nx,ny), &
              subg % lld % lont(nx,ny) )
    subg % lld % lont(:,:) = xdst(:,:)
    subg % lld % latt(:,:) = ydst(:,:)

    !---------------------------------------------------------------------------------------
    ! find the variables associated with that particular source grid
    var_grds % varstd % nomvar = 'NAME'
    var_grds % varstd % typvar ='P'
    var_grds % varstd % ip1 = 0
    var_grds % varstd % npas = 0
    var_grds % ln_mask = .false.
    var_grds % ln_vector = .false.

    call intrp_itfc_params( F_grds, F_grdr )

    ! Allocate memory and copy grid descriptors to destination file
    nxsrc = F_grds % nx
    nysrc = F_grds % ny
    nxdst = F_grdr % nx
    nydst = F_grdr % ny
    if ( F_grds % std_grd % grtyp == 'U' ) nysrc = nysrc / 2
    if ( F_grdr % std_grd % grtyp == 'U' ) nydst = nydst / 2

    allocate(Msrc(nxsrc,nysrc,1))
    allocate(Mdstw(nxdst,nydst,1,1))

! Calcul des xy de la grille source sur la destination

      subs => F_grds % subgrd(1)
      subr => F_grdr % subgrd(1)
      subgsrc(1) = subs % std_grd
      subgdst(1) = subr % std_grd

      ! do the localization of the destination points on the source grid
      if( subs % std_grd % grtyp /= 'Y' .and. &
          (intyp /= "aggr" .and. intyp(1:5) /= "scrip") ) then
        write(*,*) 'CSTINTRP: CALCUL SUR LA GRILLE DESTINATION DES XY DU REPERE SOURCE'
        call ll2xy( subs, subr, plocr )
             ! FD debug
             write(*,*) 'forward localization'
             write(*,*) 'min/max I',minval( plocr % ig(:,:) ),maxval( plocr % ig(:,:) )
             write(*,*) 'min/max J',minval( plocr % jg(:,:) ),maxval( plocr % jg(:,:) )
             write(*,*) 'min/max A',minval( plocr % ag(:,:) ),maxval( plocr % ag(:,:) )
             write(*,*) 'min/max B',minval( plocr % bg(:,:) ),maxval( plocr % bg(:,:) )
             write(*,*) 'min/max a',minval( plocr % pg(:,:) ),maxval( plocr % pg(:,:) )
      endif

      ! do reverse localization used in aggregation
      select case(intyp)
      case( 'aggr', 'mixt', 'mixtbc')
        write(*,*) 'CSTINTRP: CALCUL SUR LA GRILLE SOURCE DES XY DU REPERE DESTINATION'
        call ll2xy( subr, subs, ploci )
             ! FD debug
             write(*,*) 'backward localization'
             write(*,*) 'min/max I',minval( ploci % ig(:,:) ),maxval( ploci % ig(:,:) )
             write(*,*) 'min/max J',minval( ploci % jg(:,:) ),maxval( ploci % jg(:,:) )
             write(*,*) 'min/max A',minval( ploci % ag(:,:) ),maxval( ploci % ag(:,:) )
             write(*,*) 'min/max B',minval( ploci % bg(:,:) ),maxval( ploci % bg(:,:) )
             write(*,*) 'min/max a',minval( ploci % pg(:,:) ),maxval( ploci % pg(:,:) )
      end select

        write(*,*) 'CSTINTRP: FIRST FIELD ... FULL COMPUTATION'

        if ( subs % std_grd % grtyp == 'M') then ! unstructured grids
          SELECT CASE ( intyp )
          CASE('nearst') 
            call intrp_itfc_tri(plocr, pintr, subs, Msrc(:,:,1), nearst=.true.)
          CASE('bilin')
            call intrp_itfc_tri(plocr, pintr, subs, Msrc(:,:,1) )
          CASE('aggr') ! aggr
            call intrp_itfc_aggr(plocr, ploci, pintr, subs, subr, &
Msrc(:,:,1) )
          CASE('scrip0')
            call intrp_itfc_scrip(plocr, ploci, pintr, subs, subr, &
Msrc(:,:,1), .false.)
          END SELECT
        else ! regular grids
          SELECT CASE ( intyp )
          CASE('nearst') 
            call intrp_itfc_nearst(plocr, pintr, psprd, Msrc(:,:,1))
          CASE('bilin','bicub') ! bilin (and possibly nearest) or bicub
            call intrp_itfc_bilc(plocr, pintr, psprd, Msrc(:,:,1))
          CASE('aggr') ! aggr
            call intrp_itfc_aggr(plocr, ploci, pintr, subs, subr, &
Msrc(:,:,1))
          CASE('mixt','mixtbc') ! mixt or mixtbc
            call intrp_itfc_mixtlc(plocr, ploci, pintr, psprd, subs, &
subr, Msrc(:,:,1))
          CASE('scrip0','scrip1')
            call intrp_itfc_scrip(plocr, ploci, pintr, subs, subr, &
Msrc(:,:,1), intyp=='scrip1')
          END SELECT
        endif

        call adjust_periodic(subs, subr, pintr )

             ! FD debug
             write(*,*) 'weights', pintr % numwgt
             write(*,*) 'min/max W',minval( pintr %  wgt ),maxval( pintr %  wgt )
             write(*,*) 'min/max I',minval( pintr % iwgt ),maxval( pintr % iwgt )
             write(*,*) 'min/max J',minval( pintr % jwgt ),maxval( pintr % jwgt )
             write(*,*) 'min/max N',minval( pintr % nwgt ),maxval( pintr % nwgt )
        
! Initialize working mask
        Mdstw(:,:,1,1) = 0
        where( pintr % mwgt(:,:) ) Mdstw(:,:,1,1) = 1 ! ... Mdstw is the one written

!...  Application des poids
      allocate(psrc % u(nxsrc,nysrc)); psrc % u(:,:) = fsrc(:,:)
      allocate(pdst % u(nxdst,nydst))

      call apply_weights( pdst % u, psrc % u, pintr)
             ! FD debug
             write(*,*) 'min/max F',minval(pdst % u ),maxval(pdst % u )
      fdst(:,:) = pdst % u (:,:)

    end subroutine cstintrp_from_call

end module cstintrp_py

