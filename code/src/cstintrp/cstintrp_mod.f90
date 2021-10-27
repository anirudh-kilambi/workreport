module cstintrp_mod
    ! module declaration
    USE intrp_itfc
    USE utils
    USE mympp
    USE std
    USE grids

    IMPLICIT NONE

    PUBLIC

!! Argument variables:

    !--- file management
    character(len=1000) :: fs,fr,fd      !Source,reference,destination files,respectively
    character(len=6)    :: fstype,fdtype !Type of source destination file
    character(len=4)    :: nomvar_src    !Nomvar selection in source file (default nomvar_src=' ')
    character(len=4)    :: nomvar_ref    !Reference field defining destination grid (in fr)
    character(len=1000) :: myvtab
    integer             :: narg
    character(len=1000) :: file_status   !deprecated option

    !... variable management
    type (std_variable),  &
       dimension(:),      &
       allocatable      :: vars,    &  !All variables from source file
                           vars_g,  &  !vars for one particular grid in grds
                           vars_gs     !vars_g sorted by ip1
    integer             :: nvs, nvr, nvs_g, nvr_g
    integer             ::           ivs_g, ivr_g, ivs_g_old

    type (type_grid),     &
       dimension(:),      &
       allocatable      :: grds, & !All grids from source file
                           grdr    !All grids from reference file corresponding to nomvar=nr
    integer             :: ngs, ngr
    integer             :: dxred, dyred
    ! periodicity
    integer :: serep = 0, snrep = 0, derep = 0, dnrep = 0
    logical :: sebdy = .false., debdy = .false.
    character(len=1) :: snbdy = '', dnbdy = ''

    CONTAINS

    subroutine init_general_param
! Definition des arguments par defaut

    nomvar_src=' '
    mtrsh=0.5
    nsprd=4
    nbdyf=0
    sdfrm=.true.
    owgts=.false.
    file_status=''
    nadis=0
    npak_ou=999
    glbsrc=.false.
    naggrmax=50
    uwgdir = 'NONE'
    myvtab = ''
    splt_yy = .false.
    zland  = .false.
    dxred  = 1
    dyred  = 1
    intyp = 'mixt'
    landval = -999.9

    end subroutine init_general_param

    subroutine print_usage
    implicit none
    ! locals
    integer             :: iargc

! Traitement des arguments
    narg= iargc()
    IF ( narg < 8 .or. mod(narg,2).ne.0 ) THEN
       PRINT *, ' Usage : cstintrp '
       PRINT *, ' -fs filesrc (source file to interpolate)'
       PRINT *, ' -fr fileref (reference file defining destination grid)'
       PRINT *, ' -nr nomvar_ref (reference variable defining destination grid, in fr)'
       PRINT *, ' -fd filedst (destination file, interpolated data)'
       PRINT *, ' Optional arguments:'
       PRINT *, '   -fstype fstype (source file type: rpn or custom),'
       PRINT *, '                 deprecated now, does nothing'
       PRINT *, '   -fdtype fdtype (destination file type: rpn or custom),'
       PRINT *, '                 deprecated now, does nothing'
       PRINT *, '   -ns nomvar_src (nomvar selection in src file,default all variables)'
       PRINT *, '   -intyp intyp (regridding type, mixt (default), bilin or aggr:'
       PRINT *, '                 bilin  ==> bi-linear'
       PRINT *, '                 aggr   ==> nearest neighbouring aggregation'
       PRINT *, '                 mixt   ==> combination of bilin and aggr adjusting'
       PRINT *, '                            to local differences in resolution)'
       PRINT *, '                 nearst ==> nearest neighbour (one point)'
       PRINT *, '                 bicub  ==> bi-cubic interpolation'
       PRINT *, '                 mixtbc ==> combination of bicub and aggr adjusting'
       PRINT *, '                            to local differences in resolution)'
       PRINT *, '                 scrip0 ==> SCRIP conservative interpolation order 0'
       PRINT *, '   -nwgt  nwgt  (default=0, aggregation distance parameter'
       PRINT *, '                 for intyp=aggr or intyp=mixt, determines the distance'
       PRINT *, '                 of influence in grid indices of destination grid:'
       PRINT *, '                 if nwgt=0, aggregate from i-0.5=>i+0.5,j-0.5=>j+0.5 '
       PRINT *, '                 if nwgt=1, aggregate from i-1.5=>i+1.5,j-1.5=>j+1.5 '
       PRINT *, '                 if nwgt=2, aggregate from i-2.5=>i+2.5,j-2.5=>j+2.5  and so on)'
       PRINT *, '                 for intyp=scrip*, determines the number of points used in the kd-tree'
       PRINT *, '                 search around the source grid cell (nwgt^2 is the actual number)'
       PRINT *, '                 nwgt needs to be higher than 0 to work and some tests are needed'
       PRINT *, '                 to determine the optimal value for convergence and cost'
       PRINT *, '   -naggrmax  naggrmax  (default=50, max number of aggregation points allowed)'
       PRINT *, '   -mtrsh mtrsh (threshold applied to dst mask, data is masked if mask<mtrsh,default 0.5)'
       PRINT *, '   -nsprd nsprd (no of iter. in spreading near-coast data before interp., ...'
       PRINT *, '                 default => 4 bilin or mixt 6 for bicub and mixtbc)'
       PRINT *, '                         => 6 for bicub and mixtbc)'
       PRINT *, '   -sdfrm sdfrm (T or F, consider local deformation of lat-lon grids,default T)'
       PRINT *, '   -owgts owgts (T or F, output interpolation weights and mask)'
       PRINT *, '   -status file_status (output program status in file_status for op needs)'
       PRINT *, '   -npak npak   (compression parameter for output)'
       PRINT *, '   -glbsrc glb  (T or F, mention that source grid is global, take care of wrap around)'
       PRINT *, '   -uwgdir dir  (Use raw weight file directory, write/read weight files ...'
       PRINT *, '                 corresponding to ip1s, careful to make sure the io grids ... '
       PRINT *, '                 are the same throughout the whole process) '
       PRINT *, '   -myvtab tab  (Use your own version of vector name list file for u, v rotation)'
       PRINT *, '   -syy splt_yy (T or F, ask for seperate grid output if source grid is Yin Yang)'
       PRINT *, '   -zland zland (T or F, ask to zero field values where mask is 0 (land) )'
       PRINT *, '   -nbdyf nbdyf (Number of BDY points to extrapolate all around, ...'
       PRINT *, '                 e.g. nbdyt=2 ==> data(3) is copied to data(2) and data(1) )'
       PRINT *, '   -dxred dxred (number of points to jump in build the KDtree in x (default = 1)'
       PRINT *, '   -dyred dyred (number of points to jump in build the KDtree in y (default = 1)'
       PRINT *, '   -maskd maskd (T or F for using the reference mask as the minimum mask of destination)'
       PRINT *, '   -snbdy snbdy (optional, N-pole (degenerate grid point),'
       PRINT *, '                 T-fold or F-fold at the northern boundary of the source grid'
       PRINT *, '                 default is void. if T or F, then snrep=1 by default'
       PRINT *, '   -snrep snrep (if snbdy=T,F, defines the number of repeated rows on the grid,'
       PRINT *, '                 e.g. 1 for NEMO and 0 for CICE)'
       PRINT *, '   -sebdy sebdy (optional, T or F, True means East-West periodic condition on the source grid)'
       PRINT *, '                 default is F. if T, then serep=2 by default'
       PRINT *, '   -serep serep (if sebdy=T, defines the number of repeated rows on the grid,'
       PRINT *, '                 e.g. 1 for GEM, 2 for NEMO and 0 for CICE)'
       PRINT *, '   -dnbdy dnbdy (optional, N-pole (degenerate grid point),'
       PRINT *, '                 T-fold or F-fold at the northern boundary of the destination grid'
       PRINT *, '                 default is void. if T or F, then snrep=1 by default'
       PRINT *, '   -dnrep dnrep (if dnbdy=T,F, defines the number of repeated rows on the grid,'
       PRINT *, '                 e.g. 1 for NEMO and 0 for CICE)'
       PRINT *, '   -debdy debdy (optional, T or F, True means East-West periodic condition on the destination grid)'
       PRINT *, '                 default is F. if T, then serep=2 by default'
       PRINT *, '   -derep derep (if debdy=T, defines the number of repeated rows on the grid,'
       PRINT *, '                 e.g. 1 for GEM, 2 for NEMO and 0 for CICE)'
       PRINT *, '   -lval landval(land value, default is -999.9 (do nothing, i.e. whatever you get out'
       PRINT *, '                 of the interpolator, plus zero over uncovered regions, i.e. the original behavior),'
       PRINT *, '                 999.9 is forced all masked points to the mean of the all active points,'
       PRINT *, '                 any other value is considered the user choice for all masked points,'
       PRINT *, '                 which can compromise the rmnlib compression if not chosen carefully'
       PRINT *, '                 also interarcts with the zland value (zland=T => lval=0, lval/=-999.9 => zland=T,'
       PRINT *, '                 lval=-999.9 => zland=F)'
       CALL my_abort
    ENDIF

    end subroutine print_usage

    subroutine check_arguments
    implicit none
    ! locals
    character(len=1000) :: carg
    integer jarg

    jarg=1
    DO WHILE (jarg < narg )
      CALL getarg (jarg, carg)
      SELECT CASE ( carg)
      CASE ('-fs' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fs=TRIM(carg)
      CASE ('-fstype' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fstype=TRIM(carg)
      CASE ('-fr' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fr=TRIM(carg)
      CASE ('-nr' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; nomvar_ref=TRIM(carg)
      CASE ('-fd' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fd=TRIM(carg)
      CASE ('-fdtype' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fdtype=TRIM(carg)
      CASE ('-ns' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; nomvar_src=TRIM(carg)
      CASE ('-mtrsh' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) mtrsh
      CASE ('-nsprd' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) nsprd
      CASE ('-sdfrm' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           sdfrm=.false.
          elseif (TRIM(carg).eq.'T') then
           sdfrm=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -sdfrm T or -sdfrm F'
           CALL my_abort
         endif
      CASE ('-owgts' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           owgts=.false.
          elseif (TRIM(carg).eq.'T') then
           owgts=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -owgts T or -owgts F'
           CALL my_abort
         endif
      CASE ('-intyp' )
         jarg=jarg+1 ; CALL getarg(jarg,intyp) 
         select case(intyp)
         case('mixt','aggr','bilin','nearst','bicub','mixtbc','scrip0')
         case default
           PRINT *,' Unknown option -intyp :', TRIM(intyp)
           CALL my_abort
         end select
      CASE ('-nwgt' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) nadis
      CASE ('-naggrmax' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) naggrmax
      CASE ('-status' )
         jarg=jarg+1 ; CALL getarg(jarg,file_status)
      CASE ('-npak' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) npak_ou
      CASE ('-glbsrc' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           glbsrc=.false.
          elseif (TRIM(carg).eq.'T') then
           glbsrc=.true.
           snbdy = 'T'
           snrep = 1
           sebdy = .true.
           serep = 2
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -glbsrc T or -glbsrc F'
           CALL my_abort
         endif
      CASE ('-uwgdir' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; uwgdir=TRIM(carg)
      CASE ('-myvtab' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; myvtab=TRIM(carg)
      CASE ('-syy' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           splt_yy=.false.
          elseif (TRIM(carg).eq.'T') then
           splt_yy=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -syy T or -syy F'
           CALL my_abort
         endif
      CASE ('-zland' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           zland=.false.
          elseif (TRIM(carg).eq.'T') then
           zland=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -zland T or -zland F'
           CALL my_abort
         endif
      CASE ('-maskd' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           maskd=.false.
          elseif (TRIM(carg).eq.'T') then
           maskd=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -maskd T or -maskd F'
           CALL my_abort
         endif
      CASE ('-nbdyf' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) nbdyf
      CASE ('-dxred' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) dxred
      CASE ('-dyred' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) dyred
      CASE ('-snbdy')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) snbdy
         SELECT CASE( snbdy )
         CASE('N')
         CASE('F','T')
           snrep = 1
         CASE DEFAULT
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Accepted values: snbdy=N, T or F'
           CALL my_abort
         END SELECT
      CASE ('-snrep')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) snrep
      CASE ('-sebdy')
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           sebdy=.false.
          elseif (TRIM(carg).eq.'T') then
           sebdy=.true.
           serep = 2
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -sebdy T or -sebdy F'
           CALL my_abort
         endif
      CASE ('-serep')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) serep
      CASE ('-dnbdy')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) dnbdy
         SELECT CASE( dnbdy )
         CASE('N')
         CASE('F','T')
           dnrep = 1
         CASE DEFAULT
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Accepted values: dnbdy=N, T or F'
           CALL my_abort
         END SELECT
      CASE ('-dnrep')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) dnrep
      CASE ('-debdy')
         jarg=jarg+1 ; CALL getarg(jarg,carg) 
         if (TRIM(carg).eq.'F') then
           debdy=.false.
          elseif (TRIM(carg).eq.'T') then
           debdy=.true.
           derep = 2
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -dendy T or -dendy F'
           CALL my_abort
         endif
      CASE ('-derep')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) derep
      CASE ('-lval')
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) landval
      CASE DEFAULT
         PRINT *,' Unknown option :', TRIM(carg)
         CALL my_abort
      END SELECT
      jarg=jarg+1
    ENDDO

    iowgt=.false.
    if (uwgdir.ne.'NONE') iowgt=.true.
    if (owgts.and.iowgt) then
       write(*,*) 'CSTINTRP: INCOMPATIBLE OPTION owgts=T, and uwgdir'
       call my_abort
    endif

    if (naggrmax.lt.4) then
       write(*,*) 'naggrmax must be at least equal to 4'
       call my_abort
    endif

    ! process land value
    if (zland) landval = 0.0
    if (landval==-999.9) then
        zland = .false. ! default behavior of doing nothing specific
    else
        zland = .true.  ! enforce a specific value over land given by landval
    endif


! 6 sounds a bit far-fetch, so I will leave the user to decide
!    if (bicub.or.mixtbc) nsprd = max( nsprd, 6 )

! deprecated option, only necessary if the calling script does not stop when encountering an error
    if ( file_status /= '' ) then
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstintrp_status=ABORT'
      close(99)
    endif

    end subroutine check_arguments

    subroutine read_vector_table
    implicit none
    ! locals
    integer ivec, unitv

! Definition de la table des vecteurs

    nvec=23
    allocate(vecttbl(nvec,2))
    vecttbl(1, 1)='UU'    ; vecttbl(1, 2)='VV'
    vecttbl(2, 1)='UWIN'  ; vecttbl(2, 2)='VWIN'
    vecttbl(3, 1)='UUOR'  ; vecttbl(3, 2)='VVOR'
    vecttbl(4, 1)='UICE'  ; vecttbl(4, 2)='VICE'
    vecttbl(5, 1)='ISTX'  ; vecttbl(5, 2)='ISTY'
    vecttbl(6, 1)='TAUX'  ; vecttbl(6, 2)='TAUY'
    vecttbl(7, 1)='UU2'   ; vecttbl(7, 2)='VV2'
    vecttbl(8, 1)='UWAT'  ; vecttbl(8, 2)='VWAT'
    vecttbl(9, 1)='UWA2'  ; vecttbl(9, 2)='VWA2'
    vecttbl(10,1)='UUW'   ; vecttbl(10,2)='VVW'
    vecttbl(11,1)='UU2W'  ; vecttbl(11,2)='VV2W'
    vecttbl(12,1)='UTAW'  ; vecttbl(12,2)='VTAW'
    vecttbl(13,1)='UUI'   ; vecttbl(13,2)='VVI'
    vecttbl(14,1)='UNOC'  ; vecttbl(14,2)='VNOC'
    vecttbl(15,1)='UBOC'  ; vecttbl(15,2)='VBOC'
    vecttbl(16,1)='SIXI'  ; vecttbl(16,2)='SIYI'
    vecttbl(17,1)='SOXI'  ; vecttbl(17,2)='SOYI'
    vecttbl(18,1)='SAXI'  ; vecttbl(18,2)='SAYI'
    vecttbl(19,1)='TBU'   ; vecttbl(19,2)='TBV'
    vecttbl(20,1)='STLX'  ; vecttbl(20,2)='STLY'
    vecttbl(21,1)='SCXI'  ; vecttbl(21,2)='SCYI'
    vecttbl(22,1)='UGEO'  ; vecttbl(22,2)='VGEO'
    vecttbl(23,1)='UD'    ; vecttbl(23,2)='VD'

    if ( trim(myvtab) /= '' ) then
      deallocate(vecttbl)
      unitv=get_unit()
      open(unitv,file=trim(myvtab),status='old',form='formatted')
      read(unitv) nvec
      allocate(vecttbl(nvec,2))
      do ivec=1,nvec
        read(unitv,'(a4,x,a4)') vecttbl(ivec, 1), vecttbl(ivec, 2)
      enddo
      close(unitv)
    endif

    end subroutine read_vector_table

    subroutine open_src_ref_dst_files
    implicit none
    ! locals
    integer ierr, status_st90

! Ouverture des fichiers

    ierr = stdopen(trim(fd),'new',unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fd=',trim(fd) 
      CALL my_abort
    endif

    ierr = stdopen(trim(fr),'old',unitr,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fr=',trim(fr) 
      CALL my_abort
    endif

    ierr = stdopen(trim(fs),'old',units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fs=',trim(fs) 
      CALL my_abort
    endif

    end subroutine open_src_ref_dst_files

    subroutine read_list_var_grid_src
    implicit none
    ! locals
    integer ierr, status_st90, igs
    type(std_grid), allocatable, dimension(:) :: grdp

! Classement des grilles et des variables du fichier source

! ... toutes les variables
    nvs = getstdnvar(units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( nvs <= 0 ) then 
      write(*,*) 'CSTINTRP getstdnvar => problem reading source file variables'
      write(*,*) 'CSTINTRP getstdnvar => should have at least one variable to read'
      write(*,*) 'CSTINTRP getstdnvar => nvs,status_st90=',nvs,status_st90
      CALL my_abort
    endif
    allocate(vars(nvs), vars_g(nvs), vars_gs(nvs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP vars, vars_g, vars_gs allocation problem'
      write(*,*) 'CSTINTRP nvs=',nvs
      CALL my_abort
    endif
    ierr = fillstdvar(vars,nvs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP fillstdvar => problem filling source variable informations'
      write(*,*) 'CSTINTRP fillstdvar => ierr,status_st90=',ierr,status_st90
      CALL my_abort
    endif
! ... toutes les grilles
    write(*,*) 'CSTINTRP: SCANNING SOURCE GRIDS'
    ngs = getstdngrid(units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ngs <= 0 ) then 
      write(*,*) 'CSTINTRP getstdngrid => problem reading source file grids'
      write(*,*) 'CSTINTRP getstdngrid => should have at least one grid to read'
      write(*,*) 'CSTINTRP getstdngrid => ngs,status_st90=',ngs,status_st90
      CALL my_abort
    endif
    allocate(grds(ngs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP grds allocation problem'
      write(*,*) 'CSTINTRP ngs=',ngs
      CALL my_abort
    endif
    allocate(grdp(ngs))
    ierr = fillstdgrid(grdp,ngs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    do igs = 1, ngs
      grds(igs) % std_grd = grdp(igs)
      grds(igs) % dxtree  = dxred
      grds(igs) % dytree  = dyred
      grds(igs) % nbdy    = snbdy
      grds(igs) % nrepeat = snrep
      grds(igs) % webdy   = sebdy
      grds(igs) % werepeat= serep
      call init_grid_parameter(grds(igs),units)
      call init_sub_grid(grds(igs))
      call getll(grds(igs),units)
    enddo
    deallocate(grdp)
    if ( ierr < 0 ) then
      write(*,*) 'CSTINTRP fillstdgrid => problem filling source file grid informations'
      write(*,*) 'CSTINTRP fillstdgrid => ierr,status_st90=',ierr,status_st90
      CALL my_abort
    endif

    end subroutine read_list_var_grid_src

    subroutine read_list_var_grid_ref
    implicit none
    ! locals
    integer ierr, status_st90, igr
    type(std_grid), allocatable, dimension(:) :: grdp

! Classement des grilles du fichier reference

! ... toutes les grilles (de nomvar_ref)
    write(*,*) 'CSTINTRP: SCANNING DESTINATION GRIDS'
    ngr = getstdngrid(unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    if ( ngr <= 0 ) then 
      write(*,*) 'CSTINTRP getstdngrid => problem reading reference file grids'
      write(*,*) 'CSTINTRP getstdngrid => should have at least one grid to read'
      write(*,*) 'CSTINTRP getstdngrid => ngr,status_st90=',ngr,status_st90
      CALL my_abort
    endif
    allocate(grdr(ngr), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTINTRP grdr allocation problem'
      write(*,*) 'CSTINTRP ngr=',ngr
      CALL my_abort
    endif
    allocate(grdp(ngr))
    ierr = fillstdgrid(grdp,ngr,unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    do igr = 1, ngr
      grdr(igr) % std_grd = grdp(igr)
      grdr(igr) % dxtree  = dxred
      grdr(igr) % dytree  = dyred
      grdr(igr) % nbdy    = dnbdy
      grdr(igr) % nrepeat = dnrep
      grdr(igr) % webdy   = debdy
      grdr(igr) % werepeat= derep
      call init_grid_parameter(grdr(igr),unitr)
      call init_sub_grid(grdr(igr))
      call getll(grdr(igr),unitr)
      call init_grid_mask(grdr(igr), nomvar_ref, unitr, maskd)
    enddo
    deallocate(grdp)
    if ( ierr < 0 ) then
      write(*,*) 'CSTINTRP fillstdgrid => problem filling source file grid informations'
      write(*,*) 'CSTINTRP fillstdgrid => ierr,status_st90=',ierr,status_st90
      CALL my_abort
    endif

    end subroutine read_list_var_grid_ref

    subroutine treat_each_grid_src_dst(F_grds, F_grdr, nvartrt)
    !------------------------------------------------------------------------------
    ! main call for treating interpolation between src and dst grids`
    !------------------------------------------------------------------------------
    ! arguments
    type(type_grid), target &
                        :: F_grds, F_grdr
    integer             :: nvartrt
    ! locals
    integer             :: icnt
    integer ivs_g
    type(std_grid), pointer :: grdp
    logical             :: skip, mask_change, sortip1
    integer, dimension(:,:), allocatable &
                        :: mask_n, mask_o
    integer ierr, ivec, status_st90, id_m, i, j
    integer nx, ny
    logical vector
    character(len=4)  :: nomvar_v
    type(variable), allocatable, dimension(:) :: var_grds
    type(std_variable) :: varmask

    !---------------------------------------------------------------------------------------

    ! Initialize parameters and validate options with source and destination grid
    call intrp_itfc_params( F_grds, F_grdr )

    ! Allocate memory and copy grid descriptors to destination file
    call intrp_itfc_init( F_grds, F_grdr, nomvar_ref)

    !---------------------------------------------------------------------------------------
    ! find the variables associated with that particular source grid
    grdp => F_grds % std_grd
    nx   =  grdp % ni
    ny   =  grdp % nj
    ierr = filtergstdvar(vars_g,nvs_g,vars,nvs, &
                             grdp % grtyp,                   &
                             grdp % ig1,                     &
                             grdp % ig2,                     &
                             grdp % ig3,                     &
                             grdp % ig4,                     &
                             grdp % ni,                      &
                             grdp % nj)
    write(*,*) 'CSTINTRP: N. fields =', nvs_g
    write(*,*) 'CSTINTRP: ------------------------------------'
    if ( ierr < 0 ) then
      write(*,*) 'CSTINTRP filtergstdvar => problem filtering grid information'
      write(*,*) 'CSTINTRP filtergstdvar => grds=',grdp
      CALL my_abort
    endif

    ! reduce the number of variables by removing the mask and vector v-components
    ! sort them by ip1 so that we treat several ip1 at once
    ! (should share the same mask, therefore the same weights)
    ! introduce new type var
    allocate(var_grds(nvs_g))
    call reduce_var_list(nvs_g, vars_g, var_grds)
    do ivs_g=1,nvs_g
      write(*,*) 'CSTINTRP: remaining var',ivs_g, &
            var_grds(ivs_g) % varstd % nomvar, &
            var_grds(ivs_g) % varstd % typvar, &
            var_grds(ivs_g) % varstd % ip1, &
            var_grds(ivs_g) % varstd % npas, &
            var_grds(ivs_g) % ln_mask, & 
            var_grds(ivs_g) % ln_vector
    enddo

    ! If masks are present,
    ! treat one ip1 at a time to keep interpolation weights in memory
    ! as often as possible (avoid change of masks with depth)

    ! only if a change of mask status occurs
    mask_change=.true. !new grid
    icnt=0
    do ivs_g = 1, nvs_g
      icnt=icnt+1
      nvartrt = nvartrt + 1
      id_m=0
      if ( var_grds(ivs_g) % ln_mask) then
        ! read mask
        varmask = var_grds(ivs_g) % mask
        allocate(mask_n(nx,ny))
        write(*,*) 'read mask for variable ',varmask % nomvar, varmask % ip1, varmask % npas
        ierr = getstdvar(mask_n,nx,ny,units,varmask,ierr_st90=status_st90)
        if ( ierr < 0 ) then
          write(*,*) 'CSTINTRP getstdvar => problem reading mask_n'
          write(*,*) 'CSTINTRP getstdvar => vars_gs(id_m)=', &
                                            vars_gs(id_m)
          CALL my_abort
        endif
      else
        allocate(mask_n(nx,ny)); mask_n(:,:) = 1
      endif
      if (.not.allocated(mask_o)) then !first iteration
         allocate(mask_o(nx,ny))
         mask_o(:,:) = 0
      endif
      mask_change=.false.
      do j=1,ny
          do i=1,nx
             if (mask_n(i,j) /= mask_o(i,j)) then; mask_change = .true.; exit; endif
          enddo 
      enddo 
      mask_o(:,:) = mask_n(:,:)
      write(*,*) 'CSTINTRP HAS MASK CHANGED?: mask_change=',mask_change
      call intrp_itfc_main( var_grds(ivs_g), F_grds,F_grdr, mask_change )
      if (allocated(mask_n)) deallocate(mask_n)
    enddo ! loop over variables

    if (allocated(mask_o)) deallocate(mask_o)
    deallocate(var_grds)

    end subroutine treat_each_grid_src_dst


    subroutine close_all_files
    implicit none
    ! locals
    integer ierr, status_st90

! Fermeture des fichiers
    ierr = stdclose(units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fs=',trim(fs)
      CALL my_abort
    endif
    ierr = stdclose(unitr,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fr=',trim(fr)
      CALL my_abort
    endif
    ierr = stdclose(unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fd=',trim(fd)
      CALL my_abort
    endif

    end subroutine close_all_files

    subroutine reduce_var_list(nv, varl, varg)
    !---------------------------------------------------------------
    ! sort the variables by ip1 and remove duplicates
    ! such as mask and v-components
    !---------------------------------------------------------------
    implicit none

    ! arguments
    integer nv
    type(std_variable), intent(in)  :: varl(nv)
    type(variable), intent(out)     :: varg(nv)

    ! locals
    type(std_variable) :: var2(nv)
    logical vec(nv), lm(nv)
    character(len=4) :: nomv(nv)
    integer iv, ivec, ivnew, ierr, iv2, nv_old
    logical vector, skip, sortip1, mask
    character(len=4) :: nomvar_v
    type(std_variable) :: var_test
    integer list_glob(nv), id_var

    var2(1:nv) = varl( 1:nv)
    sortip1=.false.
    ! only ip1-sort masked variables
    do iv=1,nv
      if (varl(iv) % typvar(2:2) == '@') then; sortip1 = .true.; exit; endif
    enddo
    if (sortip1) then
      ierr = sortbyip1stdvar(var2, varl, nv)
      if ( ierr < 0 ) then
        write(*,*) 'CSTINTRP sortbyip1stdvar => problem sorting variable by ip1'
        write(*,*) 'CSTINTRP sortbyip1stdvar => vars_g(1:nvs_g)=', varl(1:nv)
        CALL my_abort
      endif
    endif
    ! now remove duplicate of mask or vector v-component
    nv_old = nv
    ivnew = 0
    do iv = 1, nv
      skip      = .false.
      vector    = .false.
      nomvar_v  = 'NONE'
      mask      = .false.
      if ( var2(iv) % typvar(2:2) == '@' ) mask = .true.
      if ( var2(iv) % typvar == '@@' )     skip = .true.
      if ( .not. skip ) then
        do ivec=1,nvec
          if ( var2(iv)%nomvar == trim(vecttbl(ivec,2)) ) skip = .true.
          if ( var2(iv)%nomvar == trim(vecttbl(ivec,1)) ) then
            vector = .true.
            nomvar_v = vecttbl(ivec,2)
          endif
        enddo
      endif
      if (.not.skip) then
         ivnew = ivnew + 1
         varg(ivnew) % varstd = var2(iv)
         vec (ivnew) = vector
         nomv(ivnew) = nomvar_v
         lm  (ivnew) = mask
         list_glob(ivnew) = iv
      endif
    enddo
    nv = ivnew

    ! retrieve records from mask or vector
    do iv = 1, nv
      if (lm(iv)) then
        id_var = list_glob(iv)
        iv2 = getstdmaskid(id_var,var2, nv_old)
        if ( iv2 > 0 ) then
          varg(iv) % ln_mask = .true.
          varg(iv) % mask = var2(iv2)
        endif
      endif
      if (vec(iv)) then
        var_test = varg(iv) % varstd
        var_test % nomvar = nomv(iv)
        do iv2 = 1, nv_old
          if( equalstdvar(var_test, var2(iv2), .false.) ) then
             varg(iv) % ln_vector = .true.
             varg(iv) % varstdv   = var2(iv2)
             exit
          endif
        enddo
      endif
    enddo

    end subroutine reduce_var_list

end module cstintrp_mod
