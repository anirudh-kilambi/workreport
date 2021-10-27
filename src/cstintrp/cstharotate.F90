    program cstharotate

!! Compute the angle of a structured grid (on each point)
!! relative to the east-north axis (regular lat-lon grid)
!! Destination angle ==> LAAN

!! Use LAAN to rotate elliptical currents or 
!! transports (in harmonic form ==> amplitude, phase (degrees))
!! (if any in the form UAMP,VAMP,UPHA,VPHA in the source file)
!! 
!! Just copy HAMP and HPHA to the destination file if any 
!! (tidal elevation, amplitude and phase)

!! Author: Francois Roy (2016)

!! Method:

!! Use convivial interface
!! Use cstintrp sub-programs
!! Write to std file

!! Arguments:
!! 

!! Modules:
    USE std
    USE utils
    USE cdf
    USE cdf_std
    USE llxy

    IMPLICIT NONE

!! Argument variables:

    character(len=1000) :: fs,fd         !Source,destination files,respectively
    character(len=4)   :: nomvar_src    !Nomvar selection in source file (default nomvar_src=' ')
    integer            :: iargc,jarg,narg, units, unitd
    character(len=1000) :: carg,file_status

!! Local variables:

    logical status , maskou, chkmsk

    !... variable management
    type (std_variable),  &
       dimension(:),      &
       allocatable      :: vars,    &  !All variables from source file
                           varsg       !vars for one particular horizontal grid
    type (std_variable) :: varwg
    integer             :: nvs, nvsg, ivsg

    type (std_grid),      &
       dimension(:),      &
       allocatable      :: hgrds    !All grids from source file
    integer             :: ngs, igs
    integer             :: ierr, status_st90, i, j, it, &
                           nxs, nys, nxd, nyd, nxw, nyw, gdidw, ns
    real(kind=8)        :: wdlat,wdlon,xnxw,xnyw,pi,torad
    real(kind=8)        :: A1u,A2u,A1v,A2v,theta1,theta2,r1,r2,Ku1,Ku2,Kv1,Kv2
    real(kind=4), parameter  :: epsilon = 1.e-5
    real(kind=4)        :: latmin,lonmin,latmax,lonmax
    integer,              &
     dimension(:,:),      &
     allocatable        :: masks
    real(kind=4),         &
     dimension(:,:),      &
     allocatable        :: lats,lons,latw,lonw,xs,ys,angles, &
                           uamp,upha,vamp,vpha
    character(len=14)   :: gdefw,gdefs

! Definition des arguments par defaut

    nomvar_src=' '
    file_status=''

! Traitement des arguments
    narg= iargc()
    IF ( narg < 4 .or. mod(narg,2).ne.0 ) THEN
       PRINT *, ' Usage : cstharotate '
       PRINT *, ' -fs filesrc (source file to compute angle)'
       PRINT *, ' -fd filedst (destination file, with angle (LAAN), ...'
       PRINT *, '              and if so in the source file the copied tidal ...'
       PRINT *, '              elevation harmonics (HAMP, HPHA), ...'
       PRINT *, '              and the rotated tidal ellipses (UAMP, UPHA, VAMP, VPHA))'
       PRINT *, ' Optional arguments:'
       PRINT *, '   -ns nomvar_src (nomvar selection in src file,default all variables)'
       PRINT *, '   -status file_status (output program status in file_status for op needs)'
       STOP
    ENDIF
    jarg=1
    DO WHILE (jarg < narg )
      CALL getarg (jarg, carg)
      SELECT CASE ( carg)
      CASE ('-fs' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fs=TRIM(carg)
      CASE ('-fd' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fd=TRIM(carg)
      CASE ('-ns' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; nomvar_src=TRIM(carg)
      CASE ('-status' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; file_status=TRIM(carg)
         status=.true.
      CASE DEFAULT
         PRINT *,' Unknown option :', TRIM(carg) ; STOP
      END SELECT
      jarg=jarg+1
    ENDDO

    if (status) then
      status=.true.
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstharotate_status=ABORT'
      close(99)
    endif

! Ouverture des fichiers

    ierr = stdopen(trim(fd),'new',unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fd=',trim(fd) 
      STOP 
    endif

    ierr = stdopen(trim(fs),'old',units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fs=',trim(fs) 
      STOP 
    endif

! Classement des grilles et des variables du fichier source

! ... toutes les variables source
    nvs = getstdnvar(units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( nvs <= 0 ) then 
      write(*,*) 'CSTHAROTATE getstdnvar => problem reading source file variables'
      write(*,*) 'CSTHAROTATE getstdnvar => should have at least one variable to read'
      write(*,*) 'CSTHAROTATE getstdnvar => nvs,status_st90=',nvs,status_st90
      STOP
    endif
    allocate(vars(nvs), varsg(nvs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTHAROTATE vars, varsg, ... allocation problem'
      write(*,*) 'CSTHAROTATE nvs=',nvs
      STOP
    endif
    ierr = fillstdvar(vars,nvs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTHAROTATE fillstdvar => problem filling source variable informations'
      write(*,*) 'CSTHAROTATE fillstdvar => ierr,status_st90=',ierr,status_st90
      STOP
    endif
! ... toutes les grilles horizontales
    write(*,*) 'CSTHAROTATE: SCANNING SOURCE HORIZONTAL GRIDS'
    ngs = getstdngrid(units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ngs <= 0 ) then 
      write(*,*) 'CSTHAROTATE getstdngrid => problem reading source file grids'
      write(*,*) 'CSTHAROTATE getstdngrid => should have at least one grid to read'
      write(*,*) 'CSTHAROTATE getstdngrid => ngs,status_st90=',ngs,status_st90
      STOP
    endif
    allocate(hgrds(ngs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTHAROTATE hgrds allocation problem'
      write(*,*) 'CSTHAROTATE ngs=',ngs
      STOP
    endif
    ierr = fillstdgrid(hgrds,ngs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ierr < 0 ) then
      write(*,*) 'CSTHAROTATE fillstdgrid => problem filling source file grid informations'
      write(*,*) 'CSTHAROTATE fillstdgrid => ierr,status_st90=',ierr,status_st90
      STOP
    endif

! Constantes
    pi=acos(-1.)
    torad=pi/180.

! Boucle principale sur les grilles source

    do igs = 1, ngs
      ierr = filtergstdvar(varsg,nvsg,vars,nvs, &
                           hgrds(igs)%grtyp,    &
                           hgrds(igs)%ig1,      &
                           hgrds(igs)%ig2,      &
                           hgrds(igs)%ig3,      &
                           hgrds(igs)%ig4,      &
                           hgrds(igs)%ni,       &
                           hgrds(igs)%nj)
      if ( ierr < 0 ) then
        write(*,*) 'CSTHAROTATE filtergstdvar => problem filtering grid information'
        write(*,*) 'CSTHAROTATE filtergstdvar => igs=',igs
        write(*,*) 'CSTHAROTATE filtergstdvar => hgrds(igs)=',hgrds(igs)
        STOP
      endif
      ierr = xferstdparpos(varsg,nvsg,units,unitd)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTHAROTATE xferstdparpos => problem transfering par pos parameters'
        write(*,*) 'CSTHAROTATE hgrds(igs) =', hgrds(igs)
        STOP
      endif
      ! Basic allocation
      nxs=hgrds(igs)%ni
      nys=hgrds(igs)%nj
      nxd=nxs
      nyd=nys

      write(*,*) 'CSTHAROTATE treating grid with dimensions nxs,nys=',nxs,nys
      write(*,*) 'CSTHAROTATE treating grid with dimensions hgrds(igs)=',hgrds(igs)

      allocate(lats(nxs,nys),lons(nxs,nys),masks(nxs,nys))
      allocate(xs(nxs,nys),ys(nxs,nys),angles(nxs,nys))

      ierr = getstdll(lats,lons,hgrds(igs)%ni,hgrds(igs)%nj, &
                      hgrds(igs)%grtyp, &
     &                hgrds(igs)%ig1,hgrds(igs)%ig2, &
                      hgrds(igs)%ig3,hgrds(igs)%ig4,units,ierr_st90=status_st90)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTHAROTATE getstdll => WARNING problem getting source lat-lon fields'
        write(*,*) 'CSTHAROTATE getstdll => WARNING ierr,status_st90=',ierr,status_st90
        STOP
      endif
      masks(:,:)=1 

      ! Build regular lat-lon work grid to compute angles based on source grid

      nyw=max(nxs,nys)*2
      nxw=nyw/2
      allocate(latw(-5:nxw+5,-5:nyw+5),lonw(-5:nxw+5,-5:nyw+5))

      latmin=90.
      lonmin=9999999.
      latmax=-90.
      lonmax=-9999999.
      do j=1,nys
      do i=1,nxs
        latmin=min(lats(i,j),latmin)
        lonmin=min(lons(i,j),lonmin)
        latmax=max(lats(i,j),latmax)
        lonmax=max(lons(i,j),lonmax)
      enddo
      enddo
      xnxw=nxw
      xnyw=nyw
      wdlat=(latmax-latmin)/(xnyw-1)
      wdlon=(lonmax-lonmin)/(xnxw-1)
      do j=-5,nyw+5
      do i=-5,nxw+5
        latw(i,j)=latmin+(j-1)*wdlat
        lonw(i,j)=lonmin+(i-1)*wdlon
      enddo
      enddo

      ! Compute angles from source grid relastive to work grid (regular lat-lon)

        !! Source grid here correspond to destination grid in ll2xy, careful!!

        !! Compute X,Y components and angle A from grid destination (lons,lats here)
        !! relative to source grid (lonw,latw here)
        !! A is relative to X direction (i==>nxs)
        !! and positive anticlockwise (radians)
        !! Return values of -9. for masked points, M=0

      ! ...

      gdefw='regular_latlon'
      gdefs='bidon'
      gdidw=-1  !bidon

      maskou=.false.
      call ll2xy(lats,lons,xs,ys,angles,masks,nxs,nys,gdefs, &
   &             latw,lonw,nxw+11,nyw+11,gdefw,gdidw,maskou)

      do j=1,nys
      do i=1,nxs
        if (masks(i,j) /= 1) then
          write(*,*) 'MASKED POINTS AT i,j=',i,j
          write(*,*) 'lats(i,j),lons(i,j)=',lats(i,j),lons(i,j)
          write(*,*) 'xs(i,j),ys(i,j)=',xs(i,j),ys(i,j)
        endif
      enddo
      enddo

      if (maskou) then
        write(*,*) 'CSTHAROTATE: ll2xy => problem maksou=T'
        STOP
      endif

      varwg=varsg(1)
      varwg%nomvar='LAAN'
      varwg%typvar='P'
      varwg%ip1=0
      varwg%nbits=32
      varwg%datyp=133
      varwg%xtype=2
      varwg%cdfscale=1.
      varwg%cdfoffset=0.

      ierr = putstdvar(unitd, angles(:,:),  &
                       nxd, nyd, varwg, ierr_st90=status_st90)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTHAROTATE putstdvar problem => angles'
        write(*,*) 'CSTHAROTATE status_st90=',status_st90
        write(*,*) 'CSTHAROTATE varwg=',varwg
        STOP
      endif

      !For rotation of vectors use the inverse of grid angle
      angles(:,:)=-angles(:,:)

      ! Grab and rotate tidal harmonics if so and needed
      do ivsg = 1, nvsg
        ! Just copy tidal elevation harmonics
        if ( varsg(ivsg)%nomvar == 'HAMP' .or. varsg(ivsg)%nomvar == 'LAMP' ) then
     
          if ( varsg(ivsg)%typvar /= '@@' ) then
            if (.not. allocated(uamp)) allocate(uamp(nxs,nys))
            if (.not. allocated(upha)) allocate(upha(nxs,nys))
            ! Note here uamp,upha used for HAMP and HPHA
            ierr = getstdvar(uamp(:,:),nxs,nys,units, &
                             varsg(ivsg),ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (HAMP) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
              STOP
            endif
            varwg = varsg(ivsg)
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'a'
            ierr = putstdvar(unitd, uamp(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (HAMP)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg = varsg(ivsg)
            if ( varsg(ivsg)%nomvar == 'LAMP' ) then
              varwg%nomvar='LPHA'
             else
              varwg%nomvar='HPHA'
            endif
            ierr = getstdvar(upha(:,:),nxs,nys,units, &
                             varwg,ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (HPHA) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
              STOP
            endif
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'p'
            ierr = putstdvar(unitd, upha(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (HPHA)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            ! Here use xs,ys as z1, z2
            do j=1,nys
            do i=1,nxs
              xs(i,j) =  uamp(i,j)*cos(upha(i,j)*torad)
              ys(i,j) = -uamp(i,j)*sin(upha(i,j)*torad) !OSU convention
            enddo
            enddo
            varwg%nomvar='Z1'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z1'
            ierr = putstdvar(unitd, xs(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z2'
            ierr = putstdvar(unitd, ys(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
           else
            ierr = getstdvar(masks(:,:),nxs,nys,units, &
                             varsg(ivsg),ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (HAMP MASKS) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
              STOP
            endif
            varwg = varsg(ivsg)
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'a'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (HAMP MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            if ( varsg(ivsg)%nomvar == 'LAMP' ) then
              varwg%nomvar='LPHA'
             else
              varwg%nomvar='HPHA'
            endif
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'p'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (HPHA MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z1'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z1'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1 MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z2'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2 MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
          endif
        endif
        ! Rotate tidal ellipses
        if ( varsg(ivsg)%nomvar == 'UAMP' ) then
          if ( varsg(ivsg)%typvar /= '@@' ) then
            ! Assumes UPHA, VAMP, VPHA present in this case 
            if (.not. allocated(uamp)) allocate(uamp(nxs,nys))
            if (.not. allocated(upha)) allocate(upha(nxs,nys))
            if (.not. allocated(vamp)) allocate(vamp(nxs,nys))
            if (.not. allocated(vpha)) allocate(vpha(nxs,nys))
            ierr = getstdvar(uamp(:,:),nxs,nys,units, &
                             varsg(ivsg),ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (UAMP) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
              STOP
            endif
            !Work var for other components with same attributes
            varwg = varsg(ivsg)
            varwg%nomvar='UPHA'
            ierr = getstdvar(upha(:,:),nxs,nys,units, &
                             varwg,ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (UPHA) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varwg%nomvar=',varwg%nomvar
              STOP
            endif
            varwg%nomvar='VAMP'
!            varwg%etiket=' '  ! Expect V in the last character
            ns=len(trim(varsg(ivsg)%etiket))
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v'
            ierr = getstdvar(vamp(:,:),nxs,nys,units, &
                             varwg,ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (VAMP) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varwg%nomvar=',varwg%nomvar
              STOP
            endif
            varwg%nomvar='VPHA'
            ierr = getstdvar(vpha(:,:),nxs,nys,units, &
                             varwg,ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (VPHA) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varwg%nomvar=',varwg%nomvar
              STOP
            endif
            ! Main action for rotation
            ! See documentation on F. Roy's WIKI PAGE
            do j=1,nys
            do i=1,nxs
              A1u=uamp(i,j)*cos(upha(i,j)*torad)
              A2u=uamp(i,j)*sin(upha(i,j)*torad)
              A1v=vamp(i,j)*cos(vpha(i,j)*torad)
              A2v=vamp(i,j)*sin(vpha(i,j)*torad)
              theta1=atan2(A1v-A2u,A1u+A2v)
              theta2=atan2(A2u+A1v,A1u-A2v)
              r1=sqrt( 0.25*( (A1u+A2v)**2 + (A1v-A2u)**2 ) )
              r2=sqrt( 0.25*( (A1u-A2v)**2 + (A2u+A1v)**2 ) )
              Ku1=r1*cos(theta1+angles(i,j))+r2*cos(theta2+angles(i,j))
              Ku2=-r1*sin(theta1+angles(i,j))+r2*sin(theta2+angles(i,j))
              Kv1=r1*sin(theta1+angles(i,j))+r2*sin(theta2+angles(i,j))
              Kv2=r1*cos(theta1+angles(i,j))-r2*cos(theta2+angles(i,j))
              ! Radians
              upha(i,j)=atan2(Ku2,Ku1)
              vpha(i,j)=atan2(Kv2,Kv1)
              if (abs(cos(upha(i,j))) < epsilon ) then
                uamp(i,j)=Ku2/sin(upha(i,j))
               else
                uamp(i,j)=Ku1/cos(upha(i,j))
              endif
              if (abs(cos(vpha(i,j))) < epsilon ) then
                vamp(i,j)=Kv2/sin(vpha(i,j))
               else
                vamp(i,j)=Kv1/cos(vpha(i,j))
              endif
              ! Degrees
              upha(i,j)=upha(i,j)/torad
              vpha(i,j)=vpha(i,j)/torad
            enddo
            enddo
            varwg = varsg(ivsg)
            varwg%nomvar='UAMP'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'a'
            ierr = putstdvar(unitd, uamp(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (UAMP)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='UPHA'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'p'
            ierr = putstdvar(unitd, upha(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (UPHA)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varwg=',varwg
              STOP
            endif
            ! Here use xs,ys as z1, z2
            do j=1,nys
            do i=1,nxs
              xs(i,j) =  uamp(i,j)*cos(upha(i,j)*torad)
              ys(i,j) = -uamp(i,j)*sin(upha(i,j)*torad)  ! OSU convention
            enddo
            enddo
            varwg%nomvar='Z1U'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z1'
            ierr = putstdvar(unitd, xs(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1 u)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2U'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z2'
            ierr = putstdvar(unitd, ys(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2 u)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='VAMP'
            ns=len(trim(varsg(ivsg)%etiket))
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'a'
            ierr = putstdvar(unitd, vamp(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (VAMP)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varwg=',varwg
              STOP
            endif
            varwg%nomvar='VPHA'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'p'
            ierr = putstdvar(unitd, vpha(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (VPHA)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varwg=',varwg
              STOP
            endif
            ! Here use xs,ys as z1, z2
            do j=1,nys
            do i=1,nxs
              xs(i,j) =  vamp(i,j)*cos(vpha(i,j)*torad)
              ys(i,j) = -vamp(i,j)*sin(vpha(i,j)*torad) ! OSU convention
            enddo
            enddo
            varwg%nomvar='Z1V'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'z1'
            ierr = putstdvar(unitd, xs(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1 v)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2V'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'z2'
            ierr = putstdvar(unitd, ys(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2 v)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
           else
            ierr = getstdvar(masks(:,:),nxs,nys,units, &
                             varsg(ivsg),ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE getstdvar (UAMP MASKS) => ierr,status_st90=',ierr,status_st90
              write(*,*) 'CSTHAROTATE getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
              STOP
            endif
            varwg=varsg(ivsg)
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'a'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (UAMP MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='UPHA'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'p'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (UPHA MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z1U'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z1'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1U MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2U'
            varwg%etiket=trim(varsg(ivsg)%etiket)//'_'//'z2'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2U MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='VAMP'
            ns=len(trim(varsg(ivsg)%etiket))
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'a'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (VAMP MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='VPHA'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'p'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (VPHA MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z1V'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'z1'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z1V MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
            varwg%nomvar='Z2V'
            varwg%etiket=trim(varsg(ivsg)%etiket(1:ns-1))//'v_'//'z2'
            ierr = putstdvar(unitd, masks(:,:),  &
                             nxd, nyd, varwg, ierr_st90=status_st90)
            if ( ierr /= 0 ) then
              write(*,*) 'CSTHAROTATE putstdvar problem (Z2V MASKS)'
              write(*,*) 'CSTHAROTATE status_st90=',status_st90
              write(*,*) 'CSTHAROTATE varsg(ivsg)=',varsg(ivsg)
              STOP
            endif
          endif
        endif
      enddo

      deallocate(lats,lons,masks,xs,ys,angles,latw,lonw)
      if (allocated(uamp)) deallocate(uamp)
      if (allocated(vamp)) deallocate(vamp)
      if (allocated(upha)) deallocate(upha)
      if (allocated(vpha)) deallocate(vpha)

    enddo    !horizontal grids source

! Fermeture des fichiers
    ierr = stdclose(units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fs=',trim(fs)
      STOP
    endif
    ierr = stdclose(unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fd=',trim(fd)
      STOP
    endif

    if (status) then
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstharotate_status=FINISHED'
      close(99)
    endif

    end program cstharotate
