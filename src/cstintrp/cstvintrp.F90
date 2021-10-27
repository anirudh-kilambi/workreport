    program cstvintrp

!! Interpolates vertically data from a source file
!! to vertical grids from reference files
!! land-sea masks

!! Author: Francois Roy (2015), inspired by other tools from F. Dupont

!! Method:

!! Use convivial interface
!! Interpolate or spread data vertically
!! Assumes horizontal grids to match
!! Assumes depth read in (ip1) corresponds to center of model cell

!! Arguments:
!! 

!! Modules:
    USE std
    USE utils
    USE cdf
    USE cdf_std
    USE intrp_oce

    IMPLICIT NONE

!! Argument variables:

    character(len=1000) :: fs,fr,fd      !Source,reference,destination files,respectively
    character(len=4)   :: nomvar_src    !Nomvar selection in source file (default nomvar_src=' ')
    character(len=4)   :: nomvar_ref    !Reference field defining destination grid (in fr)
    integer            :: iargc,jarg,narg, units, unitr, unitd, nsprdh, nsprdv
    character(len=1000) :: carg,file_status
    logical            :: upperbound0s, upperbound0r, nochkg, ignmrf, corrbt

!! Local variables:

    logical status , maskin, maskou, chkmsk, chkgrid

    !... variable management
    type (std_variable),  &
       dimension(:),      &
       allocatable      :: vars,    &  !All variables from source file
                           varsg,   &  
                            !vars for one particular vertical grid
                           varsgc,  &  
                            !vars for one particular vertical grid and calendar
                           varsgct, &  
                            !vars for one particular vertical grid and calendar and time
                           varsgctn,&
                            !vars for one particular vertical grid and calendar and time and nomvar
                           varr, varrg
    type (std_variable) :: varwv, varwm, varwb, varwmb
    integer             :: nvs, nvr, nvsg, nvrg, nvsgc, nvsgct, nvsgctn
    integer             ::           ivsg, ivrg, ivsgc, ivsgct, ivsgctn

    type (std_vgrid),     &
       dimension(:),      &
       allocatable      :: vgrds, & !All grids from source file
                           vgrdr    !All grids from reference file corresponding to nomvar=nr
    integer             :: ngs, ngr, igs, igr
    integer             :: ierr, status_st90, id_m, i, j, ks, kd, it, inom, nnom, icnt, &
                           nxs, nys, nzs, nxd, nyd, nzd, nvartrt, ncal, ntmax, nt, ical
    integer             :: stdxtype, isprd, kup, kdn
    real(kind=8)        :: d0,d1,d2, a,b,c,d, Sum
    integer,              &
     dimension(:,:,:),    &
     allocatable        :: masks, masksr, maskd, maskdr, int4s
    integer,              &
     dimension(:,:),      &
     allocatable        :: datec
    real(kind=4),         &
     dimension(:,:,:),    &
     allocatable        :: flt4s,flt4d
    real(kind=8),         &
     dimension(:,:,:),    &
     allocatable        :: flt8s,flt8d
    real(kind=8),         &
     dimension(:,:),      &
     allocatable        :: wwgt, b2ds, b2dd
    integer,              &
     dimension(:),        &
     allocatable        :: masks_n, ntcal
    real(kind=4),         &
     dimension(:),        &
     allocatable        :: depths, depthd, wdepths, wdepthd, dzs, dzd
    character(len=4),     &
     dimension(:),        &
     allocatable        :: nomvar

! Definition des arguments par defaut

    nomvar_src=' '
    file_status=''
    upperbound0s=.true.
    upperbound0r=.true.
    nsprdh=0
    nsprdv=0
    nochkg=.false.
    ignmrf=.false.
    corrbt=.false.

! Traitement des arguments
    narg= iargc()
    IF ( narg < 8 .or. mod(narg,2).ne.0 ) THEN
       PRINT *, ' Usage : cstvintrp '
       PRINT *, ' -fs filesrc (source file to vertically interpolate)'
       PRINT *, ' -fr fileref (reference file defining destination grid)'
       PRINT *, ' -nr nomvar_ref (reference variable defining destination grid, in fr)'
       PRINT *, ' -fd filedst (destination file, interpolated data)'
       PRINT *, ' Optional arguments:'
       PRINT *, '   -ns nomvar_src (nomvar selection in src file,default all variables)'
       PRINT *, '   -status file_status (output program status in file_status for op needs)'
       PRINT *, '   -upbns0 T or F ( T: upper level of source grid assumed 0, F: extrapolated ...'
       PRINT *, '    ... note depths (ip1) are assumed centered model cell, default T )'
       PRINT *, '   -upbnr0 T or F ( T: upper level of reference grid assumed 0, F: extrapolated ...'
       PRINT *, '    ... note depths (ip1) are assumed centered model cell, default T )'
       PRINT *, '   -nsprdh nsprdh (number of iterations in horizontal spread around masked data, ...'
       PRINT *, '                   default => 0, no spread, NOTE: done before vertical spread)'
       PRINT *, '   -nsprdv nsprdv (number of iterations in vertical spread above/below masked data, ...'
       PRINT *, '                   default => 0, no spread)'
       PRINT *, '   -nochkg T or F ( T: no check for matching ig1, ig2, ... grtyp, just dimension check ...'
       PRINT *, '                    to match horizontal grids, default => F)'
       PRINT *, '   -ignmrf T or F ( T: ignore reference grid mask in the destination of data, default => F ...'
       PRINT *, '                    if T the destination mask is governed only by the spreaded source mask)'
       PRINT *, '   -corrbt T or F ( T: correct quantities (e.g. currents) based on the src/dst mask ...'
       PRINT *, '                    so that src and dst are equal in barotropic mode (2D avg), default => F )'
       STOP
    ENDIF
    jarg=1
    DO WHILE (jarg < narg )
      CALL getarg (jarg, carg)
      SELECT CASE ( carg)
      CASE ('-fs' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fs=TRIM(carg)
      CASE ('-fr' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fr=TRIM(carg)
      CASE ('-nr' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; nomvar_ref=TRIM(carg)
      CASE ('-fd' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fd=TRIM(carg)
      CASE ('-ns' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; nomvar_src=TRIM(carg)
      CASE ('-status' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; file_status=TRIM(carg)
         status=.true.
      CASE ('-upbns0' )
         jarg=jarg+1 ; CALL getarg(jarg,carg)
         if (TRIM(carg).eq.'F') then
           upperbound0s=.false.
          elseif (TRIM(carg).eq.'T') then
           upperbound0s=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -upbns0 T or upbns0 F' ; STOP
         endif
      CASE ('-upbnr0' )
         jarg=jarg+1 ; CALL getarg(jarg,carg)
         if (TRIM(carg).eq.'F') then
           upperbound0r=.false.
          elseif (TRIM(carg).eq.'T') then
           upperbound0r=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -upbnr0 T or upbnr0 F' ; STOP
         endif
      CASE ('-nochkg' )
         jarg=jarg+1 ; CALL getarg(jarg,carg)
         if (TRIM(carg).eq.'F') then
           nochkg=.false.
          elseif (TRIM(carg).eq.'T') then
           nochkg=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -nochkg T or nochkg F' ; STOP
         endif
      CASE ('-ignmrf' )
         jarg=jarg+1 ; CALL getarg(jarg,carg)
         if (TRIM(carg).eq.'F') then
           ignmrf=.false.
          elseif (TRIM(carg).eq.'T') then
           ignmrf=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -ignmrf T or ignmrf F' ; STOP
         endif
      CASE ('-corrbt' )
         jarg=jarg+1 ; CALL getarg(jarg,carg)
         if (TRIM(carg).eq.'F') then
           corrbt=.false.
          elseif (TRIM(carg).eq.'T') then
           corrbt=.true.
          else
           PRINT *,' Unknown value :', TRIM(carg) ; PRINT *,'Syntax: -corrbt T or corrbt F' ; STOP
         endif
      CASE ('-nsprdh' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) nsprdh
      CASE ('-nsprdv' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; read(carg,*) nsprdv
      CASE DEFAULT
         PRINT *,' Unknown option :', TRIM(carg) ; STOP
      END SELECT
      jarg=jarg+1
    ENDDO

    if (status) then
      status=.true.
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstvintrp_status=ABORT'
      close(99)
    endif

! Ouverture des fichiers

    ierr = stdopen(trim(fd),'new',unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fd=',trim(fd) 
      STOP 
    endif

    ierr = stdopen(trim(fr),'old',unitr,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdopen => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdopen => fr=',trim(fr) 
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
      write(*,*) 'CSTVINTRP getstdnvar => problem reading source file variables'
      write(*,*) 'CSTVINTRP getstdnvar => should have at least one variable to read'
      write(*,*) 'CSTVINTRP getstdnvar => nvs,status_st90=',nvs,status_st90
      STOP
    endif
    allocate(vars(nvs), varsg(nvs), varsgc(nvs), varsgct(nvs), varsgctn(nvs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP vars, varsg, ... allocation problem'
      write(*,*) 'CSTVINTRP nvs=',nvs
      STOP
    endif
    ierr = fillstdvar(vars,nvs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP fillstdvar => problem filling source variable informations'
      write(*,*) 'CSTVINTRP fillstdvar => ierr,status_st90=',ierr,status_st90
      STOP
    endif
! ... toutes les grilles verticales
    write(*,*) 'CSTVINTRP: SCANNING SOURCE VERTICAL GRIDS'
    ngs = getstdnvgrid(units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ngs <= 0 ) then 
      write(*,*) 'CSTVINTRP getstdnvgrid => problem reading source file grids'
      write(*,*) 'CSTVINTRP getstdnvgrid => should have at least one grid to read'
      write(*,*) 'CSTVINTRP getstdnvgrid => ngs,status_st90=',ngs,status_st90
      STOP
    endif
    allocate(vgrds(ngs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP vgrds allocation problem'
      write(*,*) 'CSTVINTRP ngs=',ngs
      STOP
    endif
    ierr = fillstdvgrid(vgrds,ngs,units,nomvar=nomvar_src,ierr_st90=status_st90)
    if ( ierr < 0 ) then
      write(*,*) 'CSTVINTRP fillstdvgrid => problem filling source file grid informations'
      write(*,*) 'CSTVINTRP fillstdvgrid => ierr,status_st90=',ierr,status_st90
      STOP
    endif

! Classement des grilles du fichier reference

! ... toutes les variables reference
    nvr = getstdnvar(unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    if ( nvs <= 0 ) then 
      write(*,*) 'CSTVINTRP getstdnvar => problem reading reference file variables'
      write(*,*) 'CSTVINTRP getstdnvar => should have at least one variable to read'
      write(*,*) 'CSTVINTRP getstdnvar => nvr,status_st90=',nvr,status_st90
      STOP
    endif
    allocate(varr(nvr), varrg(nvr), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP varr, varrg allocation problem'
      write(*,*) 'CSTVINTRP nvr=',nvr
      STOP
    endif
    ierr = fillstdvar(varr,nvr,unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP fillstdvar => problem filling reference variable informations'
      write(*,*) 'CSTVINTRP fillstdvar => ierr,status_st90=',ierr,status_st90
      STOP
    endif

! ... toutes les grilles (de nomvar_ref)
    write(*,*) 'CSTVINTRP: SCANNING DESTINATION GRIDS'
    ngr = getstdnvgrid(unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    if ( ngr <= 0 ) then 
      write(*,*) 'CSTVINTRP getstdnvgrid => problem reading reference file grids'
      write(*,*) 'CSTVINTRP getstdnvgrid => should have at least one grid to read'
      write(*,*) 'CSTVINTRP getstdnvgrid => ngr,status_st90=',ngr,status_st90
      STOP
    endif
    allocate(vgrdr(ngr), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTVINTRP vgrdr allocation problem'
      write(*,*) 'CSTVINTRP ngr=',ngr
      STOP
    endif
    ierr = fillstdvgrid(vgrdr,ngr,unitr,nomvar=nomvar_ref,ierr_st90=status_st90)
    if ( ierr < 0 ) then
      write(*,*) 'CSTVINTRP fillstdvgrid => problem filling source file grid informations'
      write(*,*) 'CSTVINTRP fillstdvgrid => ierr,status_st90=',ierr,status_st90
      STOP
    endif

! Boucle principale sur les grilles source et reference

    write(*,*) 'CSTVINTRP main loop ngr=',ngr
    write(*,*) 'CSTVINTRP main loop ngs=',ngs

    nvartrt=0
    do igr = 1, ngr
      do igs = 1, ngs
        ! Check if source destination horizontal grids match
        chkgrid=(vgrds(igs)%grtyp==vgrdr(igr)%grtyp.and. &
                 vgrds(igs)%ig1  ==vgrdr(igr)%ig1  .and. &
                 vgrds(igs)%ig2  ==vgrdr(igr)%ig2  .and. &
                 vgrds(igs)%ig3  ==vgrdr(igr)%ig3  .and. &
                 vgrds(igs)%ig4  ==vgrdr(igr)%ig4  .and. &
                 vgrds(igs)%ni   ==vgrdr(igr)%ni   .and. &
                 vgrds(igs)%nj   ==vgrdr(igr)%nj)
        if ( nochkg ) then
          if ( .not. chkgrid .and. &
                (vgrds(igs)%ni   ==vgrdr(igr)%ni   .and. &
                 vgrds(igs)%nj   ==vgrdr(igr)%nj) ) then
            write(*,*) 'CSTVINTRP: WARNING MATCHING A GRID ...'
            write(*,*) 'CSTVINTRP: NOT MATCHING IN TERMS OF STD FILE PARAMS, ig1, ig2, ...'
            write(*,*) 'CSTVINTRP: CAREFUL THE GRID ONLY MATCHES IN DIMENSIONS'
            chkgrid = .true.
          endif
        endif
        if ( chkgrid ) then
          ! Check unicity of reference grid set
          ierr = filtergvstdvar(varrg,nvrg,varr,nvr,vgrdr(igr))
          if ( ierr < 0 ) then
            write(*,*) 'CSTVINTRP filtergvstdvar => problem filtering grid information'
            write(*,*) 'CSTVINTRP filtergvstdvar => igr=',igr
            write(*,*) 'CSTVINTRP filtergvstdvar => vgrdr(igr)=',vgrdr(igr)
            STOP
          endif 
          if ( nochkg ) then
            ierr = filtergvstdvar(varsg,nvsg,vars,nvs,vgrds(igs))
            if ( ierr < 0 ) then
              write(*,*) 'CSTVINTRP filtergvstdvar => problem filtering grid information'
              write(*,*) 'CSTVINTRP filtergvstdvar => igs=',igs
              write(*,*) 'CSTVINTRP filtergvstdvar => vgrds(igs)=',vgrds(igs)
              STOP
            endif
            ierr = xferstdparpos(varsg,nvsg,units,unitd)
           else
            ierr = xferstdparpos(varrg,nvrg,unitr,unitd)
          endif
          if ( ierr /= 0 ) then
            write(*,*) 'CSTVINTRP xferstdparpos => problem transfering par pos parameters'
            write(*,*) 'CSTVINTRP vgrdr(igr) =', vgrdr(igr)
            STOP
          endif
          ncal = getstdncalendars(ntmax,varrg(1:nvrg),nvrg) 
          if ( ncal /= 1 .and. ntmax /= 1 ) then
            write(*,*) 'CSTVINTRP getstdncalendars => reference variable should have only one calendar'
            write(*,*) 'CSTVINTRP getstdncalendars => reference variable should have only one time frame'
            write(*,*) 'CSTVINTRP getstdncalendars => ncal,ntmax=',ncal,ntmax
            STOP
          endif
          ! Basic allocation
          nxs=vgrds(igs)%ni
          nys=vgrds(igs)%nj
          nzs=vgrds(igs)%nip1 
          nxd=nxs
          nyd=nys
          nzd=vgrdr(igr)%nip1
          allocate(depths(nzs),depthd(nzd),wdepths(nzs+1),wdepthd(nzd+1),wwgt(nzs,nzd))
          allocate(masks(nxs,nys,nzs),masksr(nxs,nys,nzs),masks_n(nzs), &
                   maskd(nxd,nyd,nzd),maskdr(nxd,nyd,nzd))
          if ( corrbt ) then
            allocate(dzs(nzs),dzd(nzd))
            allocate(b2ds(nxs,nys),b2dd(nxd,nyd))
          endif
          masks(:,:,:)=0; masksr(:,:,:)=0; maskd(:,:,:)=0; maskdr(:,:,:)=0; masks_n(:)=0
          icnt=0
          maskou=.false.
          ! Scan ref var and determine if there is a mask to read
          do ivrg=1,nvrg
            id_m=getstdmaskid(ivrg,varrg(1:nvrg),nvrg)
            if ( id_m /= ivrg ) then
              icnt=icnt+1
              kd=0
              do i=1,vgrdr(igr)%nip1
                if (vgrdr(igr)%ip1(i)==varrg(ivrg)%ip1) kd=i
              enddo
              if ( kd == 0 ) then
                write(*,*) 'CSTVINTRP: major structural problem'
                write(*,*) 'CSTVINTRP: varrg(ivrg)=',varrg(ivrg)
                STOP
              endif
              if ( id_m > 0 ) then
                ierr = getstdvar(maskd(:,:,kd),nxd,nyd,unitr, &
                          varrg(id_m),ierr_st90=status_st90)
                if ( ierr /= 0 ) then
                  write(*,*) 'CSTVINTRP getstdvar (maskd) => ierr,status_st90=',ierr,status_st90
                  write(*,*) 'CSTVINTRP getstdvar => varrg(id_m)%nomvar=',varrg(id_m)%nomvar
                  STOP
                endif
               else
                maskd(:,:,kd) = 1
              endif
              maskdr(:,:,kd) = maskd(:,:,kd)
            endif
          enddo
          ! Deals with destination grid now
          ierr = filtergvstdvar(varsg,nvsg,vars,nvs,vgrds(igs))
          if ( ierr < 0 ) then
            write(*,*) 'CSTVINTRP filtergvstdvar => problem filtering grid information'
            write(*,*) 'CSTVINTRP filtergvstdvar => igs=',igs
            write(*,*) 'CSTVINTRP filtergvstdvar => vgrds(igs)=',vgrds(igs)
            STOP
          endif
          ! Determines vertical depths weights from ip1 list
          write(*,*) '---------------------------------------'
          write(*,*) 'CSTVINTRP LIST OF SOURCE GRID DEPTHS =>'
          write(*,*) '---------------------------------------'
          do ks=1,nzs
            ierr = rpndepth2cdf(depths(ks),vgrds(igs)%ip1(ks))
            if ( ierr < 0 ) then
              write(*,*) 'CSTVINTRP rpndepth2cdf => problem with ip1 conversion'
              write(*,*) 'CSTVINTRP rpndepth2cdf => vgrds%ip1(ks)=', vgrds%ip1(ks)
              STOP
            endif
            write(*,*) 'ks,depths(ks)=',ks,depths(ks)
          enddo
          write(*,*) '--------------------------------------------'
          write(*,*) 'CSTVINTRP LIST OF DESTINATION GRID DEPTHS =>'
          write(*,*) '--------------------------------------------'
          do kd=1,nzd
            ierr = rpndepth2cdf(depthd(kd),vgrdr(igr)%ip1(kd))
            if ( ierr < 0 ) then
              write(*,*) 'CSTVINTRP rpndepth2cdf => problem with ip1 conversion'
              write(*,*) 'CSTVINTRP rpndepth2cdf => vgrdr%ip1(kd)=', vgrdr%ip1(kd)
              STOP
            endif
            write(*,*) 'kd,depthd(kd)=',kd,depthd(kd)
          enddo
          if (upperbound0s) then
            wdepths(1)=0.
           else  !extrapolation
            d1=depths(2)-depths(1)
            d2=d1
            if ( nzs >= 3 ) d2=depths(3)-depths(2)
            d0=2.*d1-d2
            wdepths(1)=depths(1) - 0.5*d0
          endif
          do ks=2,nzs+1
            wdepths(ks)=2.*depths(ks-1)-wdepths(ks-1)
          enddo
          if (upperbound0r) then
            wdepthd(1)=0.
           else  !extrapolation
            d1=depthd(2)-depthd(1)
            d2=d1
            if ( nzd >= 3 ) d2=depthd(3)-depthd(2)
            d0=2.*d1-d2
            wdepthd(1)=depthd(1) - 0.5*d0
          endif
          do kd=2,nzd+1
            wdepthd(kd)=2.*depthd(kd-1)-wdepthd(kd-1)
          enddo
          write(*,*) '-------------------------------------------'
          write(*,*) 'CSTVINTRP DEPTH UPPER BOUNDS FOR SRC GRID  '
          write(*,*) '-------------------------------------------'
          do ks=1,nzs+1
            write(*,*) 'ks, wdepths(ks)=',ks, wdepths(ks)
          enddo
          write(*,*) '-------------------------------------------'
          write(*,*) 'CSTVINTRP DEPTH UPPER BOUNDS FOR REF GRID  '
          write(*,*) '-------------------------------------------'
          do kd=1,nzd+1
            write(*,*) 'kd, wdepthd(kd)=',kd, wdepthd(kd)
          enddo
          if ( corrbt ) then
            write(*,*) '-------------------------------------------'
            write(*,*) 'CSTVINTRP DZ FOR SRC GRID  '
            write(*,*) '-------------------------------------------'
            do ks=1,nzs
              dzs(ks)=wdepths(ks+1)-wdepths(ks)
              write(*,*) 'ks, dzs(ks)=',ks, dzs(ks)
            enddo
            write(*,*) '-------------------------------------------'
            write(*,*) 'CSTVINTRP DZ FOR REF GRID  '
            write(*,*) '-------------------------------------------'
            do kd=1,nzd
              dzd(kd)=wdepthd(kd+1)-wdepthd(kd)
              write(*,*) 'kd, dzd(kd)=',kd, dzd(kd)
            enddo
          endif
          do kd=1,nzd
            write(*,*) '----------------------------------------------------'
            write(*,*) 'CSTVINTRP WEIGHTS FOR kd=',kd
            write(*,*) '----------------------------------------------------'
            do ks=1,nzs
              a=wdepths(ks)  -wdepthd(kd)
              b=wdepths(ks+1)-wdepthd(kd+1)
              c=wdepths(ks+1)-wdepthd(kd)
              d=wdepths(ks)  -wdepthd(kd+1)
              wwgt(ks,kd)=0.
              if (nzd==1.and.nzs==1) then
                wwgt(ks,kd)=1.
               else
                if      (d >  0.              ) then  !Cell is lower
                   wwgt(ks,kd)=0.
                 elseif (c <  0.              ) then  !Cell is higher
                   wwgt(ks,kd)=0.
                 elseif (a >= 0. .and. b <= 0.) then  !Cell is within
                   wwgt(ks,kd)=(wdepths(ks+1)-wdepths(ks))/(wdepthd(kd+1)-wdepthd(kd))
                 elseif (a >= 0. .and. b >  0.) then  !Cell is partially within
                   wwgt(ks,kd)=(wdepthd(kd+1)-wdepths(ks))/(wdepthd(kd+1)-wdepthd(kd))
                 elseif (a <  0. .and. b <= 0.) then  !Cell is partially within
                   wwgt(ks,kd)=(wdepths(ks+1)-wdepthd(kd))/(wdepthd(kd+1)-wdepthd(kd))
                 elseif (a <=  0. .and. b >= 0.) then  !Cell is all included
                   wwgt(ks,kd)=1.
                endif
              endif
              if (wwgt(ks,kd) /= 0.) then
                write(*,*) 'ks,w=',ks, wwgt(ks,kd)
              endif
            enddo
            write(*,*) '----------------------------------------------------'
          enddo
          ncal = getstdncalendars(ntmax,varsg(1:nvsg),nvsg)
          write(*,*) 'CSTVINTRP: number of calendars for this grid=',ncal
          if ( ncal <= 0 ) then
            write(*,*) 'CSTVINTRP getstdncalendars => problem no calendar found for this list'
            write(*,*) 'CSTVINTRP getstdncalendars => ncal,ntmax=',ncal,ntmax
            STOP
          endif
          allocate(datec(ncal,ntmax),ntcal(ncal))
          ierr = fillstdcalendars(datec,ntcal,ncal,ntmax,varsg(1:nvsg),nvsg)
          if ( ierr < 0 ) then
            write(*,*) 'CSTVINTRP fillstdcalendars => problem filling calendar'
            write(*,*) 'CSTVINTRP fillstdcalendars => ncal,ntmax=',ncal,ntmax
            STOP
          endif
          do ical=1,ncal
            ierr = filtercalstdvar(varsgc(1:nvsg),nvsgc,varsg(1:nvsg),nvsg, &
                                   datec(ical,1:ntcal(ical)),ntcal(ical))
            if ( ierr < 0 ) then
              write(*,*) 'CSTVINTRP filtercalstdvar => problem extracting calendar sub-ensemble'
              write(*,*) 'CSTVINTRP filtercalstdvar => ierr=',ierr
              STOP
            endif
            nt=ntcal(ical)
            do it=1,nt
              ierr = filtertstdvar(varsgct(1:nvsgc),nvsgct,varsgc(1:nvsgc),nvsgc,datec(ical,it))
              nnom=getnnomvar(varsgct(1:nvsgct),nvsgct)
              if ( nnom < 0 ) then
                write(*,*) 'CSTVINTRP getnnomvar => problem getting nomvar list'
                STOP
              endif
              if ( nnom > 0 ) then
                allocate(nomvar(nnom))
                ierr=fillnomvar(nomvar,nnom,varsgct(1:nvsgct),nvsgct)
                do inom=1,nnom
                  ierr = filternomvar(varsgctn(1:nvsgct),nvsgctn, &
                                      varsgct(1:nvsgct),nvsgct,nomvar(inom))
                  if ( ierr < 0 ) then
                    write(*,*) 'CSTVINTRP filternomvar => problem isolating nomvar'
                    write(*,*) 'CSTVINTRP filternomvar => inom,nomvar(inom)=',inom,nomvar(inom)
                    STOP
                  endif
                  icnt=0
                  maskin=.false.
                  do ivsgctn=1,nvsgctn
                    id_m=getstdmaskid(ivsgctn,varsgctn(1:nvsgctn),nvsgctn)
                    if ( id_m /= ivsgctn ) then
                      icnt=icnt+1
                      ks=0
                      do i=1,vgrds(igs)%nip1
                        if (vgrds(igs)%ip1(i)==varsgctn(ivsgctn)%ip1) ks=i
                      enddo
                      if (ks==1) varwv=varsgctn(ivsgctn)
                      if ( ks == 0 ) then
                        write(*,*) 'CSTVINTRP: major structural problem'
                        write(*,*) 'CSTVINTRP: varsgctn(ivsgctn)=',varsgctn(ivsgctn)
                        STOP
                      endif
                      if ( varsgctn(ivsgctn)%xtype == 3 ) then
                        stdxtype=3
                        if (.not. allocated(flt8s) ) allocate(flt8s(nxs,nys,nzs))
                        ierr = getstdvar(flt8s(:,:,ks),nxs,nys,units, &
                                  varsgctn(ivsgctn),ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP getstdvar (flt8s) => ierr,status_st90=',ierr,status_st90
                          write(*,*) 'CSTVINTRP getstdvar => varsgctn(ivsgctn)%nomvar=',varsgctn(ivsgctn)%nomvar
                          STOP
                        endif
                       elseif ( varsgctn(ivsgctn)%xtype == 2 ) then
                        stdxtype=2
                        if (.not. allocated(flt4s)) allocate(flt4s(nxs,nys,nzs))
                        ierr = getstdvar(flt4s(:,:,ks),nxs,nys,units, &
                                  varsgctn(ivsgctn),ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP getstdvar (flt4s) => ierr,status_st90=',ierr,status_st90
                          write(*,*) 'CSTVINTRP getstdvar => varsgctn(ivsgctn)%nomvar=',varsgctn(ivsgctn)%nomvar
                          STOP
                        endif
                       else
                        stdxtype=varsgctn(ivsgctn)%xtype
                        if (.not. allocated(int4s) ) allocate(int4s(nxs,nys,nzs))
                        if (.not. allocated(flt4s) ) allocate(flt4s(nxs,nys,nzs))
                        ierr = getstdvar(int4s(:,:,ks),nxs,nys,units, &
                                  varsgctn(ivsgctn),ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP getstdvar (int4s) => ierr,status_st90=',ierr,status_st90
                          write(*,*) 'CSTVINTRP getstdvar => varsgctn(ivsgctn)%nomvar=',varsgctn(ivsgctn)%nomvar
                          STOP
                        endif
                        flt4s(:,:,ks:ks)=int4s(:,:,ks:ks)
                      endif
                    endif
                    if ( id_m > 0 ) then
                      maskin = .true.
                      ierr = getstdvar(masks(:,:,ks),nxs,nys,units, &
                                varsgctn(id_m),ierr_st90=status_st90)
                      if ( ierr /= 0 ) then
                        write(*,*) 'CSTVINTRP getstdvar (masks) => ierr,status_st90=',ierr,status_st90
                        write(*,*) 'CSTVINTRP getstdvar => varsgctn(id_m)%nomvar=',varsgctn(id_m)%nomvar
                        STOP
                      endif
                      if (ks==1) varwm=varsgctn(id_m)
                     else
                      masks(:,:,ks) = 1
                    endif
                    masksr(:,:,ks) = masks(:,:,ks)
                  enddo
                  if (icnt /= nzs) then
                    write(*,*) 'CSTVINTRP inconsistent number of records ...'
                    write(*,*) 'CSTVINTRP according to this vertical grid=> vgrds(igs)=',vgrds(igs)
                    write(*,*) 'CSTVINTRP icnt =',icnt
                    write(*,*) 'CSTVINTRP nvsgctn, nzs =',nvsgctn, nzs
                    STOP
                  endif
                  ! Compute source barotropic quantities if needed
                  if ( corrbt ) then
                    if ( stdxtype == 3 ) then
                      do j=1,nys
                      do i=1,nxs
                        b2ds(i,j)=0.
                        Sum=0.
                        do ks=1,nzs
                          b2ds(i,j) = b2ds(i,j) + dzs(ks)*flt8s(i,j,ks)*masksr(i,j,ks)
                          Sum=Sum + dzs(ks)*masksr(i,j,ks)
                        enddo
                        if ( Sum > 0. ) b2ds(i,j)=b2ds(i,j)/Sum
                      enddo
                      enddo
                     else
                      do j=1,nys
                      do i=1,nxs
                        b2ds(i,j)=0.
                        Sum=0.
                        do ks=1,nzs
                          b2ds(i,j) = b2ds(i,j) + dzs(ks)*flt4s(i,j,ks)*masksr(i,j,ks)
                          Sum=Sum + dzs(ks)*masksr(i,j,ks)
                        enddo
                        if ( Sum > 0. ) b2ds(i,j)=b2ds(i,j)/Sum
                      enddo
                      enddo
                    endif
                  endif
                  ! Main interpolation extrapolation loop
                  if ( stdxtype == 3 ) then
                    if (.not. allocated(flt8d) ) allocate(flt8d(nxd,nyd,nzd))
                    flt8d(:,:,:) = 0.
                    if (nsprdh>0) then
                      do ks=1,nzs
                        call mapospread_simple(flt8s(:,:,ks),masks(:,:,ks),nxs,nys,nsprdh)
                      enddo
                    endif
                    do j=1,nyd
                    do i=1,nxd
                      if (maskin) then
                        masks_n(:)=masks(i,j,:)
                        do isprd=1,nsprdv
                          do ks=1,nzs
                            if ( masks(i,j,ks) == 0 ) then
                              kup=max(1,ks-1)
                              kdn=min(nzs,ks+1)
                              Sum=masks(i,j,kup) + masks(i,j,kdn)
                              if ( Sum > 0. ) then
                                flt8s(i,j,ks)=(masks(i,j,kup)*flt8s(i,j,kup)+ &
                                               masks(i,j,kdn)*flt8s(i,j,kdn))/Sum
                                masks_n(ks)=1
                              endif
                            endif
                          enddo
                          masks(i,j,:)=masks_n(:)
                        enddo
                      endif
                      do kd=1,nzd
                        if ( maskd(i,j,kd) > 0 .or. ignmrf ) then
                          flt8d(i,j,kd)=0.
                          Sum=0.
                          do ks=1,nzs
                            if ( masks(i,j,ks) > 0 ) then
                              Sum=Sum+wwgt(ks,kd)
                              flt8d(i,j,kd)=flt8d(i,j,kd)+wwgt(ks,kd)*flt8s(i,j,ks)
                            endif
                          enddo
                          if (Sum.gt.0.) then
                            flt8d(i,j,kd)=flt8d(i,j,kd)/Sum
                            maskd(i,j,kd) = 1
                           else
                            maskd(i,j,kd) = 0
                            maskou = .true.
                          endif
                        endif
                      enddo
                    enddo
                    enddo
                   else  ! flt4d and int4d case
                    if (.not. allocated(flt4d) ) allocate(flt4d(nxd,nyd,nzd))
                    flt4d(:,:,:) = 0.
                    if (nsprdh>0) then
                      do ks=1,nzs
                        call mapospread_simple(flt4s(:,:,ks),masks(:,:,ks),nxs,nys,nsprdh)
                      enddo
                    endif
                    do j=1,nyd
                    do i=1,nxd
                      if (maskin) then
                        masks_n(:)=masks(i,j,:)
                        do isprd=1,nsprdv
                          do ks=1,nzs
                            if ( masks(i,j,ks) == 0 ) then
                              kup=max(1,ks-1)
                              kdn=min(nzs,ks+1)
                              Sum=masks(i,j,kup) + masks(i,j,kdn)
                              if ( Sum > 0. ) then
                                flt4s(i,j,ks)=(masks(i,j,kup)*flt4s(i,j,kup)+ &
                                               masks(i,j,kdn)*flt4s(i,j,kdn))/Sum
                                masks_n(ks)=1
                              endif
                            endif
                          enddo
                          masks(i,j,:)=masks_n(:)
                        enddo
                      endif
                      do kd=1,nzd
                        if ( maskd(i,j,kd) > 0 .or. ignmrf ) then
                          flt4d(i,j,kd)=0.
                          Sum=0.
                          do ks=1,nzs
                            if ( masks(i,j,ks) > 0 ) then
                              Sum=Sum+wwgt(ks,kd)
                              flt4d(i,j,kd)=flt4d(i,j,kd)+wwgt(ks,kd)*flt4s(i,j,ks)
                            endif
                          enddo
                          if (Sum.gt.0.) then
                            flt4d(i,j,kd)=flt4d(i,j,kd)/Sum
                            maskd(i,j,kd) = 1
                           else
                            maskd(i,j,kd) = 0
                            maskou = .true.
                          endif
                        endif
                      enddo
                    enddo
                    enddo
                  endif
                  ! Compute destination barotropic quantities if needed
                  ! and adjust destination fields
                  if ( corrbt ) then
                    if ( stdxtype == 3 ) then
                      do j=1,nyd
                      do i=1,nxd
                        b2dd(i,j)=0.
                        Sum=0.
                        do kd=1,nzd
                          b2dd(i,j) = b2dd(i,j) + dzd(kd)*flt8d(i,j,kd)*maskdr(i,j,kd)
                          Sum=Sum + dzd(kd)*maskdr(i,j,kd)
                        enddo
                        if ( Sum > 0. ) then
                          b2dd(i,j)=b2dd(i,j)/Sum
                          ! Compute correction and put it in b2dd
                          b2dd(i,j) = b2ds(i,j) - b2dd(i,j)
                          ! Apply it to flt8d 
                          if ( masksr(i,j,1) > 0 ) then
                            do kd=1,nzd
                              flt8d(i,j,kd) = flt8d(i,j,kd) + b2dd(i,j)*dzd(kd)*maskdr(i,j,kd)/Sum
                            enddo
                           else
                            b2dd(i,j) = 0.
                          endif
                        endif  
                      enddo
                      enddo
                     else
                      do j=1,nyd
                      do i=1,nxd
                        b2dd(i,j)=0.
                        Sum=0.
                        do kd=1,nzd
                          b2dd(i,j) = b2dd(i,j) + dzd(kd)*flt4d(i,j,kd)*maskdr(i,j,kd)
                          Sum=Sum + dzd(kd)*maskdr(i,j,kd)
                        enddo
                        if ( Sum > 0. ) then
                          b2dd(i,j)=b2dd(i,j)/Sum
                          ! Compute correction and put it in b2dd
                          b2dd(i,j) = b2ds(i,j) - b2dd(i,j)
                          ! Apply it to flt4d 
                          if ( masksr(i,j,1) > 0 ) then
                            do kd=1,nzd
                              flt4d(i,j,kd) = flt4d(i,j,kd) + b2dd(i,j)*dzd(kd)*maskdr(i,j,kd)/Sum
                            enddo
                           else
                            b2dd(i,j) = 0.
                          endif
                        endif  
                      enddo
                      enddo
                    endif
                  endif
                  ! Write phase
                  if ( maskou .and. .not. maskin ) then
                    varwv%typvar=varwv%typvar(1:1)//'@'
                    varwm=varwv
                    varwm%typvar='@@'
                    varwm%nbits=1
                    varwm%datyp=2
                    varwm%xtype=1
                  endif
                  if ( stdxtype == 3 ) then
                    ! write case flt8d
                    do kd=1,nzd
                      varwv%ip1=vgrdr(igr)%ip1(kd)
                      ierr = putstdvar(unitd, flt8d(:,:,kd),  &
                           nxd, nyd, varwv, ierr_st90=status_st90)
                      if ( ierr /= 0 ) then
                        write(*,*) 'CSTVINTRP putstdvar problem => flt8d'
                        write(*,*) 'CSTVINTRP status_st90=',status_st90
                        write(*,*) 'CSTVINTRP ierr,varwv%ip1=',varwv%ip1
                        write(*,*) 'CSTVINTRP varwv=',varwv
                        STOP
                      endif
                      if ( corrbt .and. kd == 1 ) then
                        varwb=varwv
                        varwb%nomvar='BCOR'
                        ierr = putstdvar(unitd, b2dd(:,:),  &
                             nxd, nyd, varwb, ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP putstdvar problem => b2dd (8)'
                          write(*,*) 'CSTVINTRP status_st90=',status_st90
                          write(*,*) 'CSTVINTRP varwb=',varwb
                          STOP
                        endif
                      endif
                    enddo
                   else
                    ! write case flt4d and int4d
                    if ( varwv%xtype == 1 ) then
                      varwv%xtype=2
                      varwv%nbits=32
                      varwv%datyp=5
                      varwv%xtype=2
                    endif
                    do kd=1,nzd
                      varwv%ip1=vgrdr(igr)%ip1(kd)
                      ierr = putstdvar(unitd, flt4d(:,:,kd),  &
                           nxd, nyd, varwv, ierr_st90=status_st90)
                      if ( ierr /= 0 ) then
                        write(*,*) 'CSTVINTRP putstdvar problem => flt4d'
                        write(*,*) 'CSTVINTRP status_st90=',status_st90
                        write(*,*) 'CSTVINTRP ierr,varwv%ip1=',varwv%ip1
                        write(*,*) 'CSTVINTRP varwv=',varwv
                        STOP
                      endif
                      if ( corrbt .and. kd == 1 ) then
                        varwb=varwv
                        varwb%nomvar='BCOR'
                        ierr = putstdvar(unitd, real(b2dd(:,:)),  &
                             nxd, nyd, varwb, ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP putstdvar problem => b2dd'
                          write(*,*) 'CSTVINTRP status_st90=',status_st90
                          write(*,*) 'CSTVINTRP varwb=',varwb
                          STOP
                        endif
                      endif
                    enddo
                  endif
                  if (maskou) then
                    do kd=1,nzd
                      varwm%ip1=vgrdr(igr)%ip1(kd)
                      ierr = putstdvar(unitd, maskd(:,:,kd),  &
                           nxd, nyd, varwm, ierr_st90=status_st90)
                      if ( ierr /= 0 ) then
                        write(*,*) 'CSTVINTRP putstdvar problem => maskd'
                        write(*,*) 'CSTVINTRP status_st90=',status_st90
                        write(*,*) 'CSTVINTRP ierr,varwm%ip1=',varwm%ip1
                        write(*,*) 'CSTVINTRP varwm=',varwm
                        STOP
                      endif
                      if ( corrbt .and. kd == 1 ) then
                        varwmb=varwm
                        varwmb%nomvar='BCOR'
                        ierr = putstdvar(unitd, maskd(:,:,kd),  &
                             nxd, nyd, varwmb, ierr_st90=status_st90)
                        if ( ierr /= 0 ) then
                          write(*,*) 'CSTVINTRP putstdvar problem => b2dd (mask)'
                          write(*,*) 'CSTVINTRP status_st90=',status_st90
                          write(*,*) 'CSTVINTRP varwmb=',varwmb
                          STOP
                        endif
                      endif
                    enddo
                    chkmsk=.false.
                    do kd=1,nzd
                      do j=1,nyd
                      do i=1,nxd
                        if (maskd(i,j,kd) /= maskdr(i,j,kd)) chkmsk=.true.
                      enddo
                      enddo
                    enddo
                    if ( chkmsk ) then
                      write(*,*) 'CSTVINTRP Warning destination mask does not match reference mask'
                      write(*,*) 'CSTVINTRP may need to increase nsprdh and/or nsprdv'
                    endif
                  endif
                  nvartrt=nvartrt+1
                enddo
                deallocate(nomvar)
              endif   ! nnom > 0
            enddo     ! nt
          enddo       ! ncal
          deallocate(depths,depthd,wdepths,wdepthd,wwgt,datec,ntcal, &
                     masks,masksr,masks_n,maskd,maskdr)
          if (allocated(int4s)) deallocate(int4s)
          if (allocated(flt4s)) deallocate(flt4s)
          if (allocated(flt8s)) deallocate(flt8s)
          if (allocated(flt4d)) deallocate(flt4d)
          if (allocated(flt8d)) deallocate(flt8d)
          if (allocated(b2ds))  deallocate(b2ds)
          if (allocated(b2dd))  deallocate(b2dd)
          if (allocated(dzs))   deallocate(dzs)
          if (allocated(dzd))   deallocate(dzd)
        endif  !Horizontal grid match
      enddo    !vertical grids source
    enddo      !vertical grids reference

! Fermeture des fichiers
    ierr = stdclose(units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fs=',trim(fs)
      STOP
    endif
    ierr = stdclose(unitr,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fr=',trim(fr)
      STOP
    endif
    ierr = stdclose(unitd,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'stdclose => ierr,status_st90=',ierr,status_st90
      write(*,*) 'stdclose => fd=',trim(fd)
      STOP
    endif

    write(*,*) 'CSTVINTRP: FINISHED, TREATED  N FIELDS = ', nvartrt
    if (status) then
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstvintrp_status=FINISHED'
      close(99)
    endif

    end program cstvintrp
