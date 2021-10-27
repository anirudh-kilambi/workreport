    program cstz1z22ap

!! Convert harmonic analysis output (z1==>x,z2==>y) from NEMO into
!! amplitudes and phases (A, phi)

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

    IMPLICIT NONE

!! Argument variables:

    character(len=1000) :: fs,fd,fl      !Source,destination files,respectively
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
    real(kind=4), parameter  :: epsilon = 1.e-5
    integer,              &
     dimension(:,:),      &
     allocatable        :: masks
    real(kind=4),         &
     dimension(:,:),      &
     allocatable        :: xs,ys,zamp,zpha
    integer, parameter  :: ncmax=1000
    character(len=4),  &
       dimension(ncmax) :: nvx, nvy, nva, nvp
    integer             :: ic, nc, ii, im
    real(kind=4)        :: pi, torad

! Definition des arguments par defaut

    file_status=''

    fs='/bidon'
    fd='/bidon'
    fl='/bidon'

! Traitement des arguments
    narg= iargc()
    IF ( narg < 6 .or. mod(narg,2).ne.0 ) THEN
       PRINT *, ' Usage : cstz1z22ap '
       PRINT *, ' -fs filesrc (source file to convert)'
       PRINT *, ' -fd filedst (destination file)'
       PRINT *, ' -fl filelst (file with list of constituent names)'
       PRINT *, ' Optional arguments:'
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
      CASE ('-fl' )
         jarg=jarg+1 ; CALL getarg(jarg,carg) ; fl=TRIM(carg)
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
      open(99,file=file_status,form='formatted', status='old')
      write(99,*) 'cstz1z22ap_status=ABORT'
      close(99)
    endif

! Ouverture des fichiers

    open(98,file=fl,form='formatted')
    nc=0
    DO
      if (nc+1 > ncmax) STOP 'Not enough constituent storage'
      read(98,'(a4,x,a4)',end=999) nvx(nc+1), nvy(nc+1)
      nc=nc+1
    ENDDO

999 continue
    if ( nc == 0 ) STOP 'No constituent to search for'

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
    nvs = getstdnvar(units,ierr_st90=status_st90)
    if ( nvs <= 0 ) then 
      write(*,*) 'CSTZ1Z22AP getstdnvar => problem reading source file variables'
      write(*,*) 'CSTZ1Z22AP getstdnvar => should have at least one variable to read'
      write(*,*) 'CSTZ1Z22AP getstdnvar => nvs,status_st90=',nvs,status_st90
      STOP
    endif
    allocate(vars(nvs), varsg(nvs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTZ1Z22AP vars, varsg, ... allocation problem'
      write(*,*) 'CSTZ1Z22AP nvs=',nvs
      STOP
    endif
    ierr = fillstdvar(vars,nvs,units,ierr_st90=status_st90)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTZ1Z22AP fillstdvar => problem filling source variable informations'
      write(*,*) 'CSTZ1Z22AP fillstdvar => ierr,status_st90=',ierr,status_st90
      STOP
    endif
! ... toutes les grilles horizontales
    write(*,*) 'CSTZ1Z22AP: SCANNING SOURCE HORIZONTAL GRIDS'
    ngs = getstdngrid(units,ierr_st90=status_st90)
    if ( ngs <= 0 ) then 
      write(*,*) 'CSTZ1Z22AP getstdngrid => problem reading source file grids'
      write(*,*) 'CSTZ1Z22AP getstdngrid => should have at least one grid to read'
      write(*,*) 'CSTZ1Z22AP getstdngrid => ngs,status_st90=',ngs,status_st90
      STOP
    endif
    allocate(hgrds(ngs), stat=ierr)
    if ( ierr /= 0 ) then
      write(*,*) 'CSTZ1Z22AP hgrds allocation problem'
      write(*,*) 'CSTZ1Z22AP ngs=',ngs
      STOP
    endif
    ierr = fillstdgrid(hgrds,ngs,units,ierr_st90=status_st90)
    if ( ierr < 0 ) then
      write(*,*) 'CSTZ1Z22AP fillstdgrid => problem filling source file grid informations'
      write(*,*) 'CSTZ1Z22AP fillstdgrid => ierr,status_st90=',ierr,status_st90
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
        write(*,*) 'CSTZ1Z22AP filtergstdvar => problem filtering grid information'
        write(*,*) 'CSTZ1Z22AP filtergstdvar => igs=',igs
        write(*,*) 'CSTZ1Z22AP filtergstdvar => hgrds(igs)=',hgrds(igs)
        STOP
      endif
      ierr = xferstdparpos(varsg,nvsg,units,unitd)
      if ( ierr /= 0 ) then
        write(*,*) 'CSTZ1Z22AP xferstdparpos => problem transfering par pos parameters'
        write(*,*) 'CSTZ1Z22AP hgrds(igs) =', hgrds(igs)
        STOP
      endif
      ! Basic allocation
      nxs=hgrds(igs)%ni
      nys=hgrds(igs)%nj
      nxd=nxs
      nyd=nys

      write(*,*) 'CSTZ1Z22AP treating grid with dimensions nxs,nys=',nxs,nys
      write(*,*) 'CSTZ1Z22AP treating grid with dimensions hgrds(igs)=',hgrds(igs)

      allocate(masks(nxs,nys))
      allocate(xs(nxs,nys),ys(nxs,nys),zamp(nxs,nys),zpha(nxs,nys))

      masks(:,:)=1 

      do ic=1,nc

        write(*,*) 'CSTZ1Z22AP TREATING ' , trim(nvx(ic)), ' ', trim(nvy(ic))

        do ivsg = 1, nvsg

          if ( varsg(ivsg)%nomvar == trim(nvx(ic)) ) then
       
            if ( varsg(ivsg)%typvar /= '@@' ) then
              if (.not. allocated(zamp)) allocate(zamp(nxs,nys))
              if (.not. allocated(zpha)) allocate(zpha(nxs,nys))
              ierr = getstdvar(xs(:,:),nxs,nys,units, &
                               varsg(ivsg),ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP getstdvar (XS) => ierr,status_st90=',ierr,status_st90
                write(*,*) 'CSTZ1Z22AP getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
                STOP
              endif
              varwg = varsg(ivsg)
              varwg%nomvar=trim(nvy(ic))
              ierr = getstdvar(ys(:,:),nxs,nys,units, &
                               varwg,ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP getstdvar (YS) => ierr,status_st90=',ierr,status_st90
                write(*,*) 'CSTZ1Z22AP getstdvar => varwg%nomvar=',varwg%nomvar
                STOP
              endif
              !ys(:,:) = -ys(:,:)  !From NEMO to OSU convention
              !NEMO tide analysis uses OSU convention
              do j=1,nys
              do i=1,nxs
                zpha(i,j)=atan2(-ys(i,j),xs(i,j))
                if (abs(cos(zpha(i,j))) < epsilon ) then
                  zamp(i,j)=-ys(i,j)/sin(zpha(i,j))
                 else
                  zamp(i,j)= xs(i,j)/cos(zpha(i,j))
                endif
                zpha(i,j)=zpha(i,j)/torad
              enddo
              enddo
              do ii=1,len(nvx(ic))
                if (nvx(ic)(ii:ii) /= nvy(ic)(ii:ii)) then
                  im=ii
                endif
              enddo
              varwg%nomvar=nvx(ic)
              varwg%nomvar(im:im)='A'
              ierr = putstdvar(unitd, zamp(:,:),  &
                               nxd, nyd, varwg, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP putstdvar problem (ZAMP)'
                write(*,*) 'CSTZ1Z22AP status_st90=',status_st90
                write(*,*) 'CSTZ1Z22AP varsg(ivsg)=',varsg(ivsg)
                STOP
              endif
              varwg%nomvar=nvx(ic)
              varwg%nomvar(im:im)='P'
              ierr = putstdvar(unitd, zpha(:,:),  &
                               nxd, nyd, varwg, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP putstdvar problem (ZPHA)'
                write(*,*) 'CSTZ1Z22AP status_st90=',status_st90
                write(*,*) 'CSTZ1Z22AP varsg(ivsg)=',varsg(ivsg)
              STOP
              endif
             else
              ierr = getstdvar(masks(:,:),nxs,nys,units, &
                               varsg(ivsg),ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP getstdvar (ZXS MASKS) => ierr,status_st90=',ierr,status_st90
                write(*,*) 'CSTZ1Z22AP getstdvar => varsg(ivsg)%nomvar=',varsg(ivsg)%nomvar
                STOP
              endif
              do ii=1,len(nvx(ic))
                if (nvx(ic)(ii:ii) /= nvy(ic)(ii:ii)) then
                  im=ii
                endif
              enddo
              varwg=varsg(ivsg)
              varwg%nomvar=nvx(ic)
              varwg%nomvar(im:im)='A'
              ierr = putstdvar(unitd, masks(:,:),  &
                               nxd, nyd, varwg, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP putstdvar problem (ZAMP MASKS)'
                write(*,*) 'CSTZ1Z22AP status_st90=',status_st90
                write(*,*) 'CSTZ1Z22AP varsg(ivsg)=',varsg(ivsg)
                STOP
              endif
              varwg%nomvar(im:im)='P'
              ierr = putstdvar(unitd, masks(:,:),  &
                               nxd, nyd, varwg, ierr_st90=status_st90)
              if ( ierr /= 0 ) then
                write(*,*) 'CSTZ1Z22AP putstdvar problem (ZPHA MASKS)'
                write(*,*) 'CSTZ1Z22AP status_st90=',status_st90
                write(*,*) 'CSTZ1Z22AP varsg(ivsg)=',varsg(ivsg)
                STOP
              endif
            endif
          endif

        enddo

      enddo

      deallocate(masks,xs,ys,zamp,zpha)

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
      write(99,*) 'cstz1z22ap_status=FINISHED'
      close(99)
    endif

    close(98)

    end program cstz1z22ap
