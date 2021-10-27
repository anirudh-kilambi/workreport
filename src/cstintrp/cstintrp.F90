    program cstintrp

!! Customized interpolator for ice-ocean grids managing 
!! land-sea masks

!! Author: Francois Roy (November 2009)
!! v_3.0.0 - Re-organized + bi-cubic and SCRIP option (2015)
!!         - Support for Yin-Yang
!! (F. Roy + F. Dupont)

!! Method:

!! Location in lat-lon and xy plans:
!! Use map projection routines from ocean grids when not standard RPN grid
!! Use ezscint xy to lat-lon routines (only) when standard RPN grid

!! Mixt regridding (intyp=mixt):
!! Adjust locally, if equal or more than 3 points per cell
!! in neighbourhood i-(nwgt+0.5),i+(nwgt+0.5),j-(nwgt+0.5),j+(nwgt+0.5),
!! will apply nearest neighbouring aggregation
!! with cell surface weighting, otherwise will
!! apply bi-linear interpolation

!! Bi-Linear interpolator (intyp=bilin):
!! Use a general 2D interpolation routine working
!! for any 4 points polygon in xy plan, therefore
!! considers local deformation of lat-lon grids
!! after locally centered PS projection

!! Bi-cubic interpolation (intyp=bicub):
!! To be complted , Fred...

!! SCRIP conservative interpolation (intyp=scrip):
!! To be complted , Fred...

!! Mixt regridding (intyp=mixtbc):
!! Aggregation combined with bi-cubic
!! To be complted , F. Roy

!! For any interpolation point, project to a PS
!! grid with center at interpolation point
!! The deformation remains very small unless
!! the resolution of the grid is poor
!! Example: For orca025 the deformation error is 
!! maximum at equator (1.000010)
!! where the resolution is 1/4 deg.

!! Mask management:
!! When there is a land-sea mask, the mask is interpolated
!! An option allows to decide of the treshold 
!! applied to the masked interpolated data

!! Data spreading before interpolation with mask:
!! A spreader program (mapospread)
!! is used to spread data around mask transition areas
!! Applying a strict treshold on interpolated mask 
!! will avoid artefacts related to a low value of
!! nsprd

!! Arguments:
!! 

!! Modules:
    USE std
    USE utils
    USE intrp_itfc
    USE mympp
    USE cstintrp_mod

    !---------------------------------------------------------------------------------------
    IMPLICIT NONE
    !! Local variables:
    integer igs, igr, nvartrt
    !---------------------------------------------------------------------------------------

    CALL init_general_param
    CALL print_usage
    CALL check_arguments
    CALL read_vector_table
    CALL open_src_ref_dst_files
    CALL read_list_var_grid_src
    CALL read_list_var_grid_ref

! Boucle principale sur les grilles source et reference

    nvartrt = 0
    do igr = 1, ngr
      do igs = 1, ngs
        CALL treat_each_grid_src_dst(grds(igs),grdr(igr),nvartrt)
      enddo
    enddo

    CALL close_all_files
    write(*,*) 'CSTINTRP: FINISHED, TREATED  N FIELDS OR VECTORS = ', nvartrt

    ! deprecated option, only necessary if the calling script does not stop when encountering an error
    if ( file_status /= '' ) then
      write(*,*) '!!! YOU ARE USING A DEPRECATED FEATURE: FILE_STATUS WILL BE REMOVED IN NEXT VERSION !!!'
      open(99,file=file_status,form='formatted')
      write(99,*) 'cstintrp_status=FINISHED'
      close(99)
    endif
    end program cstintrp
