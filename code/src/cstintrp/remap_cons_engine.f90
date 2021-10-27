module remap_cons_engine

use grids
use intrp_oce
use kdtree2_module
use mesh

implicit none

  integer :: num_wts = 1, resize_increment = 1000

  TYPE (type_mesh), pointer :: src_grid => NULL() &
                             , opp_grid => NULL()

  integer, allocatable, dimension(:) :: list_opp_ele_ref, list_src_ele_ref
  real(kind=8) :: tiny = epsilon
  type (type_mesh), target :: mesh_opp_loc, mesh_src_loc
  integer :: n_search_scrip = 20
  TYPE(type_weight_unstr) :: uwgts

contains



!**********************************************************************
!**********************************************************************

subroutine grid_integrate

!-----------------------------------------------------------------------
!
!     do the area integral for grid1 and grid2 (2 sweeps)
!
!-----------------------------------------------------------------------

integer cell_add, ilink, nlink, cell_opp
integer ncorners
integer max_link

!-----------------------------------------------------------------------
!
! find weights: grid1 as source grid, grid2 as opposed grid (first sweep)
!
!-----------------------------------------------------------------------

      allocate( src_grid % geo_p % area        (ne), &
     &          src_grid % geo_p % frac        (ne), &
     &          src_grid % geo_p % centroid_lat(ne), &
     &          src_grid % geo_p % centroid_lon(ne))

      src_grid % geo_p % area = zerd
      src_grid % geo_p % frac = zerd
      src_grid % geo_p % centroid_lat = zerd
      src_grid % geo_p % centroid_lon = zerd

 max_link = src_grid % ne * n_search_scrip / 2
 call allocate_wgts( max_link )

 do cell_add = 1, src_grid % ne
    if (mod(cell_add,1000)==0) write(*,*) 'processing cell',cell_add

    call locate_loc_mesh( cell_add ) 

    !------------------------------------------------
    ! do the formward integration
    !------------------------------------------------
    ncorners = src_grid % nconnect
    call cell_integrate( ncorners,  &
            mesh_src_loc % xgr( mesh_src_loc % in (1:ncorners,1) ), &
            mesh_src_loc % ygr( mesh_src_loc % in (1:ncorners,1) ), &
            mesh_opp_loc, cell_add, list_opp_ele_ref, .true.  )

    !------------------------------------------------
    ! then do the reverse integration
    !------------------------------------------------
    ncorners = opp_grid % nconnect
    do cell_opp = 1, size(list_opp_ele_ref)
       call cell_integrate( ncorners, &
            mesh_opp_loc % xgr( mesh_opp_loc % in (1:ncorners,cell_opp) ), &
            mesh_opp_loc % ygr( mesh_opp_loc % in (1:ncorners,cell_opp) ), &
            mesh_src_loc, list_opp_ele_ref(cell_opp), list_src_ele_ref, .false. )
    enddo

    deallocate(list_src_ele_ref, list_opp_ele_ref)
 enddo

end subroutine grid_integrate


!***********************************************************************

subroutine locate_loc_mesh( cell_add )

!---------------------------------------------------------------------
!
! locate and identify the opposing and source local mesh
!
!---------------------------------------------------------------------
      implicit none

! argumnents

      integer, intent(in) :: cell_add

! locals

      integer, pointer, dimension(:,:) :: src_in, opp_in
      integer, pointer, dimension(:,:) :: reduce_in
      integer n_series, n_opp_nodes, n_opp_ele, ncorners
      integer, allocatable, dimension(:) :: list_opp_ele, list_opp_nodes, list_loc
      type(kdtree2_result), allocatable, dimension(:) :: my_res
      real(kind=4), allocatable, dimension(:) :: my_distances, my_point
      real(kind=8), allocatable, dimension(:) :: &
         lon_src, lat_src, &
         lon_opp, lat_opp
      real(kind=8), pointer, dimension(:) :: &
         xpol_src, ypol_src
      real(kind=8) :: lon0, lat0, srchpt_x, srchpt_y
      integer i, j, k, n, i_loc, j_loc, n_loc
      integer corner
      integer oppcell_add
      logical inpoly

!---------------------------------------------------------------------

! in the following, we assume that src_grid is grid 1 and opp_src is grid 2 for the comments

      ncorners = src_grid % nconnect
      src_in => src_grid % in
      opp_in => opp_grid % in

!--------
! allocate src local mesh
!--------
      call deallocate_mesh( mesh_src_loc )
      mesh_src_loc % nconnect = ncorners
      mesh_src_loc % ne = 1
      mesh_src_loc % nn = ncorners
      allocate( mesh_src_loc % in (ncorners, 1) )
      mesh_src_loc % in(:,1) = (/ (i, i=1,ncorners) /)
      call pointer_mesh( mesh_src_loc, .false. )
      allocate(list_src_ele_ref(1)); list_src_ele_ref(1) = cell_add

!--------
! allocate nodes for source mesh
!--------
      allocate( mesh_src_loc %   xgr(ncorners) )
      allocate( mesh_src_loc %   ygr(ncorners) )
      allocate(  lon_src(ncorners),  lat_src(ncorners) )
      xpol_src => mesh_src_loc % xgr
      ypol_src => mesh_src_loc % ygr

!--------
! defines the center of the projection
!--------

      lon0 = src_grid % lonc(cell_add)
      lat0 = src_grid % latc(cell_add)

!--------
! do a polar projection for the source grid
!--------

      lon_src(:) = src_grid % longr(src_in(:,cell_add))
      lat_src(:) = src_grid % latgr(src_in(:,cell_add))

      n = ncorners
      do i = 1, n
         call ll_to_polars(oned, lat0, lon0, lat_src(i), lon_src(i), xpol_src(i), ypol_src(i))
!         write(*,*) 'corners', i, lon_src(i), lat_src(i),xpol_src(i), ypol_src(i)
      enddo
      deallocate(lon_src, lat_src)

!--------
! find the n-closest points using the KD-tree
!--------

      n_series = n_search_scrip
      allocate(list_opp_ele(n_series),my_res(n_series),my_point(3))

      call ez_lac(my_point, src_grid % lonc(cell_add), src_grid % latc(cell_add), 1)
      call kdtree2_n_nearest( opp_grid % kdtree, my_point, n_series, my_res)
      do i = 1, n_series
         list_opp_ele(i) = my_res(i) % idx
      enddo
      n_opp_ele = n_series
      deallocate(my_res, my_point)

!--------
! move the list of nodes belonging to the nearest elements to a special list
! create a local mesh using the local lists of elements
!--------

      call deallocate_mesh( mesh_opp_loc )
      mesh_opp_loc % nconnect = opp_grid % nconnect
      mesh_opp_loc % ne = n_opp_ele
      allocate( mesh_opp_loc % in (mesh_opp_loc % nconnect, mesh_opp_loc % ne) )
      reduce_in => mesh_opp_loc % in
      reduce_in(:,:) = opp_in(:, list_opp_ele(:) )

!--------
! reconstruct the list of nearest nodes
!--------
      n_loc = n_opp_ele * opp_grid % nconnect ! just a conservative number
      allocate( list_loc(n_loc), list_opp_nodes(n_loc))
      n_loc = 0
      do oppcell_add = 1, n_opp_ele
         do corner = 1, opp_grid % nconnect
            n_loc = n_loc + 1
            list_loc(n_loc) = reduce_in(corner, oppcell_add)
         enddo
!         write(*,*) 'opp cell', oppcell_add, reduce_in(:,oppcell_add)
      enddo

! re-order the list and remove duplicates
      n_opp_nodes = 0
      list_opp_nodes(:) = -1

      do i = 1, n_loc
         i_loc = list_loc(i)
         do j= i+1,n_loc
            j_loc = list_loc(j)
            if (j_loc < i_loc) then
               list_loc(i) = j_loc
               list_loc(j) = i_loc
               i_loc = j_loc
            endif
         enddo
         if ( n_opp_nodes==0 ) then
            n_opp_nodes = n_opp_nodes + 1
            list_opp_nodes(n_opp_nodes) = i_loc
         else
            if ( list_loc(i) > list_loc(i-1) ) then
               n_opp_nodes = n_opp_nodes + 1
               list_opp_nodes(n_opp_nodes) = i_loc
            endif
         endif
      enddo

! make reduce_in a pointer to the local list
      do i = 1, n_opp_nodes
         where( reduce_in(:,:) == list_opp_nodes(i)) reduce_in(:,:) = i
      enddo

      deallocate(list_loc)
      mesh_opp_loc % nn = n_opp_nodes
      allocate( mesh_opp_loc % xgr(n_opp_nodes), &
                mesh_opp_loc % ygr(n_opp_nodes) )
      allocate( lon_opp(n_opp_nodes), lat_opp(n_opp_nodes) )
      lon_opp(1:n_opp_nodes) = opp_grid % longr(list_opp_nodes(1:n_opp_nodes))
      lat_opp(1:n_opp_nodes) = opp_grid % latgr(list_opp_nodes(1:n_opp_nodes))
      do i = 1, n_opp_nodes
         call ll_to_polars(oned, lat0, lon0, lat_opp(i), lon_opp(i), mesh_opp_loc % xgr(i), mesh_opp_loc % ygr(i))
!         write(*,*) 'opp_nodes', i, mesh_opp_loc % xgr(i), mesh_opp_loc % ygr(i)
      enddo
      call pointer_mesh( mesh_opp_loc, .false. )
      deallocate( lon_opp, lat_opp )

!----------------
! pass list of opp cell to global variable
!----------------

      if (allocated(list_opp_ele_ref)) deallocate(list_opp_ele_ref)
      allocate(list_opp_ele_ref(n_opp_ele))
      list_opp_ele_ref(:) = list_opp_ele(:)
      deallocate(list_opp_ele, list_opp_nodes)

end subroutine locate_loc_mesh


!***********************************************************************

      subroutine cell_integrate(ncorners, xsrc, ysrc, mesh_opp, cell_add, list_opp, order_weight)

!-----------------------------------------------------------------------
!
!     Integrate around cell while finding intersecting with opposite
!     grid cells and finding segments of cell boundary lying in cells
!     of opposite grid
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     Input variables
!
!-----------------------------------------------------------------------
      
      integer, intent(in) :: ncorners
      real(kind=8), dimension(:), intent(in) :: &
              xsrc, ysrc
      type(type_mesh) :: mesh_opp
      integer :: cell_add   ! cell to be processed
      integer, dimension(:), intent(in) :: list_opp
      logical :: order_weight


!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer, parameter ::&
     &     max_subseg = 100     ! max number of subsegments per segment
                                ! to prevent infinite loop


      integer :: &
     &     i, inext,       &     !
     &     j, jnext,       &     ! generic counters
     &     ic, k, ns,      &     !
     &     n, next_n,      &     !
     &     nwgt, it,       &     !
     &     oppcell_add,    &     ! Cell from opposite grid we are intersecting
     &     min_add,        &     ! addresses for restricting search of
     &     max_add,        &     !   destination grid
     &     corner,         &     ! corner of cell that segment starts from
     &     next_corn,      &     ! corner of cell that segment ends on
     &     nseg,           &     ! number of segments to use to represent 
                                ! edges near the pole                   
     &     num_subseg,     &     ! number of subsegments
     &     bedgeid1,       &     !
     &     bedgeid2,       &     ! ID of edge that a point is on
     &     bedgeid3,       &     !
     &     intedge,        &     ! ID of intersected edge
     &     last_add,       &     ! Address of last cell we were in
     &     next_add,       &     ! Address of next cell we will go into
     &     adj_add              ! Address of cell adjacent to current one

      logical :: &
     &     lcoinc,         &     ! Are segments coincident?
     &     lrevers,        &     ! Are we integrating segment in reverse?
     &     lboundary1,     &
     &     lboundary2,     &     ! Is point is on cell boundary?
     &     lboundary3,     &     
     &     last_lboundary, &     ! Is last point is on cell bdry?
     &     lthresh,        &     ! Has segment crossed threshold?
     &     srch_success,   &     ! Was search for segment start successful?
     &     intrsct_success,&     ! Was intersection of segment with opposite
                                 ! grid successful?
     &     inpoly,         &     ! Is point is in polygon
     &     inface,         &     ! is intersecting a face
     &     last_endpt_inpoly,&   ! Was end point of last segment in cell
     &     lstuck,           &   ! Is the walk stuck inside a cell
     &     seg_outside,      &   ! Is segment completely outside the grid
     &     bndedge,          &   ! Is segment on the boundary of the grid
     &     search                ! Do we have to search to locate point in grid

      real(kind=8) ::   &
     &     intrsct_x,       &   ! x of next intersection point
     &     intrsct_y,       &   ! y of next intersection point
     &     begx, begy,      &   ! start point of current sub seg
     &     endx, endy,      &   ! endpoint of current seg (chg to endseg?)
     &     norm_factor          ! factor for normalizing wts

      real(kind=8), dimension(2) :: &
     &     begseg               ! begin lat/lon for full segment

! FD no gradient      real(kind=8), dimension(6) :: 
      real(kind=8), dimension(3) :: &
     &     weights              ! local wgt array

      real(kind=8) :: &
     &     vec1_x, vec1_y,     & ! vectors, products
     &     vec2_x, vec2_y,     & ! used in grid search
     &     vec1_len, dp,       &
     &     midx, midy,         & ! Midpoint of segment
     &     tmpx, tmpy,         &
     &     dist2,              & ! Square of distance b/w two points
     &     fullseg_len2,       & ! Square of full segment length
     &     partseg_len2,       & ! Square of length of segment integrated so far
     &     fullseg_dx,         & ! x diff of full segment endpoints
     &     fullseg_dy,         & ! y diff of full segment endpoints
     &     cross_product

      integer :: &
     &     tmpnseg

      integer, allocatable, dimension(:) :: list_loc
      real(kind=8), dimension( mesh_opp % nconnect ) :: srch_x_loc, srch_y_loc
      integer n_opp_ele
      integer n_opp_nodes
      real(kind=8), pointer, dimension(:) :: &
         xpol_opp, ypol_opp
      

      xpol_opp => mesh_opp % xgr
      ypol_opp => mesh_opp % ygr
      n_opp_ele = mesh_opp % ne
!--------
! back to the integration
!--------

      allocate( list_loc( mesh_opp % nconnect ) )

! FD debug
!write(*,*) 'source element',cell_add
!       do i = 1 , ncorners
!         write(*,*) xsrc(i),ysrc(i)
!         enddo
      corner_loop: do corner = 1, ncorners
         ! FD debug
!         write(*,*) 'corner',corner
         next_corn = mod(corner,ncorners) + 1

         !***
         !*** define endpoints of the current segment
         !***

         begx = xsrc(corner)
         begy = ysrc(corner)
         endx = xsrc(next_corn)
         endy = ysrc(next_corn)

         begseg(1) = begx
         begseg(2) = begy

         fullseg_dx = endx-begx
         fullseg_dy = endy-begy
         fullseg_len2 = fullseg_dx * fullseg_dx + fullseg_dy * fullseg_dy

         if ( fullseg_len2 < tiny ) cycle ! bail out if the segment is degenerate (pole?). the area calculation gives zero anyway

         partseg_len2 = 0.0 ! partial segment length (useful when nseg>1)


         !***
         !*** Is this an edge on the boundary of the grid or 
         !*** on the boundary of the active cells
         !*** 

         !***
         !*** integrate along this segment, detecting intersections 
         !*** and computing the line integral for each sub-segment
         !***

         nseg = 1

         last_add = 0
         last_lboundary = .false.
         last_endpt_inpoly = .false.
         next_add = 0
         search = .true.
         ns = 1
         inface = .false.

         tmpnseg = 0

         do while (begx /= endx .or. begy /= endy)
            
         ! FD debug
!         write(*,*) 'subseg',ns,nseg,begx,begy,endx,endy
            
            num_subseg = 0
               
               !*** 
               !*** If we integrated to the end or just past it (due to 
               !*** numerical errors), we are done with this segment
               !***

               if (partseg_len2 .ge. fullseg_len2) exit



               !******************************************************
               !*** Try to find which cell of grid 2 this segment is 
               !*** starting in and where it is exiting this cell
               !******************************************************
         
               vec1_x = endx - begx
               vec1_y = endy - begy
               vec1_len = sqrt( vec1_x * vec1_x + vec1_y * vec1_y)
               vec1_x = vec1_x / vec1_len
               vec1_y = vec1_y / vec1_len
               
               oppcell_add = 0
               intrsct_success = .false.
               lstuck = .false.
               lboundary1 = .false.
               lboundary2 = .false.
               lboundary3 = .false.
               

                  !*************************************************
                  !*** Find out which cell the segment starts in
                  !*************************************************

                  ! FD debug
!                  write(*,*) 'beginning of search of sub-segment'
!                  write(*,*) 'last_add',last_add

                  srch_success = .false.
                  if (search) then

                     !***
                     !*** Offset the start point in ever increasing 
                     !*** amounts until we are able to reliably locate 
                     !*** the point in a cell of opp_grid. Inability to locate 
                     !*** the point causes the offset amount to increase

! find whether the point under consideration is inside a polygon of the opposite grid
                        oppcell_add = 0
                        inpoly = .false.

                        do while ( .not.inpoly .and. oppcell_add < n_opp_ele )

                           oppcell_add = oppcell_add + 1

                           list_loc(:) = mesh_opp % in(:,oppcell_add)
             
                           srch_x_loc(:) = xpol_opp(list_loc(:))
                           srch_y_loc(:) = ypol_opp(list_loc(:))

                           call ptinpoly(begx, begy, mesh_opp % nconnect, &
     &                        srch_x_loc, srch_y_loc, inpoly, lboundary1, bedgeid1)

                        enddo ! end search of opposing cell

! FD debug
!write(*,*) 'find seg starts in',inpoly,oppcell_add,lboundary1,bedgeid1,begx,begy

                        ! if the beginning of the segment does not fall inside an opposing cell
                        ! find the closest interesection with a face of the opposing mesh
                        ! this will give the start of the new sub-segment
                        ! if no intersection is found, then cycle to the next corner

                        if (.not.inpoly) then

                           call find_closest_face_start_seg(mesh_opp, begx, begy, endx, endy, &
                                   oppcell_add, inface, .false.)

                           if (.not.inface) cycle corner_loop ! go to the next corner

                           srch_success = .true.
                           inpoly = .true.

                           list_loc(:) = mesh_opp % in(:,oppcell_add)
             
                           srch_x_loc(:) = xpol_opp(list_loc(:))
                           srch_y_loc(:) = ypol_opp(list_loc(:))

                        else
                           ! need to make sure that we have not prcessed that grid cell already
                           if (oppcell_add .ne. last_add) then
                              srch_success = .true.
                           endif
                        endif
                        
                  else ! no search: just use past info

                     if (last_endpt_inpoly) then

                        !*** We know the grid cell at the end of the last 
                        !*** segment (which is the beginning of this 
                        !*** segment)

                        oppcell_add = last_add
                        lboundary1 = last_lboundary

                     else if (next_add .ne. 0) then

                        !*** We know the edge of the opp_grid cell that the
                        !*** last segment intersected, so we move into
                        !*** the adjacent cell

                        oppcell_add = next_add
                        lboundary1 = .true.

                        list_loc(:) = mesh_opp % in(:,oppcell_add)
             
                        srch_x_loc(:) = xpol_opp(list_loc(:))
                        srch_y_loc(:) = ypol_opp(list_loc(:))
                        
                     else if (next_add == 0) then

                        call find_closest_face_start_seg(mesh_opp, begx, begy, endx, endy, oppcell_add, &
                                inface, .true.) ! search for next face excluding the nearest

                        if (.not.inface) cycle corner_loop ! go to the next corner

                        srch_success = .true.
                        inpoly = .true.

                        list_loc(:) = mesh_opp % in(:,oppcell_add)
             
                        srch_x_loc(:) = xpol_opp(list_loc(:))
                        srch_y_loc(:) = ypol_opp(list_loc(:))

                     endif
                     
                     srch_success = .true.
                     
                  endif
! FD debug
! write(*,*) 'search',search,corner,oppcell_add,begx,begy
! write(*,*) 'found opp cell',srch_success
                  
                  
               if (srch_success) then 

                     !*****************************************************
                     !*** Find where the segment exits this cell, if at all
                     !*****************************************************

                     !***
                     !*** First see if the segment end is in the same cell
                     !***

! FD debug
!write(*,*) 'checking element',oppcell_add
!       do i = 1, mesh_opp % nconnect
!         write(*,*) srch_x_loc(i),srch_y_loc(i)
!         enddo
!         write(*,*)
!write(*,*) 'search end of segment'
                     call ptinpoly(endx, endy, mesh_opp % nconnect, &
     &                    srch_x_loc, srch_y_loc, &
                          inpoly, lboundary2, bedgeid2)
! FD debug
!write(*,*) 'found end seg',inpoly,oppcell_add,lboundary2,bedgeid2

                     if (inpoly) then
                        intrsct_x = endx
                        intrsct_y = endy
                        intrsct_success = .true.                  
                        search = .false.
                        next_add = 0
                        last_add = oppcell_add        ! for next subseg
                        last_lboundary = lboundary2
                        last_endpt_inpoly = .true.
                        
                        if (lboundary1 .and. lboundary2) then
                  
                           !*** This is a edge on the boundary of the 
                           !*** active mesh and both of its endpoints 
                           !*** are on the boundary of the containing 
                           !*** cell. Check if the segment is also 
                           !*** on the boundary

                           midx = ( begx + endx ) / 2.0_8
                           midy = ( begy + endy ) / 2.0_8
                           
                           call ptinpoly(midx, midy, mesh_opp % nconnect, &
     &                    srch_x_loc, srch_y_loc, inpoly, lboundary3, bedgeid3) 
                           
                           if (inpoly .and. lboundary3) then
                              lcoinc = .true.
                              intedge = bedgeid3
                           endif
                           
                        else
                           lcoinc = .false.
                        endif
                        
                     else ! the end-segment is not in the same polygone

                      
                        !***
                        !*** Do an intersection to find out where the 
                        !*** segment exits the cell
                        !***

                        call intersection( begx, begy, endx, endy, &
                             begseg, bedgeid1, &
     &                       oppcell_add, mesh_opp % nconnect, &
     &                       srch_x_loc, srch_y_loc, &
     &                       intrsct_x, intrsct_y, intedge, lcoinc)
! FD debug
!write(*,*) 'intersect',corner,oppcell_add,intrsct_x,intrsct_y, intedge
                      
                        if (intedge /= 0) then
                           intrsct_success = .true.
                           last_add = oppcell_add     ! for next subseg
                           last_endpt_inpoly = .false.
                           last_lboundary = .true.

                           call find_adj_cell(mesh_opp, oppcell_add, intedge, next_add)
! FD debug
!write(*,*) 'next opp cell',next_add
                              search = .false.
                        else
                           intrsct_x = begx
                           intrsct_y = begy
                           search = .false.
                           next_add = 0

                           ! maybe an opposing corner, need to reset the search
                           call find_other_poly(begx, begy, endx, endy, last_add, oppcell_add, mesh_opp, next_add, inface)
                           oppcell_add = 0 ! reset before looping back

                        endif ! test on intersecing sucess

                     endif ! condition on end seg being inside the opposing cell
                    
                     if (lcoinc) then
                      
                        !***
                        !*** Segment is coincident with edge of other grid
                        !*** which means it could belong to one of 2 cells
                        !*** Choose the cell such that edge that is 
                        !*** coincident with the segment is in the same
                        !*** dir as the segment

                           i = intedge
                           inext  = mod(i,mesh_opp % nconnect)+1
                           vec2_x = srch_x_loc(inext) &
     &                            - srch_x_loc(i)
                           vec2_y = srch_y_loc(inext) &
     &                            - srch_y_loc(i)

                        dp = vec1_x * vec2_x + vec1_y * vec2_y
                        
                     endif

               else ! no srch_success

                     oppcell_add = 0
                     last_add = 0

                     cycle corner_loop ! go to the next corner
                     
               endif      ! if (srch_success) then ... else ....



               !**********************************************************
               !*** Compute the line integrals for this subsegment
               !**********************************************************

               cross_product = zerd
            
                  if (oppcell_add /= 0) then
                     call vector_area(cross_product, &
     &                    zerd, zerd, begx, begy, intrsct_x, intrsct_y)
                  endif

                  weights(1) = cross_product
                  weights(2:3) = zerd
! FD debug
! write(*,'(3i6,5(1x,e23.16))') corner,cell_add,oppcell_add,weights(1),begx, begy, intrsct_x, intrsct_y
             
               !***
               !*** store the appropriate addresses and weights. 
               !*** also add contributions to cell areas and centroids.
               !***

                  if (oppcell_add /= 0 .and. intrsct_success) then

! FD debug
! write(*,*) 'before storage',cell_add,list_opp(oppcell_add)
                     if (order_weight) then
                        call store_link_cnsrv(cell_add, list_opp(oppcell_add), weights)
                     else
! write(*,*) 'reverse storage',list_opp(oppcell_add),cell_add
                        call store_link_cnsrv(list_opp(oppcell_add), cell_add, weights)
                     endif

                  endif

               tmpnseg = tmpnseg + 1
               

               !***
               !*** reset beglat and beglon for next subsegment.
               !***

               begx = intrsct_x
               begy = intrsct_y
               

               !***
               !*** How far have we come from the start of the segment
               !***
             
               vec2_x = intrsct_x - begseg(1)
               vec2_y = intrsct_y - begseg(2)
               
               partseg_len2 = vec2_x * vec2_x + vec2_y * vec2_y


               !***
               !*** prevent infinite loops if integration gets stuck
               !*** near cell or threshold boundary
               !***

               num_subseg = num_subseg + 1
               if (num_subseg > max_subseg) then
                  print *, &
     &               'integration stalled: num_subseg exceeded limit'
                  print *, 'Cell ',cell_add
                  print *, 'Edge ',corner
                  print *, 'Grid ',1
                  dist2 = (endx-begseg(1))*(endx-begseg(1)) &
                        + (endy-begseg(2))*(endy-begseg(2))
                  print *, 'Fraction of segment left ', vec1_len/sqrt(dist2)
                  stop ! FD debug
                  exit       ! Give up
               endif


         end do              ! do while (beglat /= endlat ....


         !***
         !*** end of segment
         !***

      end do corner_loop                   ! do corner=....

      deallocate(list_loc)

      end subroutine cell_integrate




!***********************************************************************

      subroutine intersection(begx, begy, endx, endy, &
     &     begseg, begedge, &
     &     cell_id, ncorners, cell_corner_x, cell_corner_y, &
     &     intrsct_x, intrsct_y, intedge, lcoinc)

!-----------------------------------------------------------------------
!
!     this routine finds the intersection of a line segment given by 
!     beglon, endlon, etc. with a cell from another grid
!     A coincidence flag is returned if the segment is entirely 
!     coincident with an edge of the opposite.  
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      real(kind=8), intent(in) :: &
     &     begx, begy,  &   ! beginning lat/lon endpoints for segment
     &     endx, endy       ! ending    lat/lon endpoints for segment

      real(kind=8), dimension(2), intent(in) :: &
     &     begseg               ! begin lat/lon of full segment

      integer, intent(in) :: &
     &     begedge              ! edge that beginning point is on (can be 0)

      integer, intent(in) :: &
     &     cell_id              ! cell to intersect with
      
      integer, intent(in) :: &
     &     ncorners             ! number of corners of cell

      real(kind=8), dimension(ncorners), intent(in) :: &
     &     cell_corner_x, &   ! coordinates of cell corners
     &     cell_corner_y


!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      real(kind=8), intent(out) :: &
     &     intrsct_x, &
     &     intrsct_y          ! lat/lon coords of intersection

      integer, intent(out) :: &
     &     intedge              ! edge that is intersected

      logical, intent(out) :: &
     &     lcoinc               ! True if segment is coincident with 
                                ! a cell edge

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: &
     &     n, next_n

      logical :: found

      real(kind=8) :: &
     &     x1, x2,     &    ! local x variables for segment
     &     y1, y2,      &   ! local y variables for segment
     &     grdx1, grdx2, &  ! local x variables for grid cell
     &     grdy1, grdy2,  & ! local y variables for grid cell
     &     vec1_x, vec1_y,& 
     &     vec2_x, vec2_y,& !
     &     vec3_x, vec3_y,& ! vectors and vector products used
     &     cross_product,     & ! during grid search
     &     dot_product,       & !
     &     lensqr1, lensqr2,  & !
     &     lensqr3,           & !
     &     s1, s2, determ,    & 
     &     mat1, mat2,        & ! variables used for linear solve to
     &     mat3, mat4,        & ! find intersection
     &     rhs1, rhs2,        & !
     &     denom,             & 
     &     begsegloc(2),      & ! local copy of full segment start
     &     dist2,             & ! distance from start pt to intersection pt
     &     maxdist2,          & ! max dist from start pt to any intersection pt
     &     tmpx, tmpy

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      lcoinc = .false.
      intedge = 0
      intrsct_x = 0.d0
      intrsct_y = 0.d0

      x1 = begx
      y1 = begy
      x2 = endx
      y2 = endy

      ! No edge is allowed to span more than pi radians
      ! Accordingly transform one or the other end point

      s1 = zerd

!-----------------------------------------------------------------------
!
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

      begsegloc(1) = begseg(1)
      begsegloc(2) = begseg(2)

      intrsct_loop: do n=1,ncorners
         next_n = mod(n,ncorners) + 1
         
         grdx1 = cell_corner_x(n)
         grdy1 = cell_corner_y(n)
         grdx2 = cell_corner_x(next_n)
         grdy2 = cell_corner_y(next_n)

         lensqr2 = (grdx1-grdx2)*(grdx1-grdx2) + &
     &             (grdy1-grdy2)*(grdy1-grdy2)

         if (lensqr2 .le. tiny*tiny) cycle       ! degenerate edge

        !***
        !*** set up linear system to solve for intersection
        !***
! FD if M1=(x1,y1), M2=(x2,y2), G1=(grdx1,grdy1), G2=(grdx2,grdy2) and O the grid center
! FD the matrix problem to solve for is
! FD 
! FD                      OM1 + a M1M2 = OG1 + b G1G2
! FD which reduces to
! FD                    a M1M2 - b G1G2 = M1G1

        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident.  coincidences were detected 
        !***   above so do nothing.

        if (abs(determ) > tiny*tiny) then

          !*** if the determinant is non-zero, solve for the linear 
          !***   parameters s for the intersection point on each line 
          !***   segment.
          !*** if 0<s1,s2<1 then the segment intersects with this side.
          !***   return the point of intersection (adding a small
          !***   number so the intersection is off the grid line).
          !***

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zerd .and. s2 <= oned .and. &
     &        s1 >  zerd .and. s1 <= oned) then

             !***
             !*** recompute intersection based on full segment
             !*** so intersections are consistent for both sweeps
             !***

             mat1 = x2 - begsegloc(1)
             mat3 = y2 - begsegloc(2)
             rhs1 = grdx1 - begsegloc(1)
             rhs2 = grdy1 - begsegloc(2)

             determ = mat1*mat4 - mat2*mat3

             !***
             !*** sometimes due to roundoff, the previous 
             !*** determinant is non-zero, but the lines
             !*** are actually coincident.  if this is the
             !*** case, skip the rest.
             !***

             if (determ /= zerd) then
                s1 = (rhs1*mat4 - mat2*rhs2)/determ
                s2 = (mat1*rhs2 - rhs1*mat3)/determ

                intrsct_x = begsegloc(1) + mat1*s1
                intrsct_y = begsegloc(2) + mat3*s1

                !***
                !*** Make sure the intersection point is not within
                !*** tolerance of the starting point
                !***

                vec1_x = intrsct_x - begx
                vec1_y = intrsct_y - begy

                dist2 = vec1_x * vec1_x + vec1_y * vec1_y
                   
                if (dist2 > tiny) then
                   intedge = n
                   exit intrsct_loop
                endif

             else
                print *, 'DEBUG: zero determ'
                stop
             endif
             
          endif

       else

          !***
          !*** Coincident lines or parallel lines
          !*** 

          cross_product = mat2*rhs2 - mat4*rhs1

          !***
          !*** If area of triangle formed by endx,endy and
          !*** the gridline is negligible then the lines are coincident
          !***

          if (abs(cross_product) < tiny) then

             dot_product = mat1*(-mat2) + mat3*(-mat4)

             lensqr1 = mat1*mat1 + mat3*mat3 ! length sqrd of input segment

             if (dot_product < zerd) then

                !***
                !*** Segments oriented in the same direction
                !***
             

                tmpx = grdx2
                tmpy = grdy2
                grdx2 = grdx1
                grdy2 = grdy1
                grdx1 = tmpx
                grdy1 = tmpy

             endif


             vec2_x = grdx1 - x1
             vec2_y = grdy1 - y1

             lensqr2 = vec2_x * vec2_x + vec2_y * vec2_y

             if (vec2_x*mat1 + vec2_y*mat3 < 0) then
                lensqr2 = -lensqr2
             endif

             vec3_x = grdx2 - x1
             vec3_y = grdy2 - y1

             lensqr3 = (vec3_x * vec3_x + vec3_y * vec3_y)

             if (vec3_x*mat1 + vec3_y*mat3 < 0) then
                lensqr3 = -lensqr3
             endif
             
             found = .false.

             if (lensqr2 > 0) then
                if (lensqr2 <= lensqr1) then
                   intrsct_x = grdx1
                   intrsct_y = grdy1
                   found = .true.
                endif
             else
                if (lensqr3 > 0) then
                   if (lensqr3 > lensqr1) then
                      intrsct_x = x2
                      intrsct_y = y2
                      found = .true.
                   else
                      intrsct_x = grdx2
                      intrsct_y = grdy2
                      found = .true.
                   endif
                endif
             endif

             if (found) then

                dist2 = (intrsct_x-begx)*(intrsct_x-begx)+ &
     &                  (intrsct_y-begy)*(intrsct_y-begy)

               !*** Coincidence intersection always wins

                intedge = n
                lcoinc = .true.

                exit intrsct_loop
             endif
                
          endif

       endif

      end do intrsct_loop

!-----------------------------------------------------------------------

      end subroutine intersection

!***********************************************************************

      subroutine find_closest_face_start_seg(mesh_opp, begx, begy, endx, endy, cell_id, inface, exclude_zero_dist)
      implicit none

!-----------------------------------------------------------------------
!
!     this routine finds the intersection of a line segment given by 
!     begx, begy, endx, endy, with a face from another grid
!     test all faces, keep only the closest
!     or bail out with inface = .false.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in): 
!
!-----------------------------------------------------------------------

      type(type_mesh) :: mesh_opp

      real(kind=8), intent(in) :: &
     &     endx, endy       ! ending    lat/lon endpoints for segment
      logical, intent(in) :: &
     &     exclude_zero_dist

!-----------------------------------------------------------------------
!
!     intent(out): 
!
!-----------------------------------------------------------------------

      real(kind=8), intent(inout) :: &
     &     begx, begy   ! beginning lat/lon endpoints for segment

      integer, intent(out) :: &
     &     cell_id              ! cell to intersect with

      logical, intent(out) :: inface

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: &
     &     n1, n2, face, intedge

      logical :: found

      real(kind=8) :: &
     &     x1, x2,     &    ! local x variables for segment
     &     y1, y2,      &   ! local y variables for segment
     &     grdx1, grdx2, &  ! local x variables for grid cell
     &     grdy1, grdy2,  & ! local y variables for grid cell
     &     vec1_x, vec1_y,& 
     &     vec2_x, vec2_y,& !
     &     vec3_x, vec3_y,& ! vectors and vector products used
     &     cross_product,     & ! during grid search
     &     intrsct_x, intrsct_y, & ! intersection point
     &     lenseg,            & !
     &     lensqr1, lensqr2,  & !
     &     lensqr3,           & !
     &     s1, s2, determ,    & 
     &     mat1, mat2,        & ! variables used for linear solve to
     &     mat3, mat4,        & ! find intersection
     &     rhs1, rhs2,        & !
     &     denom,             & 
     &     begsegloc(2),      & ! local copy of full segment start
     &     dist2,             & ! distance from start pt to intersection pt
     &     mindist2,          & ! min dist from start pt to any intersection pt
     &     tmpx, tmpy

!-----------------------------------------------------------------------
!
!     initialize defaults, flags, etc.
!
!-----------------------------------------------------------------------

      inface = .false.
      cell_id = 0

      x1 = begx
      y1 = begy
      x2 = endx
      y2 = endy

      lenseg = (x1-x2) * (x1-x2) + (y1-y2) * (y1-y2)
      mindist2 = lenseg

!-----------------------------------------------------------------------
!
!     loop over sides of the cell to find intersection with side
!     must check all sides for coincidences or intersections
!
!-----------------------------------------------------------------------

      intrsct_loop: do face = 1, mesh_opp % nf
         n1 = mesh_opp % pseg(1,face)
         n2 = mesh_opp % pseg(2,face)
         
         grdx1 = mesh_opp % xgr(n1)
         grdy1 = mesh_opp % ygr(n1)
         grdx2 = mesh_opp % xgr(n2)
         grdy2 = mesh_opp % ygr(n2)

         lensqr2 = (grdx1-grdx2)*(grdx1-grdx2) + &
     &             (grdy1-grdy2)*(grdy1-grdy2)

         if (lensqr2 .le. tiny*tiny) cycle       ! degenerate edge

        !***
        !*** set up linear system to solve for intersection
        !***
! FD if M1=(x1,y1), M2=(x2,y2), G1=(grdx1,grdy1), G2=(grdx2,grdy2) and O the grid center
! FD the matrix problem to solve for is
! FD 
! FD                      OM1 + a M1M2 = OG1 + b G1G2
! FD which reduces to
! FD                    a M1M2 - b G1G2 = M1G1

        mat1 = x2 - x1
        mat2 = grdx1 - grdx2
        mat3 = y2 - y1
        mat4 = grdy1 - grdy2
        rhs1 = grdx1 - x1
        rhs2 = grdy1 - y1

        determ = mat1*mat4 - mat2*mat3

        !***
        !*** if the determinant is zero, the segments are either 
        !***   parallel or coincident.  coincidences were detected 
        !***   above so do nothing.

        if (abs(determ) > tiny*tiny) then

          !*** if the determinant is non-zero, solve for the linear 
          !***   parameters s for the intersection point on each line 
          !***   segment.
          !*** if 0<s1,s2<1 then the segment intersects with this side.
          !***   return the point of intersection (adding a small
          !***   number so the intersection is off the grid line).
          !***

          s1 = (rhs1*mat4 - mat2*rhs2)/determ
          s2 = (mat1*rhs2 - rhs1*mat3)/determ

          if (s2 >= zerd .and. s2 <= oned .and. &
     &        s1 >  zerd .and. s1 <= oned) then

                intrsct_x = begx + mat1*s1
                intrsct_y = begy + mat3*s1

                !***
                !*** Make sure the intersection point is not within
                !*** tolerance of the starting point
                !***

                vec1_x = intrsct_x - begx
                vec1_y = intrsct_y - begy

                dist2 = vec1_x * vec1_x + vec1_y * vec1_y
!                write(*,*) 'face',face,s1,s2,intrsct_x,intrsct_y,dist2 ! FD debug

                if (exclude_zero_dist .and. dist2 < tiny ) then
                   cycle intrsct_loop
                endif
                   
                inface = .true.

                if (dist2 < mindist2 ) then
                   intedge = face
                   mindist2 = dist2
                   tmpx = intrsct_x
                   tmpy = intrsct_y
                endif

          endif

        endif

      end do intrsct_loop

      ! retrieve the opposing id number
      !----------------------------------
      ! logically there should be only one cell connected to this face, as if there were another one, it would even closer

      if (inface) then

         cell_id = mesh_opp % tseg(1,intedge) 
         begx = tmpx
         begy = tmpy
!                write(*,*) 'final',intedge, mindist2, cell_id ! FD debug

      endif

!-----------------------------------------------------------------------

      end subroutine find_closest_face_start_seg

!***********************************************************************

!*************************************************************************

      subroutine vector_area(vecp, x0, y0, x1, y1, x2, y2)

      implicit none
!-----------------------------------------------------------------------
!
!     this routine computes the line integral of the flux function 
!     that results in the interpolation weights.  the line is defined
!     by the input lat/lon of the endpoints. Integration is w.r.t. lon
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     intent(in):
!
!-----------------------------------------------------------------------

      real(kind=8), intent(in) :: &
           x0, y0, & ! coord of first point
           x1, y1, & ! second
           x2, y2    ! third

!-----------------------------------------------------------------------
!
!     intent(out):
!
!-----------------------------------------------------------------------

      real(kind=8), intent(out) :: &
     &     vecp   ! signed area defined by the cross product of 01 and 02

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      real(kind=8) dxa, dya, dxb, dyb

!-----------------------------------------------------------------------
!
!     weights for the general case based on a trapezoidal approx to
!     the integrals.
!
!-----------------------------------------------------------------------

      dxa = x1 - x0
      dya = y1 - y0
      dxb = x2 - x0
      dyb = y2 - y0
      vecp = dyb*dxa - dya*dxb
      vecp = vecp * 0.5_8
      
      end subroutine vector_area

!*************************************************************************

!**********************************************************************

      subroutine ptinpoly(ptx, pty, ncorners, cell_corner_x, &
     &     cell_corner_y, inpoly, lboundary, edgeid)

!----------------------------------------------------------------------
!
!     Check if point is in (convex) polygonal cell 
!
!----------------------------------------------------------------------

      implicit none

!----------------------------------------------------------------------
!
!     Input arguments
!
!----------------------------------------------------------------------

      real(kind=8), intent(in) :: &
     &     ptx, pty             ! Point to check

      integer, intent(in) :: &
     &     ncorners             ! Number of polygon corners

      real(kind=8), dimension(ncorners), intent(in) :: &
     &     cell_corner_x,   &   ! Coordinates of cell corners
     &     cell_corner_y        ! Could be x-y or lat-lon or ...

!----------------------------------------------------------------------
!
!     Output arguments
!
!----------------------------------------------------------------------

      logical, intent(out) :: &
     &     inpoly               ! Is point in the polygon?

      logical, intent(out) :: &
     &     lboundary            ! Is point on the boundary of the polygon?

      integer, intent(out) :: &
     &     edgeid               ! if point is on boundary, which local
                                ! edge is it on? (0 otherwise)

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer :: n, next_n

      real(kind=8) :: x1, y1, x2, y2, vec1_x, vec1_y, vec2_x, vec2_y, &
     &     cross_product

      !***
      !*** here we take the cross product of the vector making 
      !*** up each cell side with the vector formed by the vertex
      !*** and search point.  if all the cross products are 
      !*** positive, the point is contained in the cell.
      !***

      inpoly = .false.
      lboundary = .false.
      edgeid = 0
        

         !*** Checking in latlon space
         !*** If the grid cell coordinates spans more than pi radians
         !*** transform the coordinates so that they don't

         do n = 1, ncorners
            next_n = MOD(n,ncorners) + 1
            
            x1 = cell_corner_x(n)
            y1 = cell_corner_y(n)
            x2 = cell_corner_x(next_n)
            y2 = cell_corner_y(next_n)
            
            vec1_x = x2 - x1
            vec1_y = y2 - y1
            vec2_x = ptx - x1
            vec2_y = pty - y1
                        
            cross_product = vec1_x * vec2_y - vec2_x * vec1_y
            
            !***
            !***   if the cross product for a side is zero, the point 
            !***   lies exactly on the side or the side is degenerate
            !***   (zero length).  if degenerate, set the cross 
            !***   product to a positive number.  
            !***

            if (abs(cross_product) < tiny) then
               lboundary = .true.
               edgeid = n
               cross_product = zerd
            endif

            !***
            !*** if cross product is less than zero, this cell
            !*** doesn't work
            !***
            !*** Should we say "if (cp < zero .and. abs(cp) > tiny)" ?

            if (cross_product < zerd) then
               inpoly = .false.
               lboundary = .false.
               return
            endif
            
         end do

      !***
      !*** if cross products all positive, we found the location
      !***

      inpoly = .true.

!----------------------------------------------------------------------

      end subroutine ptinpoly
! 
! !**********************************************************************



!----------------------------------------------------------------------
!
!     Find cell adjacent to edge (edge_id) of given cell (cell_add)
!
!----------------------------------------------------------------------


      subroutine find_adj_cell(grid_loc, cell_add, edge_id, adj_add) 

!----------------------------------------------------------------------
!
!     Input variables
!
!----------------------------------------------------------------------

      integer, intent(in) :: &
     &     cell_add,         &    ! cell whose edge we are checking
     &     edge_id                ! index of edge that we are check

      TYPE (type_mesh), intent(in), target  :: grid_loc

!----------------------------------------------------------------------
!
!     Output variables
!
!----------------------------------------------------------------------

      integer, intent(out) :: adj_add

!----------------------------------------------------------------------
!
!     Local variables
!
!----------------------------------------------------------------------

      integer :: cell1, cell2, face

      INTEGER, dimension(:,:), pointer :: tseg_loc, itm_loc

      tseg_loc => grid_loc % tseg
      itm_loc  => grid_loc % itm

!$OMP THREADPRIVATE(cell1,cell2,adj_add,face)

      adj_add = 0

      face = itm_loc(edge_id, cell_add)
      cell1 = tseg_loc(1, face)
      cell2 = tseg_loc(2, face)

      if ( cell1 == cell_add ) adj_add = cell2
      if ( cell2 == cell_add ) adj_add = cell1

      end subroutine find_adj_cell


!----------------------------------------------------------------------
!
!     !!! Panick: we have encountered a corner !!!
!     Find the cell that includes fully or partially the given segment
!
!----------------------------------------------------------------------

      subroutine find_other_poly(begx, begy, endx, endy, last_add, oppo_add, mesh_opp, next_add, found)
      ! arguments
      real(kind=8), intent(in) :: begx, begy, endx, endy
      integer, intent(in) :: last_add, oppo_add
      integer, intent(out) :: next_add
      type(type_mesh), intent(in) :: mesh_opp
      logical, intent(out) :: found
      ! locals
      integer     ncorners             ! Number of corners in cell
      real(kind=8), allocatable, dimension(:) :: &
     &     cell_corner_x,  &  ! xand y values of cell corners
     &     cell_corner_y
      integer, allocatable, dimension(:) :: list_tested
      integer n, icell, face, bedgeid, ne, it, icell2
      real(kind=8) srchx, srchy, x1, y1, x2, y2, &
              faca, facb, vec_x, vec_y, dist
      logical inpoly, lboundary, converged

      ncorners = mesh_opp % nconnect
      !
      ! we have encountered a corner along the segment going through last_add and oppo_add
      ! create a list of already tested opposing cells from neighbour cells, excluding the two precious ones
      !
      ne = mesh_opp % ne
      allocate( list_tested( ne ) )
      list_tested(:) = .false.
      icell = last_add
      if (icell > 0) then
        do n = 1, ncorners
           face = mesh_opp % itm(n, icell)
           do it = 1, 2
              icell2 = mesh_opp % tseg(it, face)
              if( icell2 > 0 ) list_tested( icell2 ) = .true.
           enddo
        enddo
        list_tested(icell) = .false.
      endif
      icell = oppo_add
      if (icell > 0) then
        do n = 1, ncorners
           face = mesh_opp % itm(n, icell)
           do it = 1, 2
              icell2 = mesh_opp % tseg(it, face)
              if( icell2 > 0 ) list_tested( icell2 ) = .true.
           enddo
        enddo
        list_tested(icell) = .false.
      endif


      ! the new cell should be a neighbour of previous cells
      found = .false.
      converged = .false.
      inpoly = .false.
      next_add = 0
      x1 = begx; y1 = begy
      x2 = endx; y2 = endy

      ! choose a point just beyond "tiny" distance from x1,y2
      vec_x = x2 - x1
      vec_y = y2 - y1
      dist = vec_x * vec_x + vec_y * vec_y
      faca = 2.0_8 * sqrt( tiny / dist )
      if (faca > oned) then ! remaining segment too small, not worth continuing
         deallocate( list_tested )
         return
      endif
      facb = oned - faca
      srchx = facb * x1 + faca * x2
      srchy = facb * y1 + faca * y2

      allocate( cell_corner_x(ncorners), &
                cell_corner_y(ncorners))

      do icell = 1, ne

         if ( list_tested(icell) ) then
                               
           cell_corner_x(:) = mesh_opp % xgr( mesh_opp % in(:, icell) )
           cell_corner_y(:) = mesh_opp % ygr( mesh_opp % in(:, icell) )

           call ptinpoly(srchx, srchy, ncorners, &
     &           cell_corner_x, cell_corner_y, &
     &           inpoly, lboundary, bedgeid) 

           if ( inpoly ) exit ! leave loop over icell
                            
         endif

      enddo ! search all tagged elements

      if ( inpoly ) next_add = icell

      deallocate( cell_corner_x, cell_corner_y, list_tested )

      end subroutine find_other_poly


!***********************************************************************



      subroutine store_link_cnsrv(add1, add2, weights)

!-----------------------------------------------------------------------
!
!     this routine stores the address and weight for this link in
!     the appropriate address and weight arrays and resizes those
!     arrays if necessary.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer, intent(in) :: &
     &        add1, & ! address on src_grid
     &        add2    ! address on opp_grid

      real(kind=8), dimension(:), intent(in) :: &
     &        weights ! array of remapping weights for this link

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: nlink, min_link, max_link ! link index

      logical, save :: first_call = .true.

      logical :: found


!-----------------------------------------------------------------------
!
!     if all weights are zero, do not bother storing the link
!
!-----------------------------------------------------------------------

      if (all(weights == zerd)) return

!-----------------------------------------------------------------------
!
!     restrict the range of links to search for existing links
!
!-----------------------------------------------------------------------

!$OMP CRITICAL(block5)
!        min_link = min(uwgts % link_add(1,add1), opp_grid % wgt_p % link_add(1,add2))
!        max_link = max(uwgts % link_add(2,add1), opp_grid % wgt_p % link_add(2,add2))
        min_link = uwgts % link_add(1,add1)
        max_link = uwgts % link_add(2,add1)
        if (min_link == 0) then
          min_link = 1
          max_link = 0
        endif
!$OMP END CRITICAL(block5)

!-----------------------------------------------------------------------
!
!     if the link already exists, add the weight to the current weight
!     arrays
!
!-----------------------------------------------------------------------
      found = .false.

      do nlink=min_link,max_link
        if ( add1 == uwgts % src_cell(nlink) .and. &
             add2 == uwgts % dst_cell(nlink)) then

!$OMP CRITICAL(block3a)
          uwgts % wts_map(:,nlink) = uwgts % wts_map(:,nlink) + weights(1:num_wts)
!          opp_grid % wgt_p % wts_map(:,nlink) = opp_grid % wgt_p % wts_map(:,nlink) + weights(1+num_wts:2*num_wts)
!$OMP END CRITICAL(block3a)
          found = .true.
          exit

        endif
      end do


      if (found) return

!-----------------------------------------------------------------------
!
!     if the link does not yet exist, increment number of links and 
!     check to see if remap arrays need to be increased to accomodate 
!     the new link.  then store the link.
!
!-----------------------------------------------------------------------

!$OMP CRITICAL(block6)

      uwgts % num_links_map  = uwgts % num_links_map + 1
      nlink = uwgts % num_links_map
      if (nlink > uwgts % max_num_links) &
     &   call resize_remap_vars(resize_increment)

      uwgts % src_cell(nlink) = add1
      uwgts % dst_cell(nlink) = add2
      uwgts % wts_map    (:,nlink) = weights(1:num_wts)

      if (uwgts % link_add(1,add1) == 0) uwgts % link_add(1:2,add1) = nlink
      uwgts % link_add(1,add1) = min( uwgts % link_add(1,add1), nlink)
      uwgts % link_add(2,add1) = max( uwgts % link_add(2,add1), nlink)

!$OMP END CRITICAL(block6)

!-----------------------------------------------------------------------

      end subroutine store_link_cnsrv

!***********************************************************************

!***********************************************************************

      subroutine resize_remap_vars(increment)

!-----------------------------------------------------------------------
!
!     this routine resizes remapping arrays by increasing(decreasing)
!     the max_links by increment
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!     input variables
!
!-----------------------------------------------------------------------

      integer, intent(in) :: increment  ! the number of links to add(subtract) to arrays

!-----------------------------------------------------------------------
!
!     local variables
!
!-----------------------------------------------------------------------

      integer :: &
     &   ierr,   &  ! error flag
     &   mxlinks,&  ! size of link arrays
     &   nlink

      integer, dimension(:), allocatable :: &
     &   add1_tmp, & ! temp array for resizing address arrays
     &   add2_tmp    ! temp array for resizing address arrays

      real(kind=8), dimension(:,:), allocatable :: &
     &   wts_tmp   ! temp array for resizing weight arrays

!-----------------------------------------------------------------------
!
!     resize map 1 arrays if required.
!
!-----------------------------------------------------------------------


        !***
        !*** allocate temporaries to hold original values
        !***

        mxlinks = size(uwgts % src_cell)
        allocate (add1_tmp(mxlinks), add2_tmp(mxlinks), &
     &            wts_tmp(num_wts,mxlinks))

        add1_tmp = uwgts % src_cell
        add2_tmp = uwgts % dst_cell
        wts_tmp  = uwgts % wts_map
        
        !***
        !*** deallocate originals and increment max_links then
        !*** reallocate arrays at new size
        !***

        deallocate (uwgts % src_cell, &
                    uwgts % dst_cell, &
                    uwgts % wts_map   )
        nlink = mxlinks + increment
        allocate (uwgts % src_cell(nlink), &
     &            uwgts % dst_cell(nlink), &
     &            uwgts % wts_map(num_wts, nlink)   )

        !***
        !*** restore original values from temp arrays and
        !*** deallocate temps
        !***

        mxlinks = min(mxlinks, uwgts % max_num_links)
        uwgts % src_cell(1:mxlinks) = add1_tmp (1:mxlinks)
        uwgts % dst_cell(1:mxlinks) = add2_tmp (1:mxlinks)
        uwgts % wts_map    (:,1:mxlinks) = wts_tmp(:,1:mxlinks)
        deallocate(add1_tmp, add2_tmp, wts_tmp)

        uwgts % max_num_links = nlink

!-----------------------------------------------------------------------

      end subroutine resize_remap_vars

!********************************************************************
!********************************************************************

subroutine allocate_wgts(max_link)
! locals
integer max_link, min_link
integer ne

      if (.not. uwgts % alloc_wgts ) then

        uwgts % alloc_wgts = .true.

        ne = src_grid % ne
        allocate( uwgts % link_add(2, ne) )
        allocate( uwgts % src_cell(max_link) )
        allocate( uwgts % dst_cell(max_link) )
        allocate( uwgts % wts_map(num_wts, max_link) )

        uwgts % num_links_map = 0
        uwgts % max_num_links = max_link
        uwgts % link_add = 0
        uwgts % src_cell = 0
        uwgts % dst_cell = 0
        uwgts % wts_map = 0

      endif

end subroutine allocate_wgts

end module remap_cons_engine
