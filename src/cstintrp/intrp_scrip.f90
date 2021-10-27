module intrp_scrip
!-----------------------------------------------------------------------
! modules
!-----------------------------------------------------------------------
      use grids
      use mesh
      use mesh_grid ! conversion between mesh (unstructured) and grid (structured)
      use remap_cons_engine, only : uwgts, num_wts, src_grid, opp_grid, &
              grid_integrate, store_link_cnsrv, allocate_wgts, &
              n_search_scrip

contains
!-------------------------------------------------------------------------!
!--------------------- generate mesh info --------------------------------!
!-------------------------------------------------------------------------!

      subroutine scrip( locr, loci, grds, grdr, intr, ln_scrip_grad, naggrmax )

      implicit none

      ! arguments
      type(type_location), intent(in), target :: locr, loci ! forward and reverse location
      type(type_grid)   , intent(inout), target :: grds, grdr ! source and reference/destination grids
      type(type_weights), intent(inout), target :: intr     ! 
      logical, intent(in) :: ln_scrip_grad
      integer, intent(in) :: naggrmax

      ! locals
      integer, pointer :: iwgt(:,:,:), jwgt(:,:,:), nwgt(:,:)
      real(kind=8), pointer :: wgt(:,:,:)
      real(kind=8) mwgt
      logical,      dimension(:,:), pointer :: mask
      integer ilink, nlink, ocn_add, id, jd, dnx, dny, is, js, snx, sny, jk

      ! preparation of the unstructured meshes
      if (grds % std_grd % grtyp /= 'M') then
          call init_mesh( grds )
      else
          call get_latlon_xyz_centroid( grds % meshd )
      endif
      if (grdr % std_grd % grtyp /= 'M') then
          call init_mesh( grdr )
      else
          call get_latlon_xyz_centroid( grdr % meshd )
      endif
      call construct_tree_mesh_center( grdr % meshd, .true. )

      ! main call to the calculation of the overlapping surfaces
      src_grid => grds % meshd
      opp_grid => grdr % meshd
      call grid_integrate

      ! post-processing of the untructured weights
      if (grdr % std_grd % grtyp == 'M') call correct_cell_to_nodes_dst( grdr % meshd ) 
      call reduce_maxaggr( naggrmax, intr % numwgt, grdr % meshd % ne, uwgts )

      ! convert to structured weights
      call convert_uwgts_to_intrp( grds, grdr, intr, naggrmax )

      ! post-processing of the structured weights
      if (grds % std_grd % grtyp == 'M') call correct_cell_to_nodes_src( grdr, grds % meshd, intr, naggrmax ) 
      call normalize_weights( grdr, intr, naggrmax ) 

      end subroutine scrip

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      subroutine reduce_maxaggr( naggrmax, naggrtot, ne_opp, wgt )

!-----------------------------------------------------------------------
! count the max number of contributing grid points

!-----------------------------------------------------------------------
      implicit none

      ! arguments
      integer naggrmax, naggrtot, ne_opp
      type(type_weight_unstr) :: wgt

      ! locals
      integer nlink, ilink, ocn_add
      integer, allocatable, dimension(:) :: countcell

      write(*,*) 'user supplied max number of neighbour: ',naggrmax

      nlink = wgt % num_links_map

      allocate(countcell(ne_opp))
      countcell(:) = 0

      do ilink=1,nlink
        ocn_add = wgt % dst_cell(ilink)
        countcell(ocn_add) = countcell(ocn_add) + 1
      enddo

      naggrtot = MAXVAL ( countcell )
      write(*,*) 'naggrtot',naggrtot

      if ( naggrtot > naggrmax ) then
        write(*,*) 'naggrtot too large !!! will try to reduce it'
        call reduce_maxneig( ne_opp, naggrmax, naggrtot, wgt, countcell )
        countcell(:) = 0
        do ilink=1,nlink
          ocn_add = wgt % dst_cell(ilink)
          if ( wgt % wts_map(1,ilink) > 0.d0 ) countcell(ocn_add) = countcell(ocn_add) + 1
        enddo
        naggrtot = MAXVAL ( countcell )
        write(*,*) 'new naggrtot',naggrtot
      endif
      deallocate( countcell )
      
      end subroutine reduce_maxaggr

!**********************************************************************

!-----------------------------------------------------------------------
      subroutine reduce_maxneig(dst_size, maxneig, num_wgts, wgts, cellcount )
!
! action: reduce the number of neighbour points
!
!-----------------------------------------------------------------------

      implicit none
! arguments
      integer dst_size, maxneig, num_wgts
      type(type_weight_unstr) :: wgts
      integer, dimension(dst_size) :: cellcount
! locals
      integer idst_maxneig
      integer link, isrc, idst, ind, idst2, jnd, ibuf, ilink, jlink
      integer, allocatable, dimension(:,:) :: neigh_list
      integer, allocatable, dimension(:) :: neigh_nb
      integer, dimension(dst_size) :: pointer_pb
      double precision :: sumtot
      integer, dimension(num_wgts) :: pclosest


!-----------------------------------------------------------------------
! find the number of points of the destination grid which causes trouble and assign a pointer value
!-----------------------------------------------------------------------

      idst_maxneig = 0
      pointer_pb(:) = 0
      do idst=1,dst_size
         if ( cellcount(idst) > maxneig ) then 
            idst_maxneig = idst_maxneig + 1
            pointer_pb(idst) = idst_maxneig
         endif
      enddo

!-----------------------------------------------------------------------
! allocate and store the matrix index which are under scrutinity
!-----------------------------------------------------------------------

      allocate(neigh_list(num_wgts,idst_maxneig),neigh_nb(idst_maxneig))
      neigh_list(:,:) = 0
      neigh_nb(:) = 0

      do link = 1 , wgts % num_links_map
         isrc = wgts % src_cell(link)
         idst = wgts % dst_cell(link)
         if ( cellcount(idst) > maxneig ) then
           ind = pointer_pb(idst)
           neigh_nb(ind) = neigh_nb(ind) + 1
           neigh_list(neigh_nb(ind),ind) = link
         endif
      enddo

!-----------------------------------------------------------------------
! find smallest weights [desired_neig] for each destination grid point
!-----------------------------------------------------------------------
      do idst2 = 1, idst_maxneig
! redorder as function of weight value
        pclosest( 1:neigh_nb(idst2) ) = 1

        do ind = 1, neigh_nb(idst2)
              ilink = neigh_list(ind,idst2)
           do jnd = 1, neigh_nb(idst2)
              jlink = neigh_list(jnd,idst2)
              if (wgts % wts_map(1,jlink) > wgts % wts_map(1,ilink) ) then
                 pclosest(ind) = pclosest(ind) + 1
              endif
           enddo
        enddo
! nullify small weights beyond the desired number of neighbours
        do ind = 1, neigh_nb(idst2)
           link = neigh_list(ind,idst2)
           if ( pclosest(ind)> maxneig ) wgts % wts_map(1,link) = 0.d0
        enddo
      enddo ! loop of suspect destination points

      end subroutine reduce_maxneig
!-----------------------------------------------------------------------

      subroutine correct_cell_to_nodes_dst( mesh_p ) 
      !
      ! convert interpolation weights cell-centered to node-centered
      !
      ! arguments
      type(type_mesh),intent(in) :: mesh_p
      !
      ! locals
      integer, allocatable, dimension(:) :: src_cell, dst_cell
      real(kind=8), allocatable, dimension(:,:) :: wgts_tmp
      integer link, nlink, isrc, idst, ncorners, corner, node, max_link
      real(kind=8) :: wtmp(num_wts)

      ncorners = mesh_p % nconnect
      nlink = size( uwgts % src_cell )

      ! copy the existing unstructured weights
      allocate( src_cell(nlink) )          ; src_cell = uwgts % src_cell
      allocate( dst_cell(nlink) )          ; dst_cell = uwgts % dst_cell
      allocate( wgts_tmp(num_wts,nlink) )  ; wgts_tmp = uwgts % wts_map

      ! re-allocate the unstructured weights for faster processing
      ! mind though that there is a memory cost to that
      nlink = uwgts % num_links_map
      max_link = nlink * 2
      deallocate( uwgts % link_add, uwgts % src_cell, uwgts % dst_cell, uwgts % wts_map )
      uwgts % alloc_wgts = .false.
      call allocate_wgts( max_link )

      ! new value after distributing to element vertices
      do link = 1 , nlink
         isrc = src_cell  (link)
         idst = dst_cell  (link)
         wtmp = wgts_tmp(:,link) / real(ncorners, kind=8)
         do corner = 1, ncorners
           node = mesh_p % in(corner, idst)
           call store_link_cnsrv(isrc, node, wtmp)
         enddo
      enddo

      deallocate( src_cell, dst_cell, wgts_tmp )

      end subroutine correct_cell_to_nodes_dst
!-----------------------------------------------------------------------

      subroutine correct_cell_to_nodes_src( grid, mesh_p, intr, naggrmax ) 
      !
      ! convert interpolation weights cell-centered to node-centered
      !
      ! arguments
      type(type_grid),intent(in) :: grid
      type(type_mesh),intent(in) :: mesh_p
      type(type_weights), intent(inout), target :: intr
      integer, intent(in) :: naggrmax
      !
      ! locals
      integer, allocatable, dimension(:) :: iwgt, iwgt2
      real(kind=8), allocatable, dimension(:) :: wgt, wgt2
      real(kind=8) mwgt
      integer nwgt
      integer isrc, ncorners, corner, node
      integer dnx, dny, id, jd, jk, jk2
      logical found

      ncorners = mesh_p % nconnect
      nwgt = intr % numwgt
      allocate( iwgt (nwgt), wgt (nwgt) )
      nwgt = ncorners * intr % numwgt
      allocate( iwgt2(nwgt), wgt2(nwgt) )

      dnx = grid % nx
      dny = grid % ny

      do jd = 1, dny
         do id = 1, dnx
            nwgt = intr % nwgt(id,jd)
            iwgt(1:nwgt) = intr % iwgt(1:nwgt,id,jd)
             wgt(1:nwgt) = intr %  wgt(1:nwgt,id,jd)
            iwgt2(:) = 0
             wgt2(:) = zerd
            nwgt = 0
            do jk = 1, intr % nwgt(id,jd)
               icell = iwgt(jk)
               do corner = 1, ncorners
                  node = mesh_p % in( corner, icell )
                  found = .false.
                  jk2 = 0
                  do while (.not.found .and. jk2 < nwgt)
                     jk2 = jk2 + 1
                     if (iwgt(jk2) == node) then
                        found = .true.
                     endif
                  enddo
                  if (found) then
                      wgt2(jk2) = wgt2(jk2) + wgt(jk)
                  else
                     nwgt = nwgt + 1
                     iwgt2(nwgt) = node
                      wgt2(nwgt) = wgt(jk)
                  endif
               enddo
            enddo
            ! need to make sure that we are not overpassing the limit
            ! order in decreasing order the weights
            do jk = 1, nwgt-1
               do jk2 = jk+1, nwgt
                  if ( wgt2(jk2) > wgt2(jk) ) then
                      mwgt =  wgt2(jk2)
                      node = iwgt2(jk2)
                       wgt2(jk2) =  wgt2(jk)
                      iwgt2(jk2) = iwgt2(jk)
                       wgt2(jk ) = mwgt
                      iwgt2(jk ) = node
                  endif
               enddo
            enddo
            !
            ! truncate the series if too long and keep the largest weights
            nwgt = min(nwgt, naggrmax)
            !
            ! copy back the final weights
            intr % iwgt(1:nwgt,id,jd) = iwgt2(1:nwgt)
            intr % iwgt(nwgt+1:naggrmax,id,jd) = 1
            intr %  wgt(1:nwgt,id,jd) =  wgt2(1:nwgt)
            intr %  wgt(nwgt+1:naggrmax,id,jd) = zerd
            intr % nwgt(id,jd) = nwgt
         enddo
      enddo

      deallocate( iwgt, wgt, iwgt2, wgt2 )

      intr % numwgt = MAXVAL( intr % nwgt )
      
      end subroutine correct_cell_to_nodes_src
!-----------------------------------------------------------------------

      subroutine convert_uwgts_to_intrp( grds, grdr, intr, naggrmax )

      implicit none

      ! arguments
      type(type_grid)   , intent(inout), target :: grds, grdr ! source and reference/destination grids
      type(type_weights), intent(inout), target :: intr     ! 
      integer, intent(in) :: naggrmax

      ! locals
      integer, pointer :: iwgt(:,:,:), jwgt(:,:,:), nwgt(:,:)
      real(kind=8), pointer :: wgt(:,:,:)
      real(kind=8) mwgt
      integer ilink, nlink, ocn_add, id, jd, dnx, dny, is, js, snx, sny, jk


      dnx = grdr % nx
      dny = grdr % ny
      snx = grds % nx
      sny = grds % ny

      nwgt => intr % nwgt; nwgt = 0
      iwgt => intr % iwgt; iwgt = 1
      jwgt => intr % jwgt; jwgt = 1
       wgt => intr %  wgt; wgt  = zerd
      nlink = uwgts % num_links_map

      do ilink = 1, nlink

         mwgt  = uwgts % wts_map(1,ilink)

         ocn_add = uwgts % dst_cell(ilink) ! find destination cell in active mesh

         if ( mwgt > 0.0_8 ) then

            if (grdr % std_grd % grtyp /= 'M') then
               ocn_add = opp_grid % ps2pm(ocn_add) ! find destination cell in unmasked grid
            endif
            id = mod(ocn_add - 1,  dnx) + 1
            jd =    (ocn_add - 1) / dnx + 1
            jk = nwgt(id,jd)
            jk = jk + 1

            ocn_add = uwgts % src_cell(ilink)! find source cell in active mesh
            if (grds % std_grd % grtyp /= 'M') then
               ocn_add = src_grid % ps2pm(ocn_add) ! find source cell in unmasked grid
            else
               snx = src_grid % ne
            endif
            is = mod(ocn_add - 1,  snx) + 1
            js =    (ocn_add - 1) / snx + 1

             wgt(jk,id,jd) = mwgt
            iwgt(jk,id,jd) = is
            jwgt(jk,id,jd) = js
            nwgt(   id,jd) = jk

        endif ! check if weight is positive

      enddo

      intr % numwgt = MAXVAL( nwgt )
      
      end subroutine convert_uwgts_to_intrp
!-----------------------------------------------------------------------

      subroutine normalize_weights( grid, intr, naggrmax ) 
      !
      ! normalize the weights so that the sum for each destination points is 1
      !
      ! arguments
      type(type_grid),intent(in) :: grid
      type(type_weights), intent(inout), target :: intr
      integer, intent(in) :: naggrmax
      !
      ! locals
      real(kind=8), allocatable, dimension(:) :: wgt
      real(kind=8) mwgt
      integer dnx, dny, id, jd

      allocate( wgt (naggrmax) )

      dnx = grid % nx
      dny = grid % ny

      do jd = 1, dny
         do id = 1, dnx
            nwgt = intr % nwgt(id,jd)
             wgt(:) = intr %  wgt(:,id,jd)
            ! rescale all weights
            mwgt = SUM( wgt )
            if (mwgt > zerd ) wgt(:) = wgt(:) / mwgt
            !
            ! copy back the final weights
            intr %  wgt(:,id,jd) =  wgt(:)
         enddo
      enddo

      deallocate( wgt )

      intr % numwgt = MAXVAL( intr % nwgt )

      end subroutine normalize_weights
!-----------------------------------------------------------------------

!***********************************************************************
!-----------------------------------------------------------------------
end module
