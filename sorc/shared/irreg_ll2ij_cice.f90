subroutine irreg_ll2ij_cice (nx, ny, grdlat, grdlon, npts, lats, lons, xi, yj)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  irreg_ll2ij_cice
!
! DESCRIPTION:  converts latitide, longitude positions to grid i,j
!               coordinates for an irregular grid
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!     Name        Type      Usage            Description
!   ---------   --------   -------   -------------------------------
!    nx         integer     input    number "x" grid positions
!    ny         integer     input    number "y" grid positions
!    grdlat     real        input    grid latitudes
!    grdlon     real        input    grid longitudes
!    lats       real        input    latitudes of points
!    lons       real        input    longitudes of points
!    npts       integer     input    number points to convert
!    xi         real        output   grid x position of points
!    yj         real        output   grid y position of points
!
!....................MAINTENANCE SECTION................................
!
! RECORD OF CHANGES:
!   Initial Installation - April 1994 -- Cummings, J.
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer, parameter :: max_neighbors=100000
!
!     ..local array dimensions
!
      integer nx, ny
      integer npts
!
      real grdlat(nx,ny), grdlon(nx,ny)
      real lons(npts), lats(npts)
      real xi(npts), yj(npts)
!
      integer ll, llon, llat
      integer i, iq
      integer i4(4), j4(4)
      integer j, k
      integer nij, nnij
      real plon, plat
!
      logical linside
!
!     ..allocatable arrays
!
      integer,  allocatable :: iref_llsort (:)
      integer,  allocatable :: index_llbin (:)
      real,     allocatable :: lcl_lat (:,:)
      real,     allocatable :: lcl_lon (:,:)
      integer,  allocatable :: ninbin (:)
      integer,  allocatable :: neighbor_index (:)
!
!...............................executable..............................
!
!     ..allocate arrays
!
      allocate (iref_llsort (nx * ny))
      allocate (index_llbin (181 * 360))
      allocate (lcl_lat (nx, ny))
      allocate (lcl_lon (nx, ny))
      allocate (ninbin (181 * 360))
      allocate (neighbor_index(0:max_neighbors))
!
!     ..create local grid lat,lon arrays since this routine destroys
!       the input grid arguments
!

!RG: PRINT *,'irreg',lats, lons
      do j = 1, ny
         do i = 1, nx
            lcl_lat(i,j) = grdlat(i,j)
            lcl_lon(i,j) = grdlon(i,j)
         enddo
      enddo
!
!     ..initialize
!
      do i = 0, max_neighbors
         neighbor_index(i) = 0
      enddo
!
      call make_llidx (nx, ny, lcl_lat, lcl_lon, iref_llsort, index_llbin, ninbin)
!
      do k=1,npts
!
        plon=lons(k)
        plat=lats(k)
!
!     define longitude 0-359.99
!
        do while (plon.lt.0.)
          plon = plon + 360.
        enddo
        do while (plon.ge.360.)
          plon = plon - 360.
        enddo
!
!     id potential surrounding grid cell corners based on lat,lon bins
!
        llon = max(int(plon-0.0001)+1,1)
        llat = max(int(plat-0.0001),-90)+91
        if (plat.lt.0.) llat = max(1,llat-1)
!
        ll = llon + (llat-1)*360
!
!     make list of indexes of potential neighbors
!
        nnij = 0
!
        if (llat.le.1.or.llat.ge.180) then
!
!       use all polar longitude bins
!
          llon = llon - 1
          if (llon.le.0) llon = 360
          llon = llon - 1
          if (llon.le.0) llon = 360
          do i = 1,360
            llon = llon + 1
            if (llon.gt.360) llon = 1
            ll = llon + (llat-1)*360
            if (ninbin(ll).gt.0) then                            
              do nij = 1,ninbin(ll)
                neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
              enddo
              nnij = nnij + ninbin(ll)
            endif
          enddo
          llon = min(int(plon-0.0001)+1,1)
!
        else
!
!       present lon,lat bin
!
          if (ninbin(ll).gt.0) then
            do nij = 1,ninbin(ll)
              neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
            enddo
            nnij = nnij + ninbin(ll)
          endif
!
!       append neighboring bins if close to line
!
          if (abs(plon-real(llon-1)).lt..4) then
!         append bin to west
            if (llon.eq.1) then
              ll = llat*360
            else
              ll = llon - 1 + (llat-1)*360
            endif
            if (ninbin(ll).gt.0) then
              do nij = 1,ninbin(ll)
                neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
              enddo
              nnij = nnij + ninbin(ll)
            endif
            if (abs(plat-real(llat-91)).lt..4.and.llat.gt.1) then
!           append bin to southwest
              if (llon.eq.1) then
                ll = (llat-1)*360
              else
                ll = llon - 1 + (llat-2)*360
              endif
              if (ninbin(ll).gt.0) then
                do nij = 1,ninbin(ll)
                  neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
                enddo
                nnij = nnij + ninbin(ll)
              endif
            elseif (abs(plat-real(llat-90)).lt..4.and.llat.lt.180) then
!           append bin to northwest
              if (llon.eq.1) then
                ll = (llat+1)*360
              else
                ll = llon - 1 + llat*360
              endif
              if (ninbin(ll).gt.0) then
                do nij = 1,ninbin(ll)
                  neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
                enddo
                nnij = nnij + ninbin(ll)
              endif
            endif
          elseif (abs(plon-real(llon)).lt..4) then
!         append bin to east
            if (llon.eq.360) then
              ll = 1 + (llat-1)*360
            else
              ll = llon + 1 + (llat-1)*360
            endif
            if (ninbin(ll).gt.0) then
              do nij = 1,ninbin(ll)
                neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
              enddo
              nnij = nnij + ninbin(ll)
            endif
            if (abs(plat-real(llat-91)).lt..4.and.llat.gt.1) then
!           append bin to southeast
              if (llon.eq.1) then
                ll = 1 + (llat-2)*360
              else
                ll = llon + 1 + (llat-2)*360
              endif
              if (ninbin(ll).gt.0) then
                do nij = 1,ninbin(ll)
                  neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
                enddo
                nnij = nnij + ninbin(ll)
              endif
            elseif (abs(plat-real(llat-90)).lt..4.and.llat.lt.180) then
!           append bin to northeast
              if (llon.eq.1) then
                ll = 1 + llat*360
              else
                ll = llon + 1 + llat*360
              endif
              if (ninbin(ll).gt.0) then
                do nij = 1,ninbin(ll)
                  neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
                enddo
                nnij = nnij + ninbin(ll)
              endif
            endif
          endif
          if (abs(plat-real(llat-91)).lt..4.and.llat.gt.1) then
!         append bin to south
            ll = llon + (llat-2)*360
            if (ninbin(ll).gt.0) then
              do nij = 1,ninbin(ll)
                neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
              enddo
              nnij = nnij + ninbin(ll)
            endif
          elseif (abs(plat-real(llat-90)).lt..4.and.llat.lt.180) then
!         append bin to north
            ll = llon + llat*360
            if (ninbin(ll).gt.0) then
              do nij = 1,ninbin(ll)
                neighbor_index(nnij+nij) = iref_llsort(index_llbin(ll)-1+nij)
              enddo
              nnij = nnij + ninbin(ll)
            endif
          endif
!
!         if (nnij.lt.1) then
!           write (6,*) 'WARNING: ',plon,plat,
!    *      ' not found in lat,lon bin'
!         endif
!
        endif 
!
!     now we have coordinates of nnij nearby points in neighbor_index
!     identify which four of these contain the point
!
        linside = .false.
!
        if(nnij.lt.1) then
!
          xi(k)=-1.
          yj(k)=-1.
!
        else
!
          do nij=1,nnij
!
!       check whether each quadrant has max,min limits which might
!       contain elon,elat
!
            j4(1) = int((neighbor_index(nij)-1)/nx)+1
            i4(1) = neighbor_index(nij)-(j4(1)-1)*nx
!
            do iq=1,4
              call insidegrid (plon, plat, nx, ny, lcl_lon, lcl_lat, iq, linside, i4, j4)
              if(linside) exit
            enddo
!
            if (linside) then
              call sort4 (i4, j4)
              call i4toxy (nx, ny, lcl_lon, lcl_lat, plon, plat, i4, j4, xi(k), yj(k))
              exit
            endif
          enddo 
        endif 
      enddo 
!
!     ..clean up
!
      deallocate (iref_llsort, index_llbin, lcl_lat, lcl_lon)
      deallocate (ninbin, neighbor_index)
!
      return
      end
      subroutine sort4 (i4, j4)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  sort4
!
! DESCRIPTION:  does some kind of a sort or something, someone said
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!    Name          Type       Usage            Description
!   -------     ----------    ------    ---------------------------
!
!....................MAINTENANCE SECTION................................
!
! METHOD:
!
! RECORD OF CHANGES:
!   Initial Installation - April 1994 -- Cummings, J.
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer i4(4),j4(4),i,j,iw,jw,n1,n2,n
!
!..............................executable...............................
!
!     sort by j
!
      do j = 2,4
        iw = i4(j)
        jw = j4(j)
        do i = j-1,1,-1
          if (j4(i).le.jw) go to 10
          i4(i+1) = i4(i)
          j4(i+1) = j4(i)
        enddo
        i = 0
 10     continue
        i4(i+1) = iw
        j4(i+1) = jw
      enddo
!
!     sort values with same j by i
!
      n1 = 0
      do while (n1.lt.4)
        n1 = n1 + 1
        n = 0
        n2 = n1+n+1
        do while (n2.le.4)
          if (j4(n2).eq.j4(n1)) then
            n = n+1
            n2 = n1+n+1
          else
            n2 = 100
          endif
        enddo
        if (n.gt.0) then
!         sort by j
          do j = n1+1,n1+n
            iw = i4(j)
            jw = j4(j)
            do i = j-1,n1,-1
              if (i4(i).le.iw) go to 20
              i4(i+1) = i4(i)
              j4(i+1) = j4(i)
            enddo
            i = n1-1
 20         continue
            i4(i+1) = iw
            j4(i+1) = jw
          enddo
        endif
        n1 = n1+n
      enddo
!
!     now switch the last two so that the set forms a continuous, simple
!     closed loop.
!
      i = 3
      iw = i4(i)
      jw = j4(i)
!
      i4(i) = i4(i+1)
      j4(i) = j4(i+1)
      i4(i+1) = iw
      j4(i+1) = jw
!
      return
      end 
      subroutine i4toxy (nx, ny, lons, lats, plon, plat, i4, j4, x, y)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  i4toxy
!
! DESCRIPTION:  does a i4toxy or something
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!    Name          Type       Usage            Description
!   -------     ----------    ------    ---------------------------
!
!....................MAINTENANCE SECTION................................
!
! METHOD:
!
! RECORD OF CHANGES:
!   Initial Installation - April 1994 -- Cummings, J.
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer nx,ny
      integer i4(4),j4(4)
      integer iter,max_iter
!
      real lons(nx,ny),lats(nx,ny)
      real plon,plat
      real x,y
      real converge
      real zero,one
!
!     parameter (max_iter=100,converge=1.0e-10)
      parameter (max_iter=10,converge=2.0e-7)
      parameter (zero = 0.0,one=1.0)
!
      real iguess, jguess
      real deli, delj
      real dth1, dth2, dth3
      real dph1, dph2, dph3
      real dthp, dphp
      real mat1, mat2, mat3, mat4
      real determinant
!
!..............................executable...............................
!
      iter = 0
!
!     iterate to find i,j for bilinear approximation
!
      dth1 = lats(i4(2),j4(2)) - lats(i4(1),j4(1))
      dth2 = lats(i4(4),j4(4)) - lats(i4(1),j4(1))
      dth3 = lats(i4(3),j4(3)) - lats(i4(2),j4(2)) - dth2
!
      dph1 = lons(i4(2),j4(2)) - lons(i4(1),j4(1))
      dph2 = lons(i4(4),j4(4)) - lons(i4(1),j4(1))
      dph3 = lons(i4(3),j4(3)) - lons(i4(2),j4(2))
!
      if (dph1 .gt.  360.) dph1 = dph1 - 360.0
      if (dph2 .gt.  360.) dph2 = dph2 - 360.0
      if (dph3 .gt.  360.) dph3 = dph3 - 360.0
      if (dph1 .lt. -350.) dph1 = dph1 + 360.0
      if (dph2 .lt. -350.) dph2 = dph2 + 360.0
      if (dph3 .lt. -350.) dph3 = dph3 + 360.0
!
      dph3 = dph3 - dph2
!
!     account for singularities in grid
!
      if(dth1.eq.0..and.dph1.eq.0.) then
        dth1=.01
        dph1=.01
      endif
      if(dth2.eq.0..and.dph2.eq.0.) then
        dth2=.01
        dph2=.01
      endif
!
      iguess = zero
      jguess = zero
      deli = 1.
      delj = 1.
!
      do while (iter.lt.max_iter.and.  converge.lt.abs(deli).and.converge.lt.abs(delj))
!
        iter = iter + 1
        dthp = plat - lats(i4(1),j4(1)) - dth1*iguess - dth2*jguess - dth3*iguess*jguess
        dphp = plon - lons(i4(1),j4(1))
!
        if (dphp .gt.  180.) dphp = dphp - 360.0
        if (dphp .lt. -180.) dphp = dphp + 360.0
!
        dphp = dphp - dph1*iguess - dph2*jguess - dph3*iguess*jguess
!
        mat1 = dth1 + dth3*jguess
        mat2 = dth2 + dth3*iguess
        mat3 = dph1 + dph3*jguess
        mat4 = dph2 + dph3*iguess
!
        determinant = mat1*mat4 - mat2*mat3
!
        deli = (dthp*mat4 - mat2*dphp)/determinant
        delj = (mat1*dphp - dthp*mat3)/determinant
!
        iguess = iguess + deli
        jguess = jguess + delj
      enddo
!
      if (iter .le. max_iter) then
!
!     successfully found i,j - compute weights
!
        x=real(i4(1))+iguess
        y=real(j4(1))+jguess
!
      else
        !orig: call error_exit('I4TOXY', 'exceed max iteration count')
        PRINT *,'I4TOXY', 'exceed max iteration count'
      endif
!
      return
      end
      subroutine make_llidx (nx, ny, lats, lons, iref_llsort, index_llbin, ninbin)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  make_llidx
!
! DESCRIPTION:  
!     reference grid points by latitude, longitude bins
!     reference by index number ij = i + (j-1)*nx
!     iref_llsort is a list of ij values corresponding to the
!       original grid sorted by longitude within 1-degree sorted
!       latitude bins
!     index_llbin is an array of starting ij values locating
!       iref_llsort indices for each bin ll
!     each bin ll is 1x1 degree with southwest corner coordinates
!       sw lat = -90+int((ll-1)/360)
!       sw lon = (ll-int(lat+90)*360-1)*360
!       llat = 1+int((ll-1)/360)
!       llon = ll-360*(llat-1)
!       sw lon = llon-1
!       sw lat = -91+llat
!       ll = llon + (llat-1)*360
!     ninbin is an array listing the number of iref_llsort points
!       in bin ll
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!    Name          Type       Usage            Description
!   -------     ----------    ------    ---------------------------
!
!....................MAINTENANCE SECTION................................
!
! RECORD OF CHANGES:  Programmer: Charlie Barron
!                     Date: 14 April 2000
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer nx,ny
      integer iref_llsort(nx*ny)
      integer index_llbin(181*360),ninbin(181*360)
      integer i,j,ij,ij1,ij2,ija,ijb
      integer ll,llat,llon,npts,nptsb
      integer ijmax
      integer nll
!
      real lats(nx,ny),lons(nx,ny)
!
!     ..allocatable arrays
!
      integer,  allocatable :: isort (:)
      integer,  allocatable :: iwrk (:)
      real,     allocatable :: r1lats (:)
      real,     allocatable :: r1lons (:)
      real,     allocatable :: rwork (:)
!
!..............................executable...............................
!
!     ..allocate arrays
!
      allocate (isort (nx * ny))
      allocate (iwrk (nx * ny))
      allocate (r1lats (nx * ny))
      allocate (r1lons (nx * ny))
      allocate (rwork (nx * ny))
!
!     ..initialize
!
      nll = 181*360
      do i = 1, nll
         index_llbin(i) = 0
         ninbin(i) = 0
      enddo
!     
      do i = 1, (nx*ny)
         iref_llsort(i) = 0
         isort(i) = 0
         iwrk(i) = 0
         r1lats(i) = 0.
         r1lons(i) = 0.
         rwork(i) = 0.
      enddo 
!
!     define longitudes 0-359.99
      do j = 1,ny
        do i = 1,nx
          do while (lons(i,j).lt.0.)
            lons(i,j) = lons(i,j) + 360.
          enddo
          do while (lons(i,j).ge.360.)
            lons(i,j) = lons(i,j) - 360.
          enddo
        enddo
      enddo
!
!     define an index number for each lat,lon at i+(j-1)*nx
!     sort these index numbers by corresponding latitude
!     within each 1-degree latitude bin, sort by longitude
!
      do j = 1,ny
        do i = 1,nx
          ij = i+(j-1)*nx
          iref_llsort(ij) = ij
          r1lons(ij) = lons(i,j)
          r1lats(ij) = lats(i,j)
        enddo
      enddo
!     sort index by latitude
      call indexx (nx*ny, nx*ny, r1lats, iref_llsort)
!     sort lats, lons by latitudes
      do ij = 1,nx*ny
        rwork(ij) = r1lons(iref_llsort(ij))
      enddo
      do ij = 1,nx*ny
        r1lons(ij) = rwork(ij)
      enddo
      do ij = 1,nx*ny
        rwork(ij) = r1lats(iref_llsort(ij))
      enddo
      do ij = 1,nx*ny
        r1lats(ij) = rwork(ij)
      enddo
!     now r1lons, r1lats and iref_llsort are all sorted by latitude
!     sort each latitude bin by longitude
!     note that when used, very high latitude bins should probably
!     be grouped together
      ij1 = 1
      ijmax = nx*ny
      do llat = 1,180
        ija = ij1
        ijb = ija
        if (llat.eq.180) ijb = ijmax
        do while (r1lats(ijb).lt.(-90+llat).and.ijb.lt.ijmax)
          ijb = ijb + 1
        enddo
        npts = ijb-ija
        if (ijb.eq.ijmax.and.r1lats(ijb).le.(-90+llat)) npts=npts+1
        if (npts.eq.0) then
!         no points or one point in this zonal band
          do llon = 1,360
            ll = llon + (llat-1)*360
            index_llbin(ll) = ija
            ninbin(ll) = 0
          enddo
        else
          if (npts.gt.1) then
!           sort by longitude within band
            call indexx (npts, npts, r1lons(ija), isort)
            do ij1 = 1,npts
              ij = ij1 + ija - 1
              iwrk(ij1) = iref_llsort(isort(ij1)+ija-1)
              rwork(ij1) = r1lons(isort(ij1)+ija-1)
            enddo
            do ij1 = 1,npts
              ij = ij1 + ija - 1
              iref_llsort(ij) = iwrk(ij1)
              r1lons(ij) = rwork(ij1)
              rwork(ij1) = r1lats(isort(ij1)+ija-1)
            enddo
            do ij1 = 1,npts
              ij = ij1 + ija - 1
              r1lats(ij) = rwork(ij1)
            enddo
          endif
          ij1 = ija
          do llon = 1,360
            ij2 = ij1
            if (llon.eq.360) ij2 = ijb
            do while (r1lons(ij2).lt.(0+llon).and.ij2.lt.ijb)
              ij2 = ij2 + 1
            enddo
            nptsb = ij2-ij1
            if (ij2.eq.ijb.and.r1lons(ij2).le.(0+llon)) nptsb=nptsb+1
            ll = llon + (llat-1)*360
            index_llbin(ll) = ij1
            ninbin(ll) = nptsb
            ij1 = ij2
          enddo
        endif
      enddo
      do j = 1,ny
        do i = 1,nx
          do while (lons(i,j).gt.180.)
            lons(i,j) = lons(i,j) - 360.
          enddo
        enddo
      enddo
!
!     ..clean up
!
      deallocate (isort, iwrk, r1lats, r1lons, rwork)
!
      return
      end
      subroutine insidegrid (plon, plat, nx, ny, lons, lats, iq, linside, i4, j4)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  insidegrid
!
! DESCRIPTION:  
!     this subroutine examines the quadrants iq around point i,j
!     to determine whether plon,plat lies within this quadrant
!     first by checking extrema of longitude and latitude, then by
!     doing a cross-product test
!          2 | 1
!          --|--
!          3 | 4
!
!     Quadrant boundary points are ordered to traverse the perimeter
!     counterclockwise, with the quadrant interior to the left
!     Assumes a right-hand logical coordinate system
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!    Name          Type       Usage            Description
!   -------     ----------    ------    ---------------------------
!
!....................MAINTENANCE SECTION................................
!
! RECORD OF CHANGES:    Programmer: Charlie Barron
!                       Date: 4-17-2000
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer nx,ny
      integer i,j,iq,n,nextn,i4(4),j4(4)
      integer :: longlobe=1
      real lonmin,lonmax,latmin,latmax
      real x(4),y(4),xtest,ytest
      real plat,plon
      real lats(nx,ny),lons(nx,ny)
      real vec1x,vec1y,vec2x,vec2y
      real cross_product,t
      logical linside,wraptest
      real rlat,rlon
      real, parameter :: pi   = 3.14159265358979323846d0
      real, parameter :: d180 = 180.d0
      real, parameter :: d2r  = pi/d180
!
!..............................executable...............................
!
      linside = .false.
      wraptest = .false.
!
!     define quadrant boundary points
!     i4(1),j4(1) assumed defined on input
!
      i = i4(1)
      j = j4(1)
      if (iq.eq.1) then
        i4(2) = i+1
        j4(2) = j
        i4(3) = i+1
        j4(3) = j+1
        i4(4) = i
        j4(4) = j+1
      elseif (iq.eq.2) then
        i4(2) = i
        j4(2) = j+1
        i4(3) = i-1
        j4(3) = j+1
        i4(4) = i-1
        j4(4) = j
      elseif (iq.eq.3) then
        i4(2) = i-1
        j4(2) = j
        i4(3) = i-1
        j4(3) = j-1
        i4(4) = i
        j4(4) = j-1
      else
        i4(2) = i
        j4(2) = j-1
        i4(3) = i+1
        j4(3) = j-1
        i4(4) = i+1
        j4(4) = j
      endif
!
!     ..ensure points index within grid limits
!
!     check at boundaries
      do n = 2,4
!       ..check j points
        if (j4(n).eq.0 .or. j4(n).gt.ny) then
          return
        endif
!       ..check i points
        if (i4(n).eq.0) then
          wraptest = .true.
          if (longlobe.eq.-1) then
!           nonoverlapping cyclic
            i4(n) = nx
          elseif (longlobe.eq.0) then
!           overlapping cyclic
            i4(n) = nx-1
          else
!           noncyclic, quadrant undefined
            return
          endif
        elseif (i4(n).gt.nx) then
          wraptest = .true.
          if (longlobe.eq.-1) then
!           nonoverlapping cyclic
            i4(n) = i4(n)-nx
          elseif (longlobe.eq.0) then
!           overlapping cyclic
            i4(n) = i4(n)-nx+1
          else
!           noncyclic, quadrant undefined
            return
          endif
        elseif (abs(lons(i4(1),j4(1))-lons(i4(n),j4(n))).gt.300.) then
          wraptest = .true.
        endif
      enddo
!
!     define longitude points, limiting changes in points to 10
!     degrees from initial point to account for global wraparound
!     and transforming near-polar points into a Lambert equivalent
!     azimuthal projection
!
      if (abs(plat).ge.89.) then
!       near polar, transform into a Lambert equivalent azimuthal
!         projection
        rlat = d2r*(45.0-0.5*plat)
        rlon = d2r*plon
        xtest = 2.*sin(rlat)*cos(rlon)
        ytest = 2.*sin(rlat)*sin(rlon)
        do n = 1,4
          rlat = d2r*(45.0-0.5*lats(i4(n),j4(n)))
          rlon = d2r*lons(i4(n),j4(n))
          x(n) = 2.*sin(rlat)*cos(rlon)
          y(n) = 2.*sin(rlat)*sin(rlon)
        enddo
      else
        xtest = plon
        ytest = plat
        do n = 1,4
          x(n) = lons(i4(n),j4(n))
          y(n) = lats(i4(n),j4(n))
        enddo
!       account for global wraparound
        if (wraptest) then
          if (abs(x(1)-xtest).gt.10.) then
            if (abs(x(1)-xtest-360.).lt.10.) then
              xtest=xtest+360.
            elseif (abs(x(1)-xtest+360.).lt.10.) then
              xtest=xtest-360.
            else
!             bad cell, return with linside false
              return
            endif
          endif
          do n = 2,4
            if (abs(x(1)-x(n)).gt.10.) then
              if (abs(x(1)-x(n)-360.).lt.10.) then
                x(n)=x(n)+360.
              elseif (abs(x(1)-x(n)+360.).lt.10.) then
                x(n)=x(n)-360.
              else
!               bad cell, return with linside false
                return
              endif
            endif
          enddo
        endif
      endif
!
!     find extrema of box
      lonmin = x(1)
      lonmax = x(1)
      latmin = y(1)
      latmax = y(1)
      do n = 2,4
        lonmin = min(lonmin,x(n))
        lonmax = max(lonmax,x(n))
        latmin = min(latmin,y(n))
        latmax = max(latmax,y(n))
      enddo
      if ((lonmax-lonmin).gt.300.) then
!       wraparound check
        t = lonmin
        lonmin = lonmax - 360.
        lonmax = t
      endif
      if (xtest.gt.lonmax) xtest = xtest-360.
!
!     test whether xtest,ytest are inside extrema
      if (.not.(lonmin.le.xtest.and.xtest.le.lonmax.and.  latmin.le.ytest.and.ytest.le.latmax)) then
!       not inside, return with linside false
        return
      endif
!
!     This grid cell is a likely candidate. Do the cross product test.
!     Note that this may fail for non-convex grid cells. See SCRIP
!     User's guide version 1.2, Philip W. Jones, Los Alamos National
!     Laboratory for more info.
!
      do n = 1,4
        nextn = mod(n,4)+1
!       !***
!       !*** here we take the cross product of the vector making
!       !*** up each box side with the vector formed by the vertex
!       !*** and search point.  if all the cross products are
!       !*** positive, the point is contained in the box.
!       !***
        vec1y = y(nextn) - y(n)
        vec1x = x(nextn) - x(n)
        vec2y = ytest - y(n)
        vec2x = xtest - x(n)
!
!       !***
!       !*** check for -180,180 crossings
!       !***
!
        if (vec1x .gt. 180.) then
          vec1x = vec1x - 360.
        else if (vec1x .lt. -180.) then
          vec1x = vec1x + 360.
        endif
        if (vec2x .gt. 180.) then
          vec2x = vec2x - 360.
        else if (vec2x .lt. -180.) then
          vec2x = vec2x + 360.
        endif
!
        cross_product = vec1x*vec2y - vec2x*vec1y
!
!       !***
!       !*** if cross product is less than zero, this cell
!       !*** doesn't work
!       !***
!
        if (cross_product .lt. 0.0 ) return
!
      enddo
!
!     if we make it here, the point is linside
      linside = .true.
!
      return
      end
      subroutine indexx (n, np, arrin, indx)
!
!.............................START PROLOGUE............................
!
! CONFIGURATION IDENTIFICATION:
!      $HeadURL$
!      @(#)$Id$
!
! MODULE NAME:  indexx
!
! DESCRIPTION:  sorting according to real array arrin
!
! LIBRARIES OF RESIDENCE:   ...ops/lib/libcoda.a
!
! PARAMETERS:
!    Name          Type       Usage            Description
!   -------     ----------    ------    ---------------------------
!
!....................MAINTENANCE SECTION................................
!
! RECORD OF CHANGES:
!   Initial Installation - April 1994 -- Cummings, J.
!
!..............................END PROLOGUE.............................
!
      implicit none
!
      integer n,np
      integer indx(np)
      integer i,j,l,indxt,ir
      real arrin(np)
      real q
!
!..............................executable...............................
!
      if (n.lt.2) return
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
!     end of indexx
      end 
