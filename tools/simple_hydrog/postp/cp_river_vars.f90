program cp_river_vars

! =========================================================================
!   program reads lat, lon, cellarea, land_frac, and tocell fields, and
!     computes the fields: subA, travel, basin, celllength, internal
! =========================================================================

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: ntilmx= 6
integer, parameter :: ncrosmx= 10000, ncrosij= 2

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, ii, jj, ip, jp, varid2, varid3, varid4, varid5
integer :: varid6, varid7, i1, j1, k, nto, ktr, rcode2, kt0, n1, np
integer :: ntiles, ndm, latgid, longid, idp2, jdp2, varid8, varid9
integer :: ilat_edge, ilon_edge, i2, ncros
real :: pi, dtr, sum, mval_cella, mval_tocell, mval_landf, b1, t1
real :: t2, a1, a2, dlat, dlon, csum, clen
logical :: scale_suba, write_rivlen
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ntilmx)            :: itw, ite, its, itn
integer, dimension (ncrosmx)           :: tcros
integer, dimension (ncrosmx,ncrosij)   :: icros, jcros

real, dimension (ni_cells,nj_cells)    :: out_flow

character(len=100), dimension (ntilmx)  :: river_input_file


integer, allocatable, dimension (:)     :: idx_to, ilo
integer, allocatable, dimension (:,:,:) :: idx_to_grd, ibas
real, allocatable, dimension (:)        :: lat_idx, lon_idx, bas_area, lat_to, &
                                           lon_to, avlat, avlon
real, allocatable, dimension (:,:,:)    :: cell_a, tocell, land_fr, suba, &
                                           travel, cell_l, basin, rivlen
real, allocatable, dimension (:,:,:)    :: lat, latb, lon, lonb, arlat, &
     sin_lat, cos_lat, drn_idx

!!Recall Fortran row-major order. Below, 8., 4., 2. occupy the 1st column.
out_flow  = reshape((/ 8.,  4.,  2., 16.,  0., 1., 32., 64., 128. /),shape(out_flow))

write_rivlen= .true.
!write_rivlen= .false.

pi= 4.*atan(1.)
dtr= pi/180.

!do j= 1,nj_cells
!   write (6,'(3f8.0)') (out_flow(i,j), i= 1,ni_cells)
!enddo

open (10, file= 'out.cp_river_vars', form= 'formatted')

read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
close (5)

! ---------------------------------------------------------------------------
!  get lon and lat dims from first file -- should be identical for all files
! ---------------------------------------------------------------------------
rcode= NF_OPEN (trim(river_input_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open netcdf file"
    write (6,*) trim(river_input_file(1))
    stop 1
endif

start= 1 ; count= 1
rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 2
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jd)

allocate (lat_idx(jd))
start= 1 ;  count= 1 ;  count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat_idx)

ilat_edge= 0
rcode= nf_inq_varid (ncid, 'latb', latid)         ! number of lat edges
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_yb', latid)
    if (rcode2 /= 0) then
!        write (6,*) "WARNING: lat edges not found, skipping"
        jdp1= jd + 1
        go to 10
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jdp1)
ilat_edge= 1
10 continue
jdp2= jd + 2
!write (6,*) 'jd= ', jd, ', jdp1= ', jdp1, ', jdp2= ', jdp2

rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 3
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)

ilon_edge= 0
rcode= nf_inq_varid (ncid, 'lonb', lonid)         ! number of lon edges
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_xb', lonid)
    if (rcode2 /= 0) then
!        write (6,*) "WARNING: lon edges not found, skipping"
        idp1= id + 1
        go to 20
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), idp1)
ilon_edge= 1
20 continue
idp2= id + 2
!write (6,*) 'id= ', id, ', idp1= ', idp1, ', idp2= ', idp2

rcode= nf_close (ncid)

allocate (lat(idp2,jdp2,ntiles), lon(idp2,jdp2,ntiles), arlat(idp2,jdp2,ntiles))
allocate (latb(id,jdp1,ntiles), lonb(idp1,jd,ntiles))

allocate (cell_a(idp2,jdp2,ntiles), land_fr(idp2,jdp2,ntiles), tocell(idp2,jdp2,ntiles))

! ----------------------------------------------------------------------
!  now get lons and lats from all input river files
! ----------------------------------------------------------------------
do n= 1,ntiles

!   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
   write (10,'(i6,3x,a)') n, trim(river_input_file(n))

   rcode= NF_OPEN (trim(river_input_file(n)), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file" ; stop 10
   endif

   start= 1 ; count= 1

   if (ntiles == 1) then
!     regular grid
       rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
       if (rcode /= 0) then
           rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
           if (rcode2 /= 0) then
               write (6,*) "ERROR: cannot find lat variable" ; stop 20
           endif
       endif
       rcode= nf_inq_vardimid (ncid, latid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= jd) then
           write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 25
       endif

       count(1)= jd
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(2,2:jdp1,n))
       do i= 3,idp1
          lat(i,:,:)= lat(2,:,:)
       enddo

       if (ilat_edge == 1) then
           rcode= nf_inq_varid (ncid, 'latb', latid)         ! number of lat edges
           if (rcode /= 0) then
               rcode2 = nf_inq_varid (ncid, 'grid_yb', latid)
           endif
           rcode= nf_inq_vardimid (ncid, latid, dimids)
           rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
           if (ndm /= jdp1) then
               write (6,*) "ERROR: inconsistent latb dimension, ", ndm, jdp1 ;  stop 28
           endif
           count(1)= jdp1
           rcode= nf_get_vara_double (ncid, latid, start, count, latb(1,:,n))
           do i= 2,idp1
              latb(i,:,:)= latb(1,:,:)
           enddo
       endif

       rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
       if (rcode /= 0) then
           rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
           if (rcode2 /= 0) then
               write (6,*) "ERROR: cannot find lon variable" ; stop 30
           endif
       endif
       rcode= nf_inq_vardimid (ncid, lonid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= id) then
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 35
       endif

       count(1)= id
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2,n))
       do j= 3,jdp1
          lon(:,j,:)= lon(:,2,:)
       enddo

       if (ilon_edge == 1) then
           rcode= nf_inq_varid (ncid, 'lonb', lonid)         ! number of lon edges
           if (rcode /= 0) then
               rcode2 = nf_inq_varid (ncid, 'grid_xb', lonid)
           endif
           rcode= nf_inq_vardimid (ncid, lonid, dimids)
           rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
           if (ndm /= idp1) then
               write (6,*) "ERROR: inconsistent lonb dimension, ", ndm, idp1 ;  stop 38
           endif
           count(1)= idp1
           rcode= nf_get_vara_double (ncid, lonid, start, count, lonb(:,1,n))
           do j= 1,jdp1
              lonb(:,j,:)= lonb(:,1,:)
           enddo
        endif

   else

!     cubic sphere -- assume no edge data
       rcode= nf_inq_varid (ncid, 'y', latid)         ! number of lats
       if (rcode /= 0) then
           write (6,*) "ERROR: cannot find lat variable (y)" ; stop 20
       endif
       rcode= nf_inq_vardimid (ncid, latid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= id) then
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 25
       endif
       rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
       if (ndm /= jd) then
           write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 25
       endif

       start= 1 ;  count(1)= id ;  count(2)= jd
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(2:idp1,2:jdp1,n))

       rcode= nf_inq_varid (ncid, 'x', lonid)         ! number of lons
       if (rcode /= 0) then
           write (6,*) "ERROR: cannot find lon variable (x)" ; stop 30
       endif
       rcode= nf_inq_vardimid (ncid, lonid, dimids)
       rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
       if (ndm /= id) then
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, id ;  stop 35
       endif
       rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
       if (ndm /= jd) then
           write (6,*) "ERROR: inconsistent lon dimension, ", ndm, jd ;  stop 35
       endif

       start= 1 ;  count(1)= id ;  count(2)= jd
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2:jdp1,n))

   endif

! ----------------------------------------------------------------------
!  now read cell area, land frac, and tocell fields
! ----------------------------------------------------------------------

!   write (6,*) 'read cellarea'
   rcode= nf_inq_varid (ncid, 'land_area', varid)
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'cellarea', varid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find land_area/cellarea variable" ; stop 22
       endif
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_cella= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_cella)
   endif

   cell_a(:,:,n)= mval_cella
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_a(2:idp1,2:jdp1,n))

!   where (cell_a(:,:,n) == mval_cella) cell_a(:,:,n)= mval_mdl
   where (cell_a(:,:,n) == mval_cella) cell_a(:,:,n)= 0.


!   write (6,*) 'read land_frac'
   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_landf= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
   endif

   land_fr(:,:,n)= mval_landf
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_fr(2:idp1,2:jdp1,n))

   where (land_fr(:,:,n) == mval_landf) land_fr(:,:,n)= mval_mdl


!   write (6,*) 'read tocell'
   rcode= nf_inq_varid (ncid, 'tocell', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_tocell= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_tocell)
   endif

   tocell(:,:,n)= mval_tocell
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(2:idp1,2:jdp1,n))

   where (tocell(:,:,n) == mval_tocell) tocell(:,:,n)= mval_mdl

   rcode= nf_close (ncid)
enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         if (cell_a(i,j,n) /= mval_mdl) sum= sum + cell_a(i,j,n)
      enddo
   enddo
enddo

scale_suba = .true.
if (sum/1.e6 < 510064460.) scale_suba = .false.


allocate (suba(idp2,jdp2,ntiles), travel(idp2,jdp2,ntiles), cell_l(idp2,jdp2,ntiles))
allocate (cos_lat(idp2,jdp2,ntiles), sin_lat(idp2,jdp2,ntiles))
allocate (ibas(idp2,jdp2,ntiles), basin(idp2,jdp2,ntiles), idx_to_grd(idp2,jdp2,ntiles))
allocate (rivlen(idp2,jdp2,ntiles))

if (ntiles == 1) then
    allocate (idx_to(idp2*jdp2/4), bas_area(idp2*jdp2/4))
    allocate (lat_to(idp2*jdp2/4), lon_to(idp2*jdp2/4), ilo(idp2*jdp2/4))
    allocate (avlat(idp2*jdp2/4), avlon(idp2*jdp2/4))
else
    allocate (idx_to(idp2*jdp2), bas_area(idp2*jdp2))
    allocate (lat_to(idp2*jdp2), lon_to(idp2*jdp2), ilo(idp2*jdp2))
    allocate (avlat(idp2*jdp2), avlon(idp2*jdp2))
endif

! okay, now set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

if (ntiles == 1) then
    itw= 1 ;  ite= 1 ;  its= 1 ;  itn= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
else
! define tiles to the west, east, south, and north
    do n= 1,ntiles
       if (mod(n,2) == 0) then
           itw(n)= mod(n+5,ntiles) ; ite(n)= mod(n+2,ntiles)
           its(n)= mod(n+4,ntiles) ; itn(n)= mod(n+1,ntiles)
       else
           itw(n)= mod(n+4,ntiles) ; ite(n)= mod(n+1,ntiles)
           its(n)= mod(n+5,ntiles) ; itn(n)= mod(n+2,ntiles)
       endif
       if (itw(n) == 0) itw(n)= ntiles
       if (ite(n) == 0) ite(n)= ntiles
       if (its(n) == 0) its(n)= ntiles
       if (itn(n) == 0) itn(n)= ntiles
!       write (6,'("tile ",i4, ", itw= ",i4, ", ite= ",i4, ", its= ",i4, ", itn= ",i4)') &
!           n, itw(n), ite(n), its(n), itn(n)
    enddo

    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
endif


if (ntiles == 1) then
    do n= 1,ntiles
! compute edge data if neccessary
       if (ilat_edge /= 1) then
           do j= 2,jd
              latb(1:id,j,n)= 0.5*(lat(2:idp1,j,n) + lat(2:idp1,j+1,n))
           enddo
           latb(1:id,1,n)=    2.*lat(2:idp1,2,n)    - latb(1:id,2,n)
           latb(1:id,jdp1,n)= 2.*lat(2:idp1,jdp1,n) - latb(1:id,jd,n)
           where (latb(:,1,n)    < lat1) latb(:,1,n)=    lat1
           where (latb(:,jdp1,n) > lat2) latb(:,jdp1,n)= lat2
       endif

       if (ilon_edge /= 1) then
          do i= 1,idp1
              lonb(i,1:jd,n)= 0.5*(lon(i,2:jdp1,n) + lon(i+1,2:jdp1,n))
           enddo
           lonb(1,1:jd,n)=    2.*lon(2,2:jdp1,n)    - lonb(2,1:jd,n)
           lonb(idp1,1:jd,n)= 2.*lon(idp1,2:jdp1,n) - lonb(id,1:jd,n)
           where (lonb(1,:,n)    < lon1) lonb(1,:,n)=    lon1
           where (lonb(idp1,:,n) > lon2) lonb(idp1,:,n)= lon2
       endif

! compute area of latitude
       do j= 1,jd
          do i= 1,id
             dlon= lonb(i+1,j,n) - lonb(i,j,n)
             if (dlon > 180.)  dlon= dlon - 360.
             if (dlon < -180.) dlon= dlon + 360.

             arlat(i+1,j+1,n)= erad*erad*abs(dlon)*dtr* &
                    (sin(latb(i,j+1,n)*dtr) - sin(latb(i,j,n)*dtr))
          enddo
       enddo
    enddo
endif

if ((ilat_edge == 1 .and. ilon_edge == 1) .or. ntiles == 1) then
    sum= 0.
    do n= 1,ntiles
       do j= 2,jdp1
          do i= 2,idp1
             sum= sum + arlat(i,j,n)
          enddo
       enddo
       write (6,*) 'sum= ', sum
    enddo
endif

deallocate (arlat)

! ----------------------------------------------------------------------
! compute sin and cos of latitudes
! ----------------------------------------------------------------------
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         sin_lat(i,j,n)= sin(lat(i,j,n)*dtr)
         cos_lat(i,j,n)= cos(lat(i,j,n)*dtr)
      enddo
   enddo
enddo

! ----------------------------------------------------------------------
! follow river downstream; get subA, travel, cell_length, and basin
! ----------------------------------------------------------------------

! first find all the grid cells where tocell=0
ktr= 0 ;  idx_to= -1 ;  idx_to_grd= -1 ;  ilo= -1
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 60
         if (tocell(i,j,n) == 0.) then
             ktr= ktr + 1
             idx_to(ktr)= ktr
             idx_to_grd(i,j,n)= ktr
             lat_to(ktr)= lat(i,j,n)
             lon_to(ktr)= lon(i,j,n)
             if (land_fr(i,j,n) >= 1.) then
                 ilo(ktr)= 0
             else
                 ilo(ktr)= 1
             endif
         endif
60       continue
      enddo
   enddo
enddo
!write (6,*) 'number of grid cells where tocell is zero= ', ktr
nto= ktr

suba= 0. ; travel= 0. ;  cell_l= mval_mdl ;  ibas= -1 ;  kt0= 0 ;  rivlen= 0.
do n= 1,ntiles
   do j= 2,jdp1
!      if (mod(j,200) == 0) write (6,*) 'j= ', j, ', tile= ', n
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl) go to 170
!             write (6,*) n, j, i, tocell(i,j,n)
             i1= i ; j1= j ; n1= n
             ktr= 0 ;  csum= 0.
!             if (i == 221 .and. j == 290 .and. n == 3) then
!                 write (6,'(a,8i6,f7.0,i7,f7.0)') '0 ', i1, j1, n1, ite(n1), i, j, n, ite(n), &
!                       tocell(i1,j1,n1), ktr, travel(i,j,n)
!             endif

             if (cell_a(i,j,n) == mval_mdl) then
                 write (6,'(a,2i5,2f10.3,2f10.0)') 'cell_a is missing value, ', j, i, &
                     lat(j,j,n), lon(i,j,n), tocell(i,j,n), cell_a(i,j,n)
                 stop 120
             endif

120          continue
             do jj= 1,nj_cells
                jp= j1+jj-2

                do ii= 1,ni_cells
                   ip=i1+ii-2

!                if (i == 712 .and. j == 247 .and. n == 1) then
!                    write (6,'(a,8i6,3f7.0)') '1 ', i1, j1, n1, ite(n1), i, j, n, ite(n), &
!                          tocell(i1,j1,n1), out_flow(ii,jj), travel(i1,j1,n1)
!                endif

                if (tocell(i1,j1,n1) == out_flow(ii,jj)) then
                    ktr= ktr + 1

                    if (ip /= i1 .or. jp /= j1) then
!                        clen= acos(sin_lat(i,j,n)*sin_lat(ip,jp,n) + &
!                         cos_lat(i,j,n)*cos_lat(ip,jp,n)*cos((abs(lon(i,j,n)-lon(ip,jp,n)))*dtr))* &
!                         erad*1000.
                        dlat= (lat(i1,j1,n)-lat(ip,jp,n))*dtr
                        dlon= (lon(i1,j1,n)-lon(ip,jp,n))*dtr
                        clen=      2.*asin( sqrt( (sin(dlat*0.5))**2. + &
                           cos_lat(i1,j1,n)*cos_lat(ip,jp,n)*(sin(dlon*0.5))**2. ) )* &
                           erad*1000.
                    else
                        clen= 0.
                    endif
                    if (ktr == 1) then
                        cell_l(i,j,n)= clen
                    endif
                    csum= csum + clen

                    if (scale_suba) then
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)*land_fr(i,j,n)
                    else
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)
                    endif

                    if (tocell(i1,j1,n1) == 0.) then
                        travel(i,j,n)= real(ktr-1)
                        ibas(i,j,n)= idx_to_grd(i1,j1,n1)
                        rivlen(i,j,n)= csum-clen
                        go to 170
                    else
                        i1= ip ; j1= jp ;  np= n1
                        if ( (ip == idp2 .and. (jp == 1 .or. jp == jdp2)) .or. &
                             (ip == 1    .and. (jp == 1 .or. jp == jdp2)) .or. &
                             (jp == 1    .and. (ip == 1 .or. ip == idp2)) .or. &
                             (jp == jdp2 .and. (ip == 1 .or. ip == idp2))) then
                               write (6,*) "WARNING: i and j edge, ", i, j, n, ip, jp
                        endif
                        if (ntiles == 1) then
                            if (ip == idp2) then
                                i1= 2
                            else if (ip == 1) then
                                i1= idp1
                            endif
                        else
                            if (mod(n1,2) == 0) then
                                if (ip == idp2) then
                                    i1= idp1-jp+2 ; j1= 2          ; n1= ite(np)
                                else if (ip == 1) then
                                    i1= idp1      ; j1= jp         ; n1= itw(np)
                                endif
                                if (jp == jdp2) then
                                    i1= ip        ; j1= 2          ; n1= itn(np)
                                else if (jp == 1) then
                                    i1= idp1      ; j1= jdp1-ip+2  ; n1= its(np)
                                endif
                            else
                                if (ip == idp2) then
                                    i1= 2         ; j1= jp         ; n1= ite(np)
                                else if (ip == 1) then
                                    i1= idp1-jp+2 ; j1= jdp1       ; n1= itw(np)
                                endif
                                if (jp == jdp2) then
                                    i1= 2         ; j1= jdp1-ip+2  ; n1= itn(np)
                                else if (jp == 1) then
                                    i1= ip        ; j1= jdp1       ; n1= its(np)
                                endif
                            endif
                        endif
                        go to 120
                    endif

                else if (tocell(i1,j1,n1) == mval_mdl) then

          ! tocell undefined, set to zero, set other vars accordingly
                    write (10,'(a,3i5,2f10.3,3i5,2f10.3,f15.1,f7.0)') "WARNING: tocell has missing value, ", &
                       j, i, n, lat(i,j,n), lon(i,j,n), j1, i1, n1, lat(i1,j1,n1), lon(i1,j1,n1), &
                       cell_a(i1,j1,n1), land_fr(i1,j1,n1)
                    kt0= kt0 + 1
                    if (cell_a(i1,j1,n1) == mval_mdl) cell_a(i1,j1,n1)= 0.
                    tocell(i1,j1,n1)= 0.
                    cell_l(i1,j1,n1)= 0.
                    travel(i1,j1,n1)= 0.
                    rivlen(i1,j1,n1)= 0.
                    if (scale_suba) then
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)*land_fr(i,j,n)
                    else
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)
                    endif
                    if (land_fr(i1,j1,n1) == mval_mdl) then
                        write (6,'(a,3i5,2f10.3)') "WARNING: land_frac has missing value, ", &
                           j1, i1, n1, lat(i1,j1,n1), lon(i1,j1,n1)
                        land_fr(i1,j1,n1)= 0.
                    endif

                    if (idx_to_grd(i1,j1,n1) == -1) then
                        nto= nto + 1
                        idx_to(nto)= nto
                        lat_to(nto)= lat(i1,j1,n1)
                        lon_to(nto)= lon(i1,j1,n1)
                        if (land_fr(i1,j1,n1) >= 1.) then
                            ilo(nto)= 0
                        else
                            ilo(nto)= 1
                        endif
                        idx_to_grd(i1,j1,n1)= nto
                        ibas(i,j,n)= nto
                        ibas(i1,j1,n1)= nto
                    else
                        write (6,*) 'here i am'
                        ibas(i,j,n)= idx_to_grd(i1,j1,n1)
                        ibas(i1,j1,n1)= idx_to_grd(i1,j1,n1)
                    endif
                    travel(i,j,n)= real(ktr)  ! ktr has not been incremented, no need to subtract 1
                    rivlen(i,j,n)= csum
                    go to 170
                endif

                enddo
             enddo

170      continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop
write (6,'(a,i6)') 'number of grid cells accepting drainage= ', nto
!write (6,'(a,i6)') 'number of grid cells where tocell has been changed from missing value to 0= ', kt0

where (tocell == mval_mdl)
   suba= mval_mdl
   travel= mval_mdl
   cell_l= mval_mdl
   rivlen= mval_mdl
endwhere

!where (tocell == 0.)
!   travel= 0.
!endwhere


! check for crossing of tocell paths
ktr= 0
if (ntiles == 1) then
   do n= 1,ntiles
      do j= 1,jdp1
         do i= 1,idp1
            if ((tocell(i,j,n) == 128. .and. tocell(i,j+1,n) ==  2.) .or. &
                (tocell(i,j,n) ==  32. .and. tocell(i,j+1,n) ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j+1
                     tcros(ktr)= n
            endif
            if ((tocell(i,j,n) ==   2. .and. tocell(i+1,j,n) ==  8.) .or. &
                (tocell(i,j,n) == 128. .and. tocell(i+1,j,n) == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i+1
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j
                     tcros(ktr)= n
            endif
         enddo
      enddo
   enddo
else
   do n= 1,ntiles
      do j= 2,jd
         do i= 2,idp1
            if ((tocell(i,j,n) == 128. .and. tocell(i,j+1,n) ==  2.) .or. &
                (tocell(i,j,n) ==  32. .and. tocell(i,j+1,n) ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j+1
                     tcros(ktr)= n
            endif
         enddo
      enddo
      do j= 2,jdp1
         do i= 2,id
            if ((tocell(i,j,n) ==   2. .and. tocell(i+1,j,n) ==  8.) .or. &
                (tocell(i,j,n) == 128. .and. tocell(i+1,j,n) == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i+1
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j
                     tcros(ktr)= n
            endif
         enddo
      enddo
   enddo

   do n= 1,ntiles
      if (mod(n,2) == 0) then
          do i= 1,idp2
             if ((tocell(i,1,n)    ==   2. .and. tocell(i,2,n)    ==  2.) .or. &
                 (tocell(i,1,n)    == 128. .and. tocell(i,2,n)    ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i
                     jcros(ktr,1)= 1 ;  jcros(ktr,2)= 2
                     tcros(ktr)= n
             endif
             if ((tocell(i,jdp1,n) == 128. .and. tocell(i,jdp2,n) ==  2.) .or. &
                 (tocell(i,jdp1,n) ==  32. .and. tocell(i,jdp2,n) ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i    ;  icros(ktr,2)= i
                     jcros(ktr,1)= jdp1 ;  jcros(ktr,2)= jdp2
                     tcros(ktr)= n
             endif
          enddo
          do j= 1,jdp2
             if ((tocell(1,j,n)    ==   2. .and. tocell(2,j,n)    ==  8.) .or. &
                 (tocell(1,j,n)    == 128. .and. tocell(2,j,n)    == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= 1 ;  icros(ktr,2)= 2
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j
                     tcros(ktr)= n
             endif
             if ((tocell(idp1,j,n) ==   2. .and. tocell(idp2,j,n) ==  2.) .or. &
                 (tocell(idp1,j,n) == 128. .and. tocell(idp2,j,n) ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= idp1 ;  icros(ktr,2)= idp2
                     jcros(ktr,1)= j    ;  jcros(ktr,2)= j
                     tcros(ktr)= n
             endif
          enddo
      else
          do i= 1,idp2
             if ((tocell(i,1,n)    == 128. .and. tocell(i,2,n)    ==  2.) .or. &
                 (tocell(i,1,n)    ==  32. .and. tocell(i,2,n)    ==  8.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i ;  icros(ktr,2)= i
                     jcros(ktr,1)= 1 ;  jcros(ktr,2)= 2
                     tcros(ktr)= n
             endif
             if ((tocell(i,jdp1,n) == 128. .and. tocell(i,jdp2,n) ==  8.) .or. &
                 (tocell(i,jdp1,n) ==  32. .and. tocell(i,jdp2,n) == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= i    ;  icros(ktr,2)= i
                     jcros(ktr,1)= jdp1 ;  jcros(ktr,2)= jdp2
                     tcros(ktr)= n
             endif
          enddo
          do j= 1,jdp2
             if ((tocell(1,j,n)    == 128. .and. tocell(2,j,n)    ==  8.) .or. &
                 (tocell(1,j,n)    ==  32. .and. tocell(2,j,n)    == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= 1 ;  icros(ktr,2)= 2
                     jcros(ktr,1)= j ;  jcros(ktr,2)= j
                     tcros(ktr)= n
             endif
             if ((tocell(idp1,j,n) ==   2. .and. tocell(idp2,j,n) ==  8.) .or. &
                 (tocell(idp1,j,n) == 128. .and. tocell(idp2,j,n) == 32.)) then
                     ktr= ktr + 1
                     icros(ktr,1)= idp1 ;  icros(ktr,2)= idp2
                     jcros(ktr,1)= j    ;  jcros(ktr,2)= j
                     tcros(ktr)= n
             endif
          enddo
      endif
   enddo
endif
ncros= ktr
write (6,'(a,i6)') "number of cases of rivers crossing= ", ncros
if (ncros > ncrosmx) then
    write (6,*) "ERROR: ncros > ncrosmx ; ncrosmx must be increased"
    write (6,*) "...exiting..."
    stop 170
endif

write (10,'(/a,i6)') "number of cases of rivers crossing= ", ncros
if (ncros > 0) then
    do l= 1,ktr
       write (10,'(2(a,2i5),a,i5)') "ERROR: rivers cross, j= ", jcros(l,1)-1, jcros(l,2)-1, &
          ", i= ", icros(l,1)-1, icros(l,2)-1, ", n= ", tcros(l)
    enddo
endif

! compute basin area
bas_area= 0. ;  avlat= 0. ;  avlon= 0.
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (ibas(i,j,n) == -1) go to 180
         do l= 1,nto
            if (ibas(i,j,n) == idx_to(l)) then
                if (scale_suba) then
                    bas_area(l)= bas_area(l) + cell_a(i,j,n)*land_fr(i,j,n)
                    avlat(l)= avlat(l) + lat(i,j,n)*cell_a(i,j,n)*land_fr(i,j,n)
                    avlon(l)= avlon(l) + lon(i,j,n)*cell_a(i,j,n)*land_fr(i,j,n)
                else
                    bas_area(l)= bas_area(l) + cell_a(i,j,n)
                    avlat(l)= avlat(l) + lat(i,j,n)*cell_a(i,j,n)
                    avlon(l)= avlon(l) + lon(i,j,n)*cell_a(i,j,n)
                endif
                go to 180
            endif
         enddo
         write (6,*) "ERROR: ibas not found (1), ", ibas(i,j,n)
         stop 180
180      continue
      enddo
   enddo
enddo

where (bas_area > 0.)
   avlat= avlat/bas_area
   avlon= avlon/bas_area
endwhere

write (10,*)
write (10,*) 'unordered list of drainage cells:'
do l= 1,nto
   write (10,'(3i6,4f10.2)') l, idx_to(l), ilo(l), lon_to(l), lat_to(l), &
       avlon(l), avlat(l)
enddo

! sort basins
do j= 1, nto
   do i= 1,nto-1
      if (bas_area(i+1) > bas_area(i)) then
          b1= bas_area(i)
          bas_area(i)= bas_area(i+1)
          bas_area(i+1)= b1

          i1= idx_to(i)
          idx_to(i)= idx_to(i+1)
          idx_to(i+1)= i1

          i2= ilo(i)
          ilo(i)= ilo(i+1)
          ilo(i+1)= i2

          t1= lat_to(i)
          lat_to(i)= lat_to(i+1)
          lat_to(i+1)= t1

          t2= lon_to(i)
          lon_to(i)= lon_to(i+1)
          lon_to(i+1)= t2

          a1= avlat(i)
          avlat(i)= avlat(i+1)
          avlat(i+1)= a1

          a2= avlon(i)
          avlon(i)= avlon(i+1)
          avlon(i+1)= a2
      endif
   enddo
enddo

write (10,*)
write (10,*) 'ordered list of drainage cells, land:'
ktr= 0
do l= 1,nto
   if (ilo(l) == 0) then
       ktr= ktr + 1
       write (10,'(2i6,f12.0,4f10.2)') ktr, l, bas_area(l)/1.0e6, &
           lon_to(l), lat_to(l), avlon(l), avlat(l)
   endif
enddo
if (ktr > 0) then
    write (6,*)
    write (6,'(a,i5,a)') "*** WARNING: ", ktr, " internal drain points have been found! ***"
    write (6,*)
endif

write (10,*)
write (10,*) 'ordered list of drainage cells, ocean:'
ktr= 0
do l= 1,nto
   if (ilo(l) == 1) then
       ktr= ktr + 1
       write (10,'(2i6,f12.0,4f10.2)') ktr, l, bas_area(l)/1.0e6, &
           lon_to(l), lat_to(l), avlon(l), avlat(l)
   endif
enddo

ktr= 0
do l= 1,nto-1
   if (bas_area(l) == bas_area(l+1)) then
!       write (10,*) 'basins have equal area, ', l
       ktr= ktr + 1
!       go to 85
   endif
enddo
85 continue
!write (6,*) 'number of basins having the same area= ', ktr

!  assign basin values ;  create drainage flag, equal to 1
!     for internal basins, 0 for external basins

allocate (drn_idx(idp2,jdp2,ntiles))

basin= mval_mdl ;  drn_idx= mval_mdl
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (ibas(i,j,n) == -1) go to 190
         do l= 1,nto
            if (ibas(i,j,n) == idx_to(l)) then
                basin(i,j,n)= real(l)
                if (ilo(l) == 0) then
                    drn_idx(i,j,n)= 1.
                else if (ilo(l) == 1) then
                    drn_idx(i,j,n)= 0.
                else
                    write (6,*) "ERROR: invalid ilo, ", ilo(l), n, j, i
                    stop 190
                endif
                go to 190
            endif
         enddo
         write (6,*) "ERROR: ibas not found (2), ", ibas(i,j,n)
         stop 190
190      continue
     enddo
   enddo
enddo


! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'river_network.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (lon, lat)
   rcode= NF_DEF_DIM (ncid, 'grid_x',  id,   londim)
   rcode= NF_DEF_DIM (ncid, 'grid_y',  jd,   latdim)

!  create coordinate variables
   rcode= NF_DEF_VAR (ncid, 'grid_x',  NF_DOUBLE, 1, (/ londim /),  lonid)
   rcode= NF_DEF_VAR (ncid, 'grid_y',  NF_DOUBLE, 1, (/ latdim /),  latid)

!  create attributes for coordinate variables
!    longitude:
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 16, 'T-cell longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'cartesian_axis', 1, 'X')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 12, 'degrees_east')

!    latitude:
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 8, 'latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 13, 'degrees_north')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'cartesian_axis', 1, 'Y')

!    create data variable and attributes
   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'subA', NF_DOUBLE, 2, ndims, varid2)
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'long_name', 13, 'subbasin area')
   rcode= NF_PUT_ATT_TEXT (ncid, varid2, 'units', 2, 'm2')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid2, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'tocell', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 28, 'direction to downstream cell')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'travel', NF_DOUBLE, 2, ndims, varid4)
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 42, &
             'cells left to travel before reaching ocean')
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'basin', NF_DOUBLE, 2, ndims, varid7)
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'long_name', 14, 'river basin id')
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid7, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'cellarea', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 9, 'cell area')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 2, 'm2')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'celllength', NF_DOUBLE, 2, ndims, varid5)
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', 11, 'cell length')
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'land_frac', NF_DOUBLE, 2, ndims, varid6)
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', 13, 'land fraction')
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'internal', NF_DOUBLE, 2, ndims, varid8)
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'long_name', 22, 'internal drainage flag')
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid8, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   if (write_rivlen) then
       rcode= NF_DEF_VAR (ncid, 'rivlen', NF_DOUBLE, 2, ndims, varid9)
       rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'long_name', 12, 'river length')
       rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'units', 1, 'm')
       rcode= NF_PUT_ATT_DOUBLE (ncid, varid9, 'missing_value', NF_DOUBLE, 1, mval_mdl)
   endif

   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 12, 'degrees_east')

   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 13, 'degrees_north')

!  leave define mode
   rcode= NF_ENDDEF (ncid)

!  write coordinate data
   start= 1 ;  count= 1

   count(1)= id
   rcode= NF_PUT_VARA_DOUBLE (ncid, lonid, start, count, lon_idx)

   count(1)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latid, start, count, lat_idx)

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lon(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, lat(2:idp1,2:jdp1,n))

!    basin data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, basin(2:idp1,2:jdp1,n))

!    cell area data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, cell_a(2:idp1,2:jdp1,n))

!    subA data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, suba(2:idp1,2:jdp1,n))

!    tocell data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, tocell(2:idp1,2:jdp1,n))

!    travel data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, travel(2:idp1,2:jdp1,n))

!    cell length data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, cell_l(2:idp1,2:jdp1,n))

!    land fraction data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, land_fr(2:idp1,2:jdp1,n))

!    drainage flag data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid8, start, count, drn_idx(2:idp1,2:jdp1,n))

   if (write_rivlen) then
!    river length data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid9, start, count, rivlen(2:idp1,2:jdp1,n))
   endif

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

close (10)

deallocate (lat_idx, lon_idx)
deallocate (lat, latb, lon, lonb, sin_lat, cos_lat)
deallocate (idx_to, idx_to_grd, ilo)
deallocate (lat_to, lon_to, avlat, avlon)
deallocate (cell_a, tocell, land_fr, rivlen, drn_idx)
deallocate (suba, travel, cell_l, basin, ibas, bas_area)

contains


! ----------------------------------------------------------------------
subroutine create_halo (ntl, id, jd, itw, ite, its, itn, field)
! ----------------------------------------------------------------------

implicit none

integer :: n, i, j, ip1, ip2, jp1, jp2

integer, intent(in)     :: ntl, id, jd
integer, intent(in)     :: itw(:), ite(:), its(:), itn(:)

real, intent(inout)     :: field(:,:,:)

ip1= id + 1  ;  jp1= jd + 1
ip2= id + 2  ;  jp2= jd + 2


if (ntl == 1) then
    field(1,:,ntl)= field(ip1,:,ntl)
    field(ip2,:,ntl)= field(2,:,ntl)
else
    do n= 1,ntl
! define tiles to the west, east, south, and north
       if (mod(n,2) == 0) then
           do j= 2,jp1
              field(1,j,n)=     field(ip1,j,itw(n))  ! western edge
           enddo

           do j= 2,jp1
              i= ip1-j+2
              field(ip2,j,n)=  field(i,2,ite(n))     ! eastern edge
           enddo

           do i= 2,ip1
              j= jp1-i+2
              field(i,1,n)=     field(ip1,j,its(n))  ! southern edge
           enddo

           do i= 2,ip1
              field(i,jp2,n)=  field(i,2,itn(n))     ! northern edge
           enddo
       else
           do j= 2,jp1
              i= ip1-j+2
              field(1,j,n)=     field(i,jp1,itw(n))  ! western edge
           enddo

           do j= 2,jp1
              field(ip2,j,n)=  field(2,j,ite(n))     ! eastern edge
           enddo

           do i= 2,ip1
              field(i,1,n)=     field(i,jp1,its(n))  ! southern edge
           enddo

           do i= 2,ip1
              j= jp1-i+2
              field(i,jp2,n)=  field(2,j,itn(n))     ! northern edge
           enddo
       endif
    enddo

endif

end subroutine create_halo

end program cp_river_vars
