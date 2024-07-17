!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
!*
!* FRE-NCTools is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FRE-NCTools..  If not, see
!* <http://www.gnu.org/licenses/>.
!***********************************************************************
#define CHECK_NF_ERRSTAT(ierr) call nfu_check_err(ierr,__FILE__,__LINE__)

program rmv_parallel_rivers

! ===============================================================================
!   program reads lat, lon, cellarea, land_frac, tocell, subA, and basin fields.
!   runs iterative procedure using subA to eliminate cases of parallel rivers.
!   Tocell is changed so that it points to the grid-cell-neighbor having the
!     largest subA value, IF:
!       -- the basin does not change, OR IF
!       -- the grid cell pointed to is ocean or part-land, OR IF
!       -- the grid cell pointed to by the cell pointed to is ocean or part-land
!   River fields then are updated and iteration continues until fields (subA,
!   travel, basin, celllength) no longer change.
! ===============================================================================
use nfu_mod

implicit none

integer, parameter :: maxdims= 2
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: ntilmx= 6
integer, parameter :: nitermx= 100
integer, parameter :: ncrosmx= 10000, ncrosij= 2

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.
real, parameter :: add_const_to_ocean= 5.e15

include 'netcdf.inc'

integer :: i, j, n, l, m, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latdim, londim
integer :: ii, jj, ip, jp, varid2, varid3, varid4, varid5
integer :: varid6, varid7, i1, j1, k, nto, ktr, rcode2, kt0, n1, np
integer :: ntiles, ndm, latgid, longid, idp2, jdp2, ktr2, varid8, i2
integer :: niter, ichk, idp3, jdp3, idp4, jdp4, iii, jjj, iip, jjp
integer :: ncros, ipnew, jpnew
real :: pi, dtr, sum, mval_in, b1, t1, t2, a1, a2, dlat, dlon, smax
real :: smax1, sbamx
logical :: scale_suba, add_ocean_const
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname, cover_type_file

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (ntilmx)            :: itw, ite, its, itn
integer, dimension (ncrosmx)           :: tcros
integer, dimension (ncrosij,ncrosmx)   :: icros, jcros
integer, dimension (ncrosij)           :: ipv, jpv

real, dimension (ncrosij)              :: sba

real, dimension (ni_cells,nj_cells)    :: out_flow

character(len=4) :: vname_cvr= "frac"

character(len=100), dimension (ntilmx)  :: river_input_file

integer, allocatable, dimension (:)     :: idx_to, ilo
integer, allocatable, dimension (:,:,:) :: idx_to_grd, ibas
real, allocatable, dimension (:)        :: lat_idx, lon_idx, bas_area, lat_to, &
                                           lon_to, avlat, avlon
real, allocatable, dimension (:,:,:)    :: cell_a, tocell, land_fr, suba, &
                                           travel, cell_l, basin, suba_temp, &
                                           travel_temp, celll_temp, basin_temp,&
                                           suba_wo_const, drn_idx
real, allocatable, dimension (:,:,:)    :: lat, latb, lon, lonb, arlat, &
     sin_lat, cos_lat, tocell_new

!!Recall Fortran row-major order. Below, 8., 4., 2. occupy the 1st column.
out_flow  = reshape((/ 8.,  4.,  2., 16.,  0., 1., 32., 64., 128. /),shape(out_flow))

pi= 4.*atan(1.)
dtr= pi/180.

add_ocean_const= .false.

open (10, file= 'out.rmv_parallel_rivers', form= 'formatted')

read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
read (5,'(l1)') add_ocean_const
close (5)

!write (6,*) "add_ocean_const= ", add_ocean_const

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

jdp1= jd + 1
jdp2= jd + 2
jdp3= jd + 3
jdp4= jd + 4
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

idp1= id + 1
idp2= id + 2
idp3= id + 3
idp4= id + 4
!write (6,*) 'id= ', id, ', idp1= ', idp1, ', idp2= ', idp2

CHECK_NF_ERRSTAT(nf_close(ncid))

allocate (lat(idp4,jdp4,ntiles), lon(idp4,jdp4,ntiles), arlat(idp4,jdp4,ntiles))

allocate (cell_a(idp4,jdp4,ntiles), land_fr(idp4,jdp4,ntiles), tocell(idp4,jdp4,ntiles))
allocate (suba(idp4,jdp4,ntiles), basin(idp4,jdp4,ntiles))
allocate (suba_wo_const(idp4,jdp4,ntiles))

! ----------------------------------------------------------------------
!  now get lons and lats from all input river files
! ----------------------------------------------------------------------
do n= 1,ntiles

!   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
!   write (10,'(i6,3x,a)') n, trim(river_input_file(n))

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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(3,3:jdp2,n))
       do i= 4,idp2
          lat(i,:,:)= lat(3,:,:)
       enddo

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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(3:idp2,3,n))
       do j= 4,jdp2
          lon(:,j,:)= lon(:,3,:)
       enddo

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
       rcode= nf_get_vara_double (ncid, latid, start, count, lat(3:idp2,3:jdp2,n))

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
       rcode= nf_get_vara_double (ncid, lonid, start, count, lon(3:idp2,3:jdp2,n))

   endif

! ----------------------------------------------------------------------
!  now read cell area, land frac, tocell, suba, and basin fields
! ----------------------------------------------------------------------

   rcode= nf_inq_varid (ncid, 'land_area', varid)
   if (rcode /= 0) then
       rcode2 = nf_inq_varid (ncid, 'cellarea', varid)
       if (rcode2 /= 0) then
           write (6,*) "ERROR: cannot find land_area/cellarea variable" ; stop 22
       endif
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif

   cell_a(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_a(3:idp2,3:jdp2,n))

   where (cell_a(:,:,n) == mval_in) cell_a(:,:,n)= 0.


   rcode= nf_inq_varid (ncid, 'land_frac', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif

   land_fr(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_fr(3:idp2,3:jdp2,n))

   where (land_fr(:,:,n) == mval_in) land_fr(:,:,n)= mval_mdl


   rcode= nf_inq_varid (ncid, 'tocell', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif

   tocell(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(3:idp2,3:jdp2,n))

   where (tocell(:,:,n) == mval_in) tocell(:,:,n)= mval_mdl


   rcode= nf_inq_varid (ncid, 'subA', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif

   suba(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, suba(3:idp2,3:jdp2,n))

   where (suba(:,:,n) == mval_in) suba(:,:,n)= mval_mdl


   rcode= nf_inq_varid (ncid, 'basin', varid)
   if (rcode /= 0) then
       print *, "rcode in ncvid is ", rcode ;  stop 10
   endif

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_in= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_in)
   endif

   basin(:,:,n)= mval_in
   start= 1; count= 1; count(1)= id ; count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, basin(3:idp2,3:jdp2,n))

   where (basin(:,:,n) == mval_in) basin(:,:,n)= mval_mdl

   CHECK_NF_ERRSTAT(nf_close(ncid))
enddo

sum= 0.
do n= 1,ntiles
   do j= 1,jdp4
      do i= 1,idp4
         if (cell_a(i,j,n) /= mval_mdl) sum= sum + cell_a(i,j,n)
      enddo
   enddo
enddo

scale_suba = .true.
if (sum/1.e6 < 510064460.) scale_suba = .false.



allocate (travel(idp4,jdp4,ntiles), cell_l(idp4,jdp4,ntiles))
allocate (cos_lat(idp4,jdp4,ntiles), sin_lat(idp4,jdp4,ntiles))
allocate (tocell_new(idp4,jdp4,ntiles))
allocate (ibas(idp4,jdp4,ntiles), idx_to_grd(idp4,jdp4,ntiles))

allocate (idx_to(idp4*jdp4), bas_area(idp4*jdp4))
allocate (lat_to(idp4*jdp4), lon_to(idp4*jdp4), ilo(idp4*jdp4))
allocate (avlat(idp4*jdp4), avlon(idp4*jdp4))

allocate (suba_temp(idp4,jdp4,ntiles), travel_temp(idp4,jdp4,ntiles))
allocate (celll_temp(idp4,jdp4,ntiles), basin_temp(idp4,jdp4,ntiles))
allocate (drn_idx(idp4,jdp4,ntiles))

!  adjust tocell to follow max subA, except at coastal cells
!    reiterate until suba, celllength, travel, and basin no longer change

220 continue

suba_temp= suba ; travel_temp= mval_mdl ; celll_temp= mval_mdl ; basin_temp= basin
niter= 0
250 continue

niter= niter + 1
!write (6,*) "ITERATION ", niter
if (niter > nitermx) then
    write (6,*) 'ERROR: too many iterations, niter= ', niter
    go to 260
endif

!  add very large value to ocean cells
suba_wo_const= suba_temp
if ( add_ocean_const ) then
    where (land_fr < 1.) suba_temp= suba_temp + add_const_to_ocean
endif

! find all the grid cells where tocell=0 and land_frac /= 1
    ktr= 0 ;  ; ktr2=  0
    do n= 1,ntiles
       do j= 3,jdp2
          do i= 3,idp2
             if (tocell(i,j,n) == mval_mdl) go to 449
             if (tocell(i,j,n) == 0. .and. land_fr(i,j,n) < 1.) then
                 ktr= ktr + 1
             else if (tocell(i,j,n) == 0. .and. land_fr(i,j,n) == 1.) then
                 ktr2= ktr2 + 1
             endif
449           continue
          enddo
       enddo
    enddo

!    write (6,*) 'number of externally draining basins= ', ktr
!    write (6,*) 'number of internally draining basins= ', ktr2

! okay, now set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

if (ntiles == 1) then
    itw= 1 ;  ite= 1 ;  its= 1 ;  itn= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, suba_temp)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, basin_temp)
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
    enddo

    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, land_fr)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, cell_a)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, suba_temp)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, basin_temp)
endif


! ----------------------------------------------------------------------
! compute sin and cos of latitudes
! ----------------------------------------------------------------------
do n= 1,ntiles
   do j= 1,jdp4
      do i= 1,idp4
         sin_lat(i,j,n)= sin(lat(i,j,n)*dtr)
         cos_lat(i,j,n)= cos(lat(i,j,n)*dtr)
      enddo
   enddo
enddo

! ----------------------------------------------------------------------
! compute new tocell value using max suba
! ----------------------------------------------------------------------
tocell_new= mval_mdl
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2

         if (tocell(i,j,n) == mval_mdl) go to 270

!         if (i == 151 .and. j == 3 .and. n == 4) then
!            write (6,'("1.0",3i6,f7.0,f10.3,f7.0)') i, j, n, tocell(i,j,n), land_fr(i,j,n), basin_temp(i,j,n)
!         endif

         if (tocell(i,j,n) == 0.) then
             tocell_new(i,j,n)= 0.
             go to 270
         endif

         smax= -99999. ;  i1= -1 ;  j1= -1
         do jj= 1,nj_cells
            jp= j+jj-2

            do ii= 1,ni_cells
               ip= i+ii-2

!               if (i == 151 .and. j == 3 .and. n == 4) then
!                  write (6,'("2.0",5i6,2f7.0,f10.3,f12.1)') i, j, n, ii, jj, out_flow(ii,jj), &
!                        tocell(ip,jp,n), land_fr(ip,jp,n), suba_temp(ip,jp,n)/1.e6
!               endif
               if (ii == 2 .and. jj == 2) go to 269  ! without this check, can create internal drain cells (18 sep 2013)

               smax1= smax
               smax= max(smax,suba_temp(ip,jp,n))
               if (smax == suba_temp(ip,jp,n)) then
                   if (basin_temp(ip,jp,n) == basin_temp(i,j,n)) then
                       i1= ii  ;  j1= jj
                   else if (tocell(ip,jp,n) == 0. .and. land_fr(ip,jp,n) < 1.) then
                       i1= ii  ;  j1= jj
                   else
                      do jjj= 1,nj_cells
                         jjp= jp+jjj-2
                         do iii= 1,ni_cells
                            iip= ip+iii-2
                            if (tocell(iip,jjp,n) == 0. .and. land_fr(iip,jjp,n) < 1.) then
                                i1= ii ;  j1= jj
                                go to 268
                            endif
                         enddo
                      enddo
                      smax= smax1
                   endif
268                continue
!                   if (i == 151 .and. j == 3 .and. n == 4) then
!                      write (6,'("3.0",7i6,3f15.3,f7.0)') i, j, n, ii, jj, i1, j1, smax/1.e6, &
!                           suba_temp(ip,jp,n)/1.e6, basin_temp(ip,jp,n)
!                   endif
               endif
269            continue

            enddo
         enddo
         if (i1 > 0 .and. j1 > 0) then
             tocell_new(i,j,n)= out_flow(i1,j1)
         else
             tocell_new(i,j,n)= tocell(i,j,n)
         endif
270      continue
!         if (i == 3 .and. j == 263 .and. n == 1) then
!            write (6,'("4.0",5i6,f7.0,f10.3)') i, j, n, i1, j1, tocell_new(i,j,n)
!         endif

      enddo
   enddo
enddo

tocell= tocell_new

! find all the grid cells where tocell=0 and land_frac /= 1
    ktr= 0 ;  ; ktr2=  0
    do n= 1,ntiles
       do j= 3,jdp2
          do i= 3,idp2
             if (tocell(i,j,n) == mval_mdl) go to 49
             if (tocell(i,j,n) == 0. .and. land_fr(i,j,n) < 1.) then
                 ktr= ktr + 1
             else if (tocell(i,j,n) == 0. .and. land_fr(i,j,n) == 1.) then
                 ktr2= ktr2 + 1
             endif
49           continue
          enddo
       enddo
    enddo

!    write (6,*) 'number of externally draining basins= ', ktr
!    write (6,*) 'number of internally draining basins= ', ktr2


! ----------------------------------------------------------------------
! follow river downstream; get subA, travel, cell_length, and basin
! ----------------------------------------------------------------------

! first find all the grid cells where tocell=0
ktr= 0 ;  idx_to= -1 ;  idx_to_grd= -1 ;  ilo= -1
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
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

suba= 0. ; travel= 0. ;  cell_l= mval_mdl ;  ibas= -1 ;  kt0= 0
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
         if (tocell(i,j,n) == mval_mdl) go to 170
             i1= i ; j1= j ; n1= n
             ktr= 0
!             if (i == 17 .and. j == 50 .and. n == 1) then
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

!                if (i == 17 .and. j == 50 .and. n == 1) then
!                    write (6,'(a,8i6,3f7.0)') '1 ', i1, j1, n1, ite(n1), i, j, n, ite(n), &
!                          tocell(i1,j1,n1), out_flow(ii,jj), travel(i1,j1,n1)
!                endif

                if (tocell(i1,j1,n1) == out_flow(ii,jj)) then
                    ktr= ktr + 1

                    if (ktr == 1) then
                        if (ip /= i .or. jp /= j) then
                            dlat= (lat(i,j,n)-lat(ip,jp,n))*dtr
                            dlon= (lon(i,j,n)-lon(ip,jp,n))*dtr
                            cell_l(i,j,n)=      2.*asin( sqrt( (sin(dlat*0.5))**2. + &
                               cos_lat(i,j,n)*cos_lat(ip,jp,n)*(sin(dlon*0.5))**2. ) )* &
                               erad*1000.
                        else
                            cell_l(i,j,n)= 0.
                        endif
                    endif

                    if (scale_suba) then
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)*land_fr(i,j,n)
                    else
                        suba(i1,j1,n1)= suba(i1,j1,n1) + cell_a(i,j,n)
                    endif

                    if (tocell(i1,j1,n1) == 0.) then
                        travel(i,j,n)= real(ktr-1)
                        ibas(i,j,n)= idx_to_grd(i1,j1,n1)
                        go to 170
                    else
                        i1= ip ; j1= jp ;  np= n1
                        if ( (ip == idp3 .and. (jp == 2 .or. jp == jdp3)) .or. &
                             (ip == 2    .and. (jp == 2 .or. jp == jdp3)) .or. &
                             (jp == 2    .and. (ip == 2 .or. ip == idp3)) .or. &
                             (jp == jdp3 .and. (ip == 2 .or. ip == idp3))) then
!                               write (6,*) "WARNING: i and j edge, ", i, j, n, ip, jp
                        endif
                        if (ntiles == 1) then
                            if (ip == idp3) then
                                i1= 3
                            else if (ip == 2) then
                                i1= idp2
                            endif
                        else
                            if (mod(n1,2) == 0) then
                                if (ip == idp3) then
                                    i1= idp2-jp+3 ; j1= 3          ; n1= ite(np)
                                else if (ip == 2) then
                                    i1= idp2      ; j1= jp         ; n1= itw(np)
                                endif
                                if (jp == jdp3) then
                                    i1= ip        ; j1= 3          ; n1= itn(np)
                                else if (jp == 2) then
                                    i1= idp2      ; j1= jdp2-ip+3  ; n1= its(np)
                                endif
                            else
                                if (ip == idp3) then
                                    i1= 3         ; j1= jp         ; n1= ite(np)
                                else if (ip == 2) then
                                    i1= idp2-jp+3 ; j1= jdp2       ; n1= itw(np)
                                endif
                                if (jp == jdp3) then
                                    i1= 3         ; j1= jdp2-ip+3  ; n1= itn(np)
                                else if (jp == 2) then
                                    i1= ip        ; j1= jdp2       ; n1= its(np)
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
                        ibas(i,j,n)= idx_to_grd(i1,j1,n1)
                        ibas(i1,j1,n1)= idx_to_grd(i1,j1,n1)
                    endif
                    travel(i,j,n)= real(ktr)  ! ktr has not been incremented, no need to subtract 1
                    go to 170
                endif

                enddo
             enddo

170      continue
      enddo     ! end of i loop
   enddo        ! end of j loop
enddo           ! end of ntiles loop
!write (6,'(a,i6)') 'number of grid cells accepting drainage= ', nto
!write (6,'(a,i6)') 'number of grid cells where tocell has been changed from missing value to 0= ', kt0

if (scale_suba) then
    where (tocell == 0. .and. land_fr == 0. .and. suba == cell_a*land_fr) tocell= mval_mdl
else
    where (tocell == 0. .and. land_fr == 0. .and. suba == cell_a) tocell= mval_mdl
endif

where (tocell == mval_mdl)
   suba= mval_mdl
   travel= mval_mdl
   cell_l= mval_mdl
   ibas= -1
endwhere

suba_temp= suba_wo_const

ichk= 0
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
         if (suba(i,j,n) /= suba_temp(i,j,n) .or.     &
             travel(i,j,n) /= travel_temp(i,j,n) .or. &
             basin(i,j,n) /= basin_temp(i,j,n)   .or. &
             cell_l(i,j,n) /= celll_temp(i,j,n)) then
!                write (6,*) "FIELDS DIFFER, ", i, j, n
!                write (6,*) 'suba ',   suba(i,j,n),   suba_temp(i,j,n)
!                write (6,*) 'travel ', travel(i,j,n), travel_temp(i,j,n)
!                write (6,*) 'basin ',   basin(i,j,n), basin_temp(i,j,n)
!                write (6,*) 'cell_l ', cell_l(i,j,n), celll_temp(i,j,n)
                ichk= 1
                go to 240
         endif
      enddo
   enddo
enddo
240 continue

! compute basin area
bas_area= 0. ;  avlat= 0. ;  avlon= 0.
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
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

!write (10,*)
!write (10,*) 'ordered list of drainage cells, land:'
!ktr= 0
!do l= 1,nto
!   if (ilo(l) == 0) then
!       ktr= ktr + 1
!       write (10,'(2i6,f12.0,4f10.2)') ktr, idx_to(l), bas_area(l)/1.0e6, &
!           lon_to(l), lat_to(l), avlon(l), avlat(l)
!   endif
!enddo

!write (10,*)
!write (10,*) 'ordered list of drainage cells, ocean:'
!ktr= 0
!do l= 1,nto
!   if (ilo(l) == 1) then
!       ktr= ktr + 1
!       write (10,'(2i6,f12.0,4f10.2)') ktr, idx_to(l), bas_area(l)/1.0e6, &
!           lon_to(l), lat_to(l), avlon(l), avlat(l)
!   endif
!enddo

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

basin= mval_mdl ;  drn_idx= mval_mdl
do n= 1,ntiles
   do j= 3,jdp2
      do i= 3,idp2
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

if (ichk == 1) then
    suba_temp= suba
    travel_temp= travel
    basin_temp= basin
    celll_temp= cell_l
    go to 250
endif

260 continue

! check for crossing of tocell paths
!ktr= 0
!if (ntiles == 1) then
!   do n= 1,ntiles
!      do j= 2,jdp2
!         do i= 2,idp2
!            if ((tocell(i,j,n) == 128. .and. tocell(i,j+1,n) ==  2.) .or. &
!                (tocell(i,j,n) ==  32. .and. tocell(i,j+1,n) ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j+1
!                     tcros(ktr)= n
!            endif
!            if ((tocell(i,j,n) ==   2. .and. tocell(i+1,j,n) ==  8.) .or. &
!                (tocell(i,j,n) == 128. .and. tocell(i+1,j,n) == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i+1
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!            endif
!         enddo
!      enddo
!   enddo
!else
!   do n= 1,ntiles
!      do j= 3,jdp1
!         do i= 3,idp2
!            if ((tocell(i,j,n) == 128. .and. tocell(i,j+1,n) ==  2.) .or. &
!                (tocell(i,j,n) ==  32. .and. tocell(i,j+1,n) ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j+1
!                     tcros(ktr)= n
!            endif
!         enddo
!      enddo
!      do j= 3,jdp2
!         do i= 3,idp1
!            if ((tocell(i,j,n) ==   2. .and. tocell(i+1,j,n) ==  8.) .or. &
!                (tocell(i,j,n) == 128. .and. tocell(i+1,j,n) == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i+1
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!            endif
!         enddo
!      enddo
!   enddo
!
!   do n= 1,ntiles
!      if (mod(n,2) == 0) then
!          do i= 2,idp3
!             if ((tocell(i,2,n)    ==   2. .and. tocell(i,3,n)    ==  2.) .or. &
!                 (tocell(i,2,n)    == 128. .and. tocell(i,3,n)    ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i
!                     jcros(1,ktr)= 2 ;  jcros(2,ktr)= 3
!                     tcros(ktr)= n
!             endif
!             if ((tocell(i,jdp2,n) == 128. .and. tocell(i,jdp3,n) ==  2.) .or. &
!                 (tocell(i,jdp2,n) ==  32. .and. tocell(i,jdp3,n) ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i    ;  icros(2,ktr)= i
!                     jcros(1,ktr)= jdp2 ;  jcros(2,ktr)= jdp3
!                     tcros(ktr)= n
!             endif
!          enddo
!          do j= 2,jdp3
!             if ((tocell(2,j,n)    ==   2. .and. tocell(3,j,n)    ==  8.) .or. &
!                 (tocell(2,j,n)    == 128. .and. tocell(3,j,n)    == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= 2 ;  icros(2,ktr)= 3
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!             endif
!             if ((tocell(idp2,j,n) ==   2. .and. tocell(idp3,j,n) ==  2.) .or. &
!                 (tocell(idp2,j,n) == 128. .and. tocell(idp3,j,n) ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= idp2 ;  icros(2,ktr)= idp3
!                     jcros(1,ktr)= j    ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!            endif
!          enddo
!      else
!          do i= 2,idp3
!             if ((tocell(i,2,n)    == 128. .and. tocell(i,3,n)    ==  2.) .or. &
!                 (tocell(i,2,n)    ==  32. .and. tocell(i,3,n)    ==  8.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i ;  icros(2,ktr)= i
!                     jcros(1,ktr)= 2 ;  jcros(2,ktr)= 3
!                     tcros(ktr)= n
!             endif
!             if ((tocell(i,jdp2,n) == 128. .and. tocell(i,jdp3,n) ==  8.) .or. &
!                 (tocell(i,jdp2,n) ==  32. .and. tocell(i,jdp3,n) == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= i    ;  icros(2,ktr)= i
!                     jcros(1,ktr)= jdp2 ;  jcros(2,ktr)= jdp3
!                     tcros(ktr)= n
!             endif
!          enddo
!          do j= 2,jdp3
!             if ((tocell(2,j,n)    == 128. .and. tocell(3,j,n)    ==  8.) .or. &
!                 (tocell(2,j,n)    ==  32. .and. tocell(3,j,n)    == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= 2 ;  icros(2,ktr)= 3
!                     jcros(1,ktr)= j ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!             endif
!             if ((tocell(idp2,j,n) ==   2. .and. tocell(idp3,j,n) ==  8.) .or. &
!                 (tocell(idp2,j,n) == 128. .and. tocell(idp3,j,n) == 32.)) then
!                     ktr= ktr + 1
!                     icros(1,ktr)= idp2 ;  icros(2,ktr)= idp3
!                     jcros(1,ktr)= j    ;  jcros(2,ktr)= j
!                     tcros(ktr)= n
!             endif
!          enddo
!      endif
!   enddo
!endif
!ncros= ktr
!write (6,*) "number of cases of rivers crossing= ", ncros
!if (ncros > 0) then
!    do l= 1,ncros
!       write (6,'(4(a,2i5))') "ERROR: tocell paths cross, j= ", jcros(1,l)-2, jcros(2,l)-2, &
!          ", i= ", icros(1,l)-2, icros(2,l)-2, ", n= ", tcros(l)
!    enddo
!endif
!ncros= 0

call find_crossing_rivers (ntiles, id, jd, tocell, ncros, icros, jcros, tcros)
write (6,*) "number of cases of rivers crossing= ", ncros

if (ncros > 0) then
    do l= 1,ncros
       write (6,'(4(a,2i5))') "FIX: rivers cross, j= ", jcros(1,l)-2, jcros(2,l)-2, &
          ", i= ", icros(1,l)-2, icros(2,l)-2, ", n= ", tcros(l)
    enddo

    do l= 1,ncros
       n1=  tcros(l)
       do m= 1,ncrosij
          i1= icros(m,l) ; j1= jcros(m,l)
          do jj= 1,nj_cells
             jp= j1+jj-2
             do ii= 1,ni_cells
                ip=i1+ii-2
                if (tocell(i1,j1,n1) == out_flow(ii,jj)) then
                    sba(m)= suba(ip,jp,n1)
                    ipv(m)= ip
                    jpv(m)= jp
                endif
             enddo
          enddo
       enddo
       sbamx= max(sba(1),sba(2))
       ipnew= -1 ; jpnew= -1
!       write (6,*) "sbamx= ", sbamx
       do m= 1,ncrosij
!          write (6,*) m, ipv(m), jpv(m), sba(m)
          if (sba(m) == sbamx) then
              ipnew= ipv(m)
              jpnew= jpv(m)
          endif
       enddo
       do m= 1,ncrosij
          if (sba(m) /= sbamx) then
              i1= icros(m,l) ; j1= jcros(m,l)
              do jj= 1,nj_cells
                 jp= j1+jj-2
                 do ii= 1,ni_cells
                    ip=i1+ii-2
                    if (ip == ipnew .and. jp == jpnew) then
                        tocell(i1,j1,n1)= out_flow(ii,jj)
                        go to 290
                    endif
                 enddo
              enddo
290           continue
          endif
       enddo
    enddo

endif

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
   rcode= NF_PUT_VARA_DOUBLE (ncid, longid, start, count, lon(3:idp2,3:jdp2,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, latgid, start, count, lat(3:idp2,3:jdp2,n))

!    basin data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, basin(3:idp2,3:jdp2,n))

!    cell area data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, cell_a(3:idp2,3:jdp2,n))

!    subA data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid2, start, count, suba(3:idp2,3:jdp2,n))

!    tocell data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, tocell(3:idp2,3:jdp2,n))

!    travel data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, travel(3:idp2,3:jdp2,n))

!    cell length data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, cell_l(3:idp2,3:jdp2,n))

!    land fraction data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, land_fr(3:idp2,3:jdp2,n))

!    drainage flag data
   count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid8, start, count, drn_idx(3:idp2,3:jdp2,n))

!  close netcdf file
   CHECK_NF_ERRSTAT(nf_close(ncid))
enddo

deallocate (lat, lon, sin_lat, cos_lat)
deallocate (idx_to, idx_to_grd, ilo)
deallocate (lat_to, lon_to, avlat, avlon)
deallocate (cell_a, tocell, land_fr, tocell_new)
deallocate (suba, travel, cell_l, basin, ibas, bas_area)
deallocate (suba_temp, travel_temp, celll_temp, basin_temp)
deallocate (suba_wo_const)


contains

! ----------------------------------------------------------------------
subroutine create_halo (ntl, id, jd, itw, ite, its, itn, field)
! ----------------------------------------------------------------------

implicit none

integer :: n, i, j, ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4

integer, intent(in)     :: ntl, id, jd
integer, intent(in)     :: itw(:), ite(:), its(:), itn(:)

real, intent(inout)     :: field(:,:,:)

ip1= id + 1  ;  jp1= jd + 1
ip2= id + 2  ;  jp2= jd + 2
ip3= id + 3  ;  jp3= jd + 3
ip4= id + 4  ;  jp4= jd + 4


if (ntl == 1) then
    field(1,:,ntl)= field(ip1,:,ntl)
    field(2,:,ntl)= field(ip2,:,ntl)
    field(ip3,:,ntl)= field(3,:,ntl)
    field(ip4,:,ntl)= field(4,:,ntl)
else
    do n= 1,ntl
! define tiles to the west, east, south, and north
       if (mod(n,2) == 0) then
           do j= 3,jp2
              field(1,j,n)=     field(ip1,j,itw(n))  ! western edge
              field(2,j,n)=     field(ip2,j,itw(n))  ! western edge
           enddo

           do j= 3,jp2
              i= ip2-j+3
              field(ip3,j,n)=  field(i,3,ite(n))     ! eastern edge
              field(ip4,j,n)=  field(i,4,ite(n))     ! eastern edge
           enddo

           do i= 3,ip2
              j= jp2-i+3
              field(i,1,n)=     field(ip1,j,its(n))  ! southern edge
              field(i,2,n)=     field(ip2,j,its(n))  ! southern edge
           enddo

           do i= 3,ip2
              field(i,jp3,n)=  field(i,3,itn(n))     ! northern edge
              field(i,jp4,n)=  field(i,4,itn(n))     ! northern edge
           enddo
       else
           do j= 3,jp2
              i= ip2-j+3
              field(1,j,n)=     field(i,jp1,itw(n))  ! western edge
              field(2,j,n)=     field(i,jp2,itw(n))  ! western edge
           enddo

           do j= 3,jp2
              field(ip3,j,n)=  field(3,j,ite(n))     ! eastern edge
              field(ip4,j,n)=  field(4,j,ite(n))     ! eastern edge
           enddo

           do i= 3,ip2
              field(i,1,n)=     field(i,jp1,its(n))  ! southern edge
              field(i,2,n)=     field(i,jp2,its(n))  ! southern edge
           enddo

           do i= 3,ip2
              j= jp2-i+3
              field(i,jp3,n)=  field(3,j,itn(n))     ! northern edge
              field(i,jp4,n)=  field(4,j,itn(n))     ! northern edge
           enddo
       endif
    enddo

endif

end subroutine create_halo


! ----------------------------------------------------------------------
subroutine find_crossing_rivers (ntl, id, jd, tcl, ncr, icr, jcr, tcr)
! ----------------------------------------------------------------------

implicit none

integer :: n, i, j, ktr, ip1, ip2, ip3, ip4, jp1, jp2, jp3, jp4

integer, intent(in)     :: ntl, id, jd
integer, intent(out)    :: ncr
integer, intent(out)    :: icr(:,:), jcr(:,:), tcr(:)

real, intent(in)        :: tcl(:,:,:)


ip1= id + 1  ;  jp1= jd + 1
ip2= id + 2  ;  jp2= jd + 2
ip3= id + 3  ;  jp3= jd + 3
ip4= id + 4  ;  jp4= jd + 4


! check for crossing of tocell paths
ktr= 0 ; icr= 0 ; jcr= 0 ; tcr= 0
if (ntl == 1) then
   do n= 1,ntl
      do j= 2,jp2
         do i= 2,ip2
            if ((tcl(i,j,n) == 128. .and. tcl(i,j+1,n) ==  2.) .or. &
                (tcl(i,j,n) ==  32. .and. tcl(i,j+1,n) ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j+1
                     tcr(ktr)= n
            endif
            if ((tcl(i,j,n) ==   2. .and. tcl(i+1,j,n) ==  8.) .or. &
                (tcl(i,j,n) == 128. .and. tcl(i+1,j,n) == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i+1
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j
                     tcr(ktr)= n
            endif
         enddo
      enddo
   enddo
else
   do n= 1,ntl
      do j= 3,jp1
         do i= 3,ip2
            if ((tcl(i,j,n) == 128. .and. tcl(i,j+1,n) ==  2.) .or. &
                (tcl(i,j,n) ==  32. .and. tcl(i,j+1,n) ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j+1
                     tcr(ktr)= n
            endif
         enddo
      enddo
      do j= 3,jp2
         do i= 3,ip1
            if ((tcl(i,j,n) ==   2. .and. tcl(i+1,j,n) ==  8.) .or. &
                (tcl(i,j,n) == 128. .and. tcl(i+1,j,n) == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i+1
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j
                     tcr(ktr)= n
            endif
         enddo
      enddo
   enddo

   do n= 1,ntl
      if (mod(n,2) == 0) then
          do i= 2,ip3
             if ((tcl(i,2,n)    ==   2. .and. tcl(i,3,n)    ==  2.) .or. &
                 (tcl(i,2,n)    == 128. .and. tcl(i,3,n)    ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i
                     jcr(1,ktr)= 2 ;  jcr(2,ktr)= 3
                     tcr(ktr)= n
             endif
             if ((tcl(i,jp2,n) == 128. .and. tcl(i,jp3,n) ==  2.) .or. &
                 (tcl(i,jp2,n) ==  32. .and. tcl(i,jp3,n) ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i   ;  icr(2,ktr)= i
                     jcr(1,ktr)= jp2 ;  jcr(2,ktr)= jp3
                     tcr(ktr)= n
             endif
          enddo
          do j= 2,jp3
             if ((tcl(2,j,n)    ==   2. .and. tcl(3,j,n)    ==  8.) .or. &
                 (tcl(2,j,n)    == 128. .and. tcl(3,j,n)    == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= 2 ;  icr(2,ktr)= 3
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j
                     tcr(ktr)= n
             endif
             if ((tcl(ip2,j,n) ==   2. .and. tcl(ip3,j,n) ==  2.) .or. &
                 (tcl(ip2,j,n) == 128. .and. tcl(ip3,j,n) ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= ip2 ;  icr(2,ktr)= ip3
                     jcr(1,ktr)= j   ;  jcr(2,ktr)= j
                     tcr(ktr)= n
             endif
          enddo
      else
          do i= 2,ip3
             if ((tcl(i,2,n)    == 128. .and. tcl(i,3,n)    ==  2.) .or. &
                 (tcl(i,2,n)    ==  32. .and. tcl(i,3,n)    ==  8.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i ;  icr(2,ktr)= i
                     jcr(1,ktr)= 2 ;  jcr(2,ktr)= 3
                     tcr(ktr)= n
             endif
             if ((tcl(i,jp2,n) == 128. .and. tcl(i,jp3,n) ==  8.) .or. &
                 (tcl(i,jp2,n) ==  32. .and. tcl(i,jp3,n) == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= i   ;  icr(2,ktr)= i
                     jcr(1,ktr)= jp2 ;  jcr(2,ktr)= jp3
                     tcr(ktr)= n
             endif
          enddo
          do j= 2,jp3
             if ((tcl(2,j,n)    == 128. .and. tcl(3,j,n)    ==  8.) .or. &
                 (tcl(2,j,n)    ==  32. .and. tcl(3,j,n)    == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= 2 ;  icr(2,ktr)= 3
                     jcr(1,ktr)= j ;  jcr(2,ktr)= j
                     tcr(ktr)= n
             endif
             if ((tcl(ip2,j,n) ==   2. .and. tcl(ip3,j,n) ==  8.) .or. &
                 (tcl(ip2,j,n) == 128. .and. tcl(ip3,j,n) == 32.)) then
                     ktr= ktr + 1
                     icr(1,ktr)= ip2 ;  icr(2,ktr)= ip3
                     jcr(1,ktr)= j   ;  jcr(2,ktr)= j
                     tcr(ktr)= n
             endif
          enddo
      endif
   enddo
endif
ncr= ktr

end subroutine find_crossing_rivers

end
