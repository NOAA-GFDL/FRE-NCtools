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
!* License along with FRE-NCTools (LICENSE.md) .  If not, see
!* <http://www.gnu.org/licenses/>.
!***********************************************************************
program cr_lake_files

use horiz_interp_mod

implicit none

integer, parameter :: maxdims= 3
integer, parameter :: idl= 2880, jdl= 1440
integer, parameter :: idlp1= idl+1, jdlp1= jdl+1
integer, parameter :: nvar_glcc= 2
integer, parameter :: iwbd= 1, ipwt= 2
integer, parameter :: ido= 360, jdo= 180
integer, parameter :: idop1= ido+1, jdop1= jdo+1
integer, parameter :: ncell_cnv= 8
integer, parameter :: ni_cells= 3, nj_cells= 3
integer, parameter :: ntilmx= 6

real, parameter :: lat1= -90., lat2= 90., lon1= 0., lon2= 360.
real, parameter :: erad= 6.371e+3
real, parameter :: mval_mdl= -9999.
real, parameter :: lake_depth_small= 2.
real, parameter :: lake_tau_large= 86400.*100.
real, parameter :: lake_tau_small= 1200.
real, parameter :: reso= 1.
real, parameter :: connect_to_next_0= 0.
!real, parameter :: travel_thresh= 3.
real, parameter :: max_slope_to_next = 0.1

include 'netcdf.inc'

integer :: i, j, n, l, rcode, varid, ncid, attnum, id, jd, idp1, jdp1
integer :: latid, lonid, latbid, lonbid, latdim, latbdim, londim
integer :: lonbdim, k, ktr, rcode2, ndm, latgid, longid, ntiles
integer :: varid2, varid3, varid4, varid5, varid6, varid7, i2, j2
integer :: ip, jp, fill_val, j1, i1, idx, jdx, ii, jj, idp2, jdp2
integer :: n1, np, varid8, varid9
real :: pi, dtr, sum, mval_frac, mval_tocell, mval_landf, mval_cella
real :: scale, sum1, sum2, asum1, asum2, wmax, wmin, mval_travel
real :: travel_thresh
character(len=8)   :: var
character(len=40)  :: var_units
character(len=100) :: fname, glcc_file
logical :: interp_glcc_to_1deg

integer, dimension (maxdims)           :: start, count, dimids, ndims
integer, dimension (idl,jdl)           :: idat
integer, dimension (ntilmx)            :: itw, ite, its, itn

real, dimension (idl)                  :: lonl, lonin
real, dimension (jdl)                  :: latl, latin
real, dimension (ido)                  :: lono
real, dimension (jdo)                  :: lato
real, dimension (idlp1)                :: lonlb
real, dimension (jdlp1)                :: latlb
real, dimension (idop1)                :: lonob
real, dimension (jdop1)                :: latob
real, dimension (idl,jdl)              :: arlat
real, dimension (ido,jdo)              :: arlato
real, dimension (idl,jdl,nvar_glcc)    :: wbdat
real, dimension (ido,jdo,nvar_glcc)    :: wbd_cnv


character(len=8), dimension (nvar_glcc) :: vname_glcc= &
(/ 'WaterBod', 'PWetland' /)

character(len=100), dimension (ntilmx) :: lname_glcc, river_input_file

real, allocatable, dimension (:)       :: lat_idx, lon_idx
real, allocatable, dimension (:,:)     :: interp_out, interp_mask, data_in, data_out
real, allocatable, dimension (:,:,:)   :: lat, lon, tocell, land_frac, cell_area, lake_depth, &
                                          lake_tau, wbd, pwt, cnct_next, whole_lake, travel, &
                                          max_slp2nxt

real*4, allocatable, dimension (:,:)     :: adat2

pi= 4.*atan(1.)
dtr= pi/180.

interp_glcc_to_1deg= .false.

!open (10, file= 'out.cr_lake_files', form= 'formatted')

read (5,*) ntiles
do n= 1,ntiles
   read (5,'(a)') river_input_file(n)
enddo
!read (5,'(l1)') interp_glcc_to_1deg
read (5,'(a)') glcc_file
read (5,*) travel_thresh
close (5)

!write (6,*) 'interp_glcc_to_1deg= ', interp_glcc_to_1deg


! ----------------------------------------------------------------------
!  read glcc file -- get lats, lons, waterbod, pwetland
! ----------------------------------------------------------------------
rcode= NF_OPEN (trim(glcc_file), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open glcc netcdf file"
    write (6,*) trim(glcc_file)
    stop 100
endif

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lat variable" ; stop 102
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= jdl) then
    write (6,*) "ERROR: inconsistent lat dimension, glcc, jdl= ", n
    stop 103
endif
start= 1 ;  count= 1 ;  count(1)= jdl
rcode= nf_get_vara_double (ncid, latid, start, count, latl)

if (latl(1) > -89.) then
!    write (6,*) "WARNING: glcc lats start at ", latl(1)
    if (latl(1) > 89.) then
!        write (6,*) "  -- flip lats to start at 90S -- "
        latin= latl
        do j= 1,jdl
           j1= jdlp1-j
           latl(j)= latin(j1)
        enddo
    else
        stop 2
    endif

endif

latlb(1)= lat1
latlb(jdlp1)= -(lat1)
do j= 2,jdl
   latlb(j)= 0.5*(latl(j)+latl(j-1))
enddo


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    write (6,*) "ERROR: cannot find glcc lon variable" ; stop 104
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), n)
if (n /= idl) then
    write (6,*) "ERROR: inconsistent lon dimension, glcc, idl= ", n
    stop 105
endif
start= 1 ;  count= 1 ;  count(1)= idl
rcode= nf_get_vara_double (ncid, lonid, start, count, lonl)

if (lonl(1) < -179.) then
!    write (6,*) "WARNING: glcc lons start at ", lonl(1)
!    write (6,*) "  -- shift lons to start at zero -- "
    lonin= lonl
    do i= 1,idl/2
       lonl(i)= lonin(i+idl/2)
    enddo
    do i= idl/2+1,idl
       lonl(i)= lon2 + lonin(i-idl/2)
    enddo
endif

lonlb(1)= lon1
lonlb(idlp1)= lon2
do i= 2,idl
   lonlb(i)= 0.5*(lonl(i)+lonl(i-1))
enddo

! compute area of latitude
arlat= 0.
do j= 1,jdl
   do i= 1,idl
      arlat(i,j)= erad*erad*(lonlb(i+1)-lonlb(i))*dtr* &
           (sin(latlb(j+1)*dtr) - sin(latlb(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdl
   do i= 1,idl
      sum= sum + arlat(i,j)
   enddo
enddo
!write (6,*) 'sum of glcc area= ', sum

!write (10,'(/a)') 'glcc lats'
!do j= 1,jdl
!   write (10,'(i6,5f10.3)') j, latl(j), latlb(j), latlb(j+1), arlat(1,j)
!enddo

!write (10,'(/a)') 'glcc lons'
!do i= 1,idl
!   write (10,'(i6,3f10.3)') i, lonl(i), lonlb(i), lonlb(i+1)
!enddo


wbdat= mval_mdl ;   lname_glcc=''
do n= 1,nvar_glcc
!   write (6,*) 'read glcc data: ', trim(vname_glcc(n))
   rcode= nf_inq_varid (ncid, trim(vname_glcc(n)), varid)
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find glcc variable ",  trim(vname_glcc(n))
       stop 110
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), i)
   rcode= nf_inq_dimlen (ncid, dimids(2), l)
   if (i /= idl .or. l /= jdl) then
       write (6,*) "ERROR: inconsistent WaterBod dimensions, glcc, id= ", i, ', jd= ', l
       stop 111
   endif
   start= 1 ;  count= 1 ;  count(1)= idl ;  count(2)= jdl
   rcode= nf_get_vara_int (ncid, varid, start, count, idat)

   fill_val= -999999
   rcode= nf_inq_attid (ncid, varid, '_FillValue', attnum)
   if (rcode == 0) rcode= nf_get_att_int (ncid, varid, '_FillValue', fill_val)
!   write (6,*) 'fill_val= ', fill_val

   scale= -99999999.
   rcode= nf_inq_attid (ncid, varid, 'scale_factor', attnum)
   if (rcode == 0) rcode= nf_get_att_double (ncid, varid, 'scale_factor', scale)
!   write (6,*) 'scale= ', scale

   rcode= nf_inq_attid (ncid, varid, 'long_name', attnum)
   if (rcode == 0) rcode= nf_get_att_text (ncid, varid, 'long_name', lname_glcc(n))
!   write (6,*) 'long_name= ', trim(lname_glcc(n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)
!   write (6,*) 'units= ', var_units

! flip lats and lons, scale by scale factor
   do j= 1,jdl
      j1= jdlp1-j
!      write (6,*) 'j= ', j, j1
      do i= 1,idl/2
         if (idat(i+idl/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(idat(i+idl/2,j1))*scale
         endif
      enddo
      do i= idl/2+1,idl
         if (idat(i-idl/2,j1) /= fill_val) then
             wbdat(i,j,n)= real(idat(i-idl/2,j1))*scale
         endif
      enddo
   enddo
enddo

rcode= nf_close (ncid)

do j= 1,jdl
   do i= 1,idl
      if (wbdat(i,j,1) == mval_mdl .or. wbdat(i,j,2) == mval_mdl) then
          write (6,*) "WARNING: wbdat has missing value, ", wbdat(i,j,1), wbdat(i,j,2)
      endif
   enddo
enddo

sum= 0.
do j= 1,jdl
   if (latl(j) >= 12. .and. latl(j) <= 14.) then
       do i= 1,idl
          if (lonl(i) >= 13. .and. lonl(i) <= 16.) then
              sum= sum + arlat(i,j)*wbdat(i,j,1)
          endif
       enddo
   endif
enddo

! ----------------------------------------------------------------------
! try averaging 0.125-degree data to 1-degree
! ----------------------------------------------------------------------
do j= 1,jdo
   lato(j)= lat1 + reso*0.5 + real(j-1)*reso
enddo

latob(1)= lat1
latob(jdop1)= -(lat1)
do j= 2,jdo
   latob(j)= 0.5*(lato(j)+lato(j-1))
enddo

do i= 1,ido
   lono(i)= lon1 + reso*0.5 + real(i-1)*reso
enddo

lonob(1)= lon1
lonob(idop1)= lon2
do i= 2,ido
   lonob(i)= 0.5*(lono(i)+lono(i-1))
enddo

! compute area of latitude
arlato= 0.
do j= 1,jdo
   do i= 1,ido
      arlato(i,j)= erad*erad*(lonob(i+1)-lonob(i))*dtr* &
           (sin(latob(j+1)*dtr) - sin(latob(j)*dtr))
   enddo
enddo

sum= 0.
do j= 1,jdo
   do i= 1,ido
      sum= sum + arlato(i,j)
   enddo
enddo
!write (6,*) 'sum of 1-degree area= ', sum

!write (10,'(/a)') '1-degree lats'
!do j= 1,jdo
!   write (10,'(i6,5f10.3)') j, lato(j), latob(j), latob(j+1), arlato(1,j)
!enddo

!write (10,'(/a)') '1-degree lons'
!do i= 1,ido
!   write (10,'(i6,3f10.3)') i, lono(i), lonob(i), lonob(i+1)
!enddo

if ( interp_glcc_to_1deg ) then
    wbd_cnv= mval_mdl ; sum= 0.
    do j= 1,jdo
       j1= (j-1)*ncell_cnv+1
       j2= j*ncell_cnv
       do i= 1,ido
          i1= (i-1)*ncell_cnv+1
          i2= i*ncell_cnv
          asum1= 0. ;  asum2= 0.
          sum1= 0. ;  sum2= 0.
          do jj= j1,j2
             do ii= i1,i2
                if (wbdat(ii,jj,1) /= mval_mdl) then
                    sum= sum + arlat(ii,jj)
                    asum1= asum1 + arlat(ii,jj)
                    sum1= sum1 + wbdat(ii,jj,1)*arlat(ii,jj)
                endif
                if (wbdat(ii,jj,2) /= mval_mdl) then
                    asum2= asum2 + arlat(ii,jj)
                    sum2= sum2 + wbdat(ii,jj,2)*arlat(ii,jj)
                endif
             enddo
          enddo
          if (abs(asum1 - arlato(i,j)) < 1.0e-10) then
              wbd_cnv(i,j,1)= sum1/asum1
          else
              write (6,*) "ERROR: wbd 1-deg areas do not agree, ", asum1, arlato(i,j)
              stop 17
          endif
          if (abs(asum2 - arlato(i,j)) < 1.0e-10) then
              wbd_cnv(i,j,2)= sum2/asum2
          else
              write (6,*) "ERROR: pwt 1-deg areas do not agree, ", asum2, arlato(i,j)
              stop 17
          endif
       enddo
    enddo
    write (6,*) 'sum of conversion area= ', sum
endif

! ----------------------------------------------------------------------
!  get lon and lat dims from first river file -- should be identical
!    for all files
! ----------------------------------------------------------------------

rcode= NF_OPEN (trim(river_input_file(1)), NF_NOWRITE, ncid)
if (rcode /= 0) then
    write (6,*) "ERROR: cannot open river netcdf file"
    write (6,*) trim(river_input_file(1))
    stop 1
endif

rcode= nf_inq_varid (ncid, 'lat', latid)         ! number of lats
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_y', latid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lat variable" ; stop 2
    endif
endif
rcode= nf_inq_vardimid (ncid, latid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), jd)
!write (6,*) 'jd= ', jd
jdp1= jd + 1
jdp2= jd + 2

allocate (lat_idx(jd))
start= 1 ;  count= 1 ;  count(1)= jd
rcode= nf_get_vara_double (ncid, latid, start, count, lat_idx)


rcode= nf_inq_varid (ncid, 'lon', lonid)         ! number of lons
if (rcode /= 0) then
    rcode2 = nf_inq_varid (ncid, 'grid_x', lonid)
    if (rcode2 /= 0) then
        write (6,*) "ERROR: cannot find lon variable" ; stop 3
    endif
endif
rcode= nf_inq_vardimid (ncid, lonid, dimids)
rcode= nf_inq_dimlen (ncid, dimids(1), id)
!write (6,*) 'id= ', id
idp1= id + 1
idp2= id + 2

allocate (lon_idx(id))
start= 1 ;  count(1)= id
rcode= nf_get_vara_double (ncid, lonid, start, count, lon_idx)

rcode= nf_close (ncid)



! ----------------------------------------------------------------------
! now open river files -- read lat,lon grids, tocell, land_frac,
!   cellarea, and travel
! ----------------------------------------------------------------------

allocate (lat(idp2,jdp2,ntiles), lon(idp2,jdp2,ntiles), travel(idp2,jdp2,ntiles))
allocate (tocell(idp2,jdp2,ntiles), land_frac(idp2,jdp2,ntiles), cell_area(idp2,jdp2,ntiles))

do n= 1,ntiles
!   write (6,'(i6,3x,a)')  n, trim(river_input_file(n))
!   write (10,'(/i6,3x,a)') n, trim(river_input_file(n))

   rcode= NF_OPEN (trim(river_input_file(n)), NF_NOWRITE, ncid)
   if  (rcode /= 0) then
       write (6,*) "ERROR: cannot open netcdf file"  ; stop 12
   endif

   start= 1 ;  count= 1

! ----------------------------------------------------------------------
! get latitudes and longitudes of input river file
! ----------------------------------------------------------------------

   rcode= nf_inq_varid (ncid, 'y', latid)         ! lat field
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


   rcode= nf_inq_varid (ncid, 'x', lonid)         ! lon field
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
       write (6,*) "ERROR: inconsistent lat dimension, ", ndm, jd ;  stop 35
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, lonid, start, count, lon(2:idp1,2:jdp1,n))


!   write (6,*) 'read tocell'
   rcode= nf_inq_varid (ncid, 'tocell', varid)     ! tocell field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find tocell" ; stop 40
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent tocell dimension, ", ndm, id ;  stop 45
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent tocell dimension, ", ndm, jd ;  stop 45
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, tocell(2:idp1,2:jdp1,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_tocell= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_tocell)
   endif

   where (tocell(:,:,n) == mval_tocell) tocell(:,:,n)= mval_mdl


!   write (6,*) 'read land_frac'
   rcode= nf_inq_varid (ncid, 'land_frac', varid)     ! land fraction field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find land_frac" ; stop 50
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent land_frac dimension, ", ndm, id ;  stop 55
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent land_frac dimension, ", ndm, jd ;  stop 55
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, land_frac(2:idp1,2:jdp1,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_landf= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_landf)
   endif

   where (land_frac(:,:,n) == mval_landf) land_frac(:,:,n)= mval_mdl


!   write (6,*) 'read cellarea'
   rcode= nf_inq_varid (ncid, 'cellarea', varid)     ! cell area field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find cellarea" ; stop 60
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent cellarea dimension, ", ndm, id ;  stop 65
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent cellarea dimension, ", ndm, jd ;  stop 65
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, cell_area(2:idp1,2:jdp1,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_cella= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_cella)
   endif

   where (cell_area(:,:,n) == mval_cella) cell_area(:,:,n)= 0.


!   write (6,*) 'read travel'
   rcode= nf_inq_varid (ncid, 'travel', varid)     ! travel field
   if (rcode /= 0) then
       write (6,*) "ERROR: cannot find travel" ; stop 70
   endif
   rcode= nf_inq_vardimid (ncid, varid, dimids)
   rcode= nf_inq_dimlen (ncid, dimids(1), ndm)
   if (ndm /= id) then
       write (6,*) "ERROR: inconsistent travel dimension, ", ndm, id ;  stop 72
   endif
   rcode= nf_inq_dimlen (ncid, dimids(2), ndm)
   if (ndm /= jd) then
       write (6,*) "ERROR: inconsistent travel dimension, ", ndm, jd ;  stop 72
   endif

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= nf_get_vara_double (ncid, varid, start, count, travel(2:idp1,2:jdp1,n))

   var_units= ' '
   rcode= nf_get_att_text (ncid, varid, "units", var_units)

   mval_travel= 1.e+20
   rcode= nf_inq_attid (ncid, varid, 'missing_value', attnum)
   if (rcode == 0) then
       rcode= nf_get_att_double (ncid, varid, 'missing_value', mval_travel)
   endif

   where (travel(:,:,n) == mval_travel) travel(:,:,n)= mval_mdl

   rcode= nf_close (ncid)

enddo


! ----------------------------------------------------------------------
! interpolate glcc water bodies and permanent wetlands to model grid
! ----------------------------------------------------------------------

allocate (wbd(idp2,jdp2,ntiles), pwt(idp2,jdp2,ntiles), interp_out(id,jd), data_out(id,jd))

wbd= mval_mdl
pwt= mval_mdl

call horiz_interp_init

 if (interp_glcc_to_1deg) then
    allocate (interp_mask(ido,jdo))
    allocate (data_in(ido,jdo))
    interp_mask= 1.
    where (wbd_cnv(:,:,1) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
!       write (6,*) 'WaterBod, tile= ', n
       data_in= wbd_cnv(:,:,1)
       call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       wbd(2:idp1,2:jdp1,n)= data_out
       where (interp_out == 0.) wbd(2:idp1,2:jdp1,n)= mval_mdl
    enddo

    interp_mask= 1.
    where (wbd_cnv(:,:,2) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
!       write (6,*) 'PWetland, tile= ', n
       data_in= wbd_cnv(:,:,2)
       call horiz_interp (data_in, lonob*dtr, latob*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       pwt(2:idp1,2:jdp1,n)= data_out
       where (interp_out == 0.) pwt(2:idp1,2:jdp1,n)= mval_mdl
    enddo

else
    allocate (interp_mask(idl,jdl))
    allocate (data_in(idl,jdl))
    interp_mask= 1.
    where (wbdat(:,:,1) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
!       write (6,*) 'WaterBod, tile= ', n
       data_in= wbdat(:,:,1)
       call horiz_interp (data_in, lonlb*dtr, latlb*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       wbd(2:idp1,2:jdp1,n)= data_out
       where (interp_out == 0.) wbd(2:idp1,2:jdp1,n)= mval_mdl
    enddo

    interp_mask= 1.
    where (wbdat(:,:,2) == mval_mdl) interp_mask= 0.

    do n= 1,ntiles
!       write (6,*) 'PWetland, tile= ', n
       data_in= wbdat(:,:,2)
       call horiz_interp (data_in, lonlb*dtr, latlb*dtr, lon(2:idp1,2:jdp1,n)*dtr, &
            lat(2:idp1,2:jdp1,n)*dtr, data_out, verbose=0, mask_in=interp_mask, &
            mask_out=interp_out, interp_method="bilinear")
       pwt(2:idp1,2:jdp1,n)= data_out
       where (interp_out == 0.) pwt(2:idp1,2:jdp1,n)= mval_mdl
    enddo
endif

deallocate (interp_out, interp_mask)
deallocate (data_in, data_out)


allocate (lake_depth(idp2,jdp2,ntiles), lake_tau(idp2,jdp2,ntiles))
allocate (whole_lake(idp2,jdp2,ntiles), cnct_next(idp2,jdp2,ntiles))
allocate (max_slp2nxt(idp2,jdp2,ntiles))

lake_depth= lake_depth_small
lake_tau= 0.
cnct_next= connect_to_next_0
max_slp2nxt= max_slope_to_next

! define whole_lake_area field and lake depth field
do n= 1,ntiles
   do j= 1,jdp2
      do i= 1,idp2
         whole_lake(i,j,n)= wbd(i,j,n)*cell_area(i,j,n)
      enddo
   enddo
enddo

! set ocean values to missing
where (land_frac == 0.)
   wbd= mval_mdl
   pwt= mval_mdl
   lake_depth= mval_mdl
   lake_tau= mval_mdl
   cnct_next= mval_mdl
   whole_lake= mval_mdl
   max_slp2nxt= mval_mdl
endwhere

! set part-land values to zero
where (land_frac > 0 .and. land_frac < 1.)
   wbd= 0.
   pwt= 0.
   whole_lake= 0.
endwhere


! set up the 'halo' for each tile
!   get edge data for tocell, lat, and lon

if (ntiles == 1) then
    itw(n)= 1 ;  ite(n)= 1 ;  its(n)= 1 ;  itn(n)= 1
    call create_halo (ntiles, id, jd, itw, ite, its, itn, tocell)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
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
    call create_halo (ntiles, id, jd, itw, ite, its, itn, travel)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lat)
    call create_halo (ntiles, id, jd, itw, ite, its, itn, lon)
endif

! set lake value to zero at grid cells in antarctic, arctic, and greenland, and also where travel < 2

do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (tocell(i,j,n) == mval_mdl .or. wbd(i,j,n) == mval_mdl) go to 109
         if (lat(i,j,n) < -60.) then                    ! antarctica
             wbd(i,j,n)= 0.
         else if (lat(i,j,n) > 70. ) then               ! arctic
             wbd(i,j,n)= 0.
         else if (lat(i,j,n) > 58. .and. &              ! rest of greenland
             (lon(i,j,n) > 302. .and. lon(i,j,n) < 348.)) then
             wbd(i,j,n)= 0.
         else
             if (wbd(i,j,n) > 0.1 .and. travel(i,j,n) <= travel_thresh) then
                 wbd(i,j,n)= 0.
             endif
         endif
109      continue
      enddo
   enddo
enddo

where (wbd == 0.)
   whole_lake= 0.
endwhere


wmax= -9999999. ;  wmin= 9999999.
do n= 1,ntiles
   do j= 2,jdp1
      do i= 2,idp1
         if (wbd(i,j,n) /= mval_mdl) then
             wmax= max(wmax,wbd(i,j,n))
             wmin= min(wmin,wbd(i,j,n))
         endif
      enddo
   enddo
enddo
!write (6,*) 'wmax= ', wmax, ', wmin= ', wmin



! ----------------------------------------------------------------------
!  create new netCDF data set; overwrite existing dataset
! ----------------------------------------------------------------------

do n= 1,ntiles
   write (fname, '(a,i1,a)') 'lake_frac.tile', n, '.nc'
   rcode= NF_CREATE (trim(fname), NF_CLOBBER, ncid)
   rcode= NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'filename', len_trim(fname), trim(fname))

! ----------------------------------------------------------------------
!  create dimensions, coordinate variables, coordinate attributes for
!    mean files
! ----------------------------------------------------------------------

!  create dimensions (grid_x, grid_y)
   rcode= NF_DEF_DIM (ncid, 'grid_x',  id,   londim)
   rcode= NF_DEF_DIM (ncid, 'grid_y',  jd,   latdim)

!  create coordinate variables
   rcode= NF_DEF_VAR (ncid, 'grid_x',  NF_DOUBLE, 1, (/ londim /),  lonid)
   rcode= NF_DEF_VAR (ncid, 'grid_y',  NF_DOUBLE, 1, (/ latdim /),  latid)

!  create attributes for coordinate variables
!    longitude:
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'long_name', 16, 'T-cell longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'axis', 1, 'X')
   rcode= NF_PUT_ATT_TEXT (ncid, lonid, 'units', 9, 'degrees_E')

!    latitude:
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'long_name', 15, 'T-cell latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'axis', 1, 'Y')
   rcode= NF_PUT_ATT_TEXT (ncid, latid, 'units', 9, 'degrees_N')

!    create data variable and attributes
   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'x', NF_DOUBLE, 2, ndims, longid)
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'long_name', 20, 'Geographic longitude')
   rcode= NF_PUT_ATT_TEXT (ncid, longid, 'units', 9, 'degrees_E')

   rcode= NF_DEF_VAR (ncid, 'y', NF_DOUBLE, 2, ndims, latgid)
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'long_name', 19, 'Geographic latitude')
   rcode= NF_PUT_ATT_TEXT (ncid, latgid, 'units', 9, 'degrees_N')

   ndims(1)= londim ;  ndims(2)= latdim
   rcode= NF_DEF_VAR (ncid, 'lake_frac', NF_DOUBLE, 2, ndims, varid)
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'long_name', 13, 'lake_fraction')
   rcode= NF_PUT_ATT_TEXT (ncid, varid, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'lake_depth_sill', NF_DOUBLE, 2, ndims, varid3)
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'long_name', 15, 'lake_depth_sill')
   rcode= NF_PUT_ATT_TEXT (ncid, varid3, 'units', 1, 'm')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid3, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'lake_tau', NF_DOUBLE, 2, ndims, varid4)
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'long_name', 8, 'lake_tau')
   rcode= NF_PUT_ATT_TEXT (ncid, varid4, 'units', 1, 's')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid4, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, trim(vname_glcc(1)), NF_DOUBLE, 2, ndims, varid5)
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'long_name', len_trim(lname_glcc(1)), &
          trim(lname_glcc(1)))
   rcode= NF_PUT_ATT_TEXT (ncid, varid5, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid5, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, trim(vname_glcc(2)), NF_DOUBLE, 2, ndims, varid6)
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'long_name', len_trim(lname_glcc(2)), &
          trim(lname_glcc(2)))
   rcode= NF_PUT_ATT_TEXT (ncid, varid6, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid6, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'connected_to_next', NF_DOUBLE, 2, ndims, varid7)
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'long_name', 20, 'lake connection flag')
   rcode= NF_PUT_ATT_TEXT (ncid, varid7, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid7, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'whole_lake_area', NF_DOUBLE, 2, ndims, varid8)
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'long_name', 18, 'total area of lake')
   rcode= NF_PUT_ATT_TEXT (ncid, varid8, 'units', 2, 'm2')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid8, 'missing_value', NF_DOUBLE, 1, mval_mdl)

   rcode= NF_DEF_VAR (ncid, 'max_slope_to_next', NF_DOUBLE, 2, ndims, varid9)
   rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'long_name', 26, 'max value of slope_to_next')
   rcode= NF_PUT_ATT_TEXT (ncid, varid9, 'units', 4, 'none')
   rcode= NF_PUT_ATT_DOUBLE (ncid, varid9, 'missing_value', NF_DOUBLE, 1, mval_mdl)

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

!    lake fraction data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid, start, count, wbd(2:idp1,2:jdp1,n))

!    lake depth data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid3, start, count, lake_depth(2:idp1,2:jdp1,n))

!    lake tau data
   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid4, start, count, lake_tau(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid5, start, count, wbd(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid6, start, count, pwt(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid7, start, count, cnct_next(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid8, start, count, whole_lake(2:idp1,2:jdp1,n))

   start= 1 ;  count(1)= id ;  count(2)= jd
   rcode= NF_PUT_VARA_DOUBLE (ncid, varid9, start, count, max_slp2nxt(2:idp1,2:jdp1,n))

!  close netcdf file
   rcode= NF_CLOSE (ncid)
enddo

!close (10)

deallocate (lat_idx, lon_idx)
deallocate (lat, lon, tocell, land_frac, cell_area)
deallocate (lake_depth, lake_tau, max_slp2nxt)
deallocate (wbd, pwt, cnct_next, whole_lake)


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


end
