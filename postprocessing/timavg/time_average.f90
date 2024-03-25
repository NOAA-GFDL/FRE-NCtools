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
!-----------------------------------------------------------------------
! Copyright (C) 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
!
!  This program averages variables stored in netCDF format over the time axis
!-----------------------------------------------------------------------

program time_average

use  netcdf
implicit none

integer, parameter  :: MAX_FILES = 200
integer, parameter  :: MAX_ATT_NAME = 256
integer, parameter  :: header_buffer_val = 16384
integer, parameter  :: BSIZE = 65536
integer            ::  blksz = BSIZE  ! blksz must be writable for nf__create

!-----------------------------------------------------------------------
character(len=2048) :: file_names(MAX_FILES), file_name_out
logical :: use_end_time = .true.
logical :: verbose = .false.
logical :: add_cell_methods = .false.
logical :: skip_tavg_errors = .false.
logical :: suppress_warnings = .false.
real    :: frac_valid_data = 0.
integer :: user_deflation = -1, user_shuffle = -1   ! user requests, -1 means use defaults
integer :: deflate_in, deflation_in, shuffle_in, deflation_out, shuffle_out
integer :: max_input_deflation = 0, max_input_shuffle = 0   ! to use for new variables

namelist /input/   file_names, file_name_out, use_end_time, verbose, &
                   add_cell_methods, skip_tavg_errors, frac_valid_data, &
                   suppress_warnings, user_deflation, user_shuffle
!-----------------------------------------------------------------------
! data structure for storing information about individual variables
 type variable_type
   character(len=NF90_MAX_NAME) :: name, avg_info
   integer :: ndim, natt, xtype, varid_out, shape(4), axes(NF90_MAX_VAR_DIMS)
   logical :: static, do_missval, do_avg, use_data_wt, skip, uses_tod
   real(8)        :: tavg(3)
   real           :: missval, data_wt_max
   real, pointer  :: data_sum(:,:,:,:), data_wt(:,:,:,:)
 end type

   type(variable_type), allocatable :: Vars(:)

   ! default precision (32/64-bit) real variables
   real, allocatable :: ddata(:,:,:,:)
   real    :: half = 0.50
   real    :: small = 1.e-6
   real    :: data_wt_min
   ! double precision (64-bit) time variables
   real(8) :: tavg(3), time
   real(8) :: time_bounds(3), last_time, time_out, climo_time

   integer :: ncid_in, ncid_out, numvars, numdims, numatts
   integer :: istat, ifile, itime, ivar, idim, nvar, ndim, natt, dimlen, varid,      &
              i, ix, jy, kz, numtime, recdim, varid_recdim, varid_in, attnum, dimid, &
              nvar_out, nv_dimid, nc, nca, ncu, dimid_tod, mt, nv_varid
   integer :: max_shape(4), shape(4), tid(0:3), start(5)
   integer :: axes(NF90_MAX_VAR_DIMS), dimids_out(NF90_MAX_VAR_DIMS)
   logical :: do_avg, first, time_bounds_present, do_time_bounds, add_time_bounds, &
              do_climo_avg, change_bounds_to_climo, climo_bounds_present, do_period_out

   character(len=NF90_MAX_NAME) :: name, name_recdim, attname, version, name_time_bounds
   character(len=NF90_MAX_NAME) :: period_in, period_out
   character(len=NF90_MAX_NAME) :: time_att_units, time_units
   character(len=NF90_MAX_NAME) :: cell_methods, cell_methods_tavg, name_bnds
   integer :: time_bnds_id
   integer :: in_format, cmode

!-----------------------------------------------------------------------

  ! initialize blank input strings
    do ifile = 1,MAX_FILES
       do i=1,len(file_names); file_names(ifile)(i:i) = ' '; enddo
    enddo
       do i=1,len(file_name_out); file_name_out(i:i) = ' '; enddo

  ! create version string (may replace with CVS $Id: time_average.f90,v 20.0 2013/12/14 00:30:08 fms Exp $)
    version = 'FMS time averaging, version 3.0'
    if (precision(ddata) == precision(time)) then
        version = trim(version)//', precision=double'
    else
        version = trim(version)//', precision=float'
    endif

  ! read namelist
    read  (*,input,end=1,err=1)
 1  continue
    if (verbose) then
       ! write (*,input)
       write (*,'(a,a)') 'NetCDF Version = ', NF90_INQ_LIBVERS()
    endif


  ! error checks
    if (file_names(1)(1:1) == ' ')  &
              call error_handler ('No input file names specified')
    if (file_name_out(1:1) == ' ')  &
              call error_handler ('No output file name specified')
    if (user_deflation < -1 .or. user_deflation > 9) &
              call error_handler ('Deflation level must be 0 to 9')
    if (user_shuffle < -1 .or. user_shuffle > 1) &
              call error_handler ('Shuffle must be 0 or 1')

!-----------------------------------------------------------------------
!----------------- Loop through input files ----------------------------
!-----------------------------------------------------------------------

 do ifile = 1, MAX_FILES
      if ( file_names(ifile)(1:1) == ' ' ) exit
      if (verbose) print *, 'Processing file ',trim(file_names(ifile))

    !--- opening input files ---
     blksz = BSIZE  ! reset since this may have been altered by nf90_open
     istat = NF90_OPEN (trim(file_names(ifile)), NF90_NOWRITE, ncid_in, chunksize=blksz)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

    !--- get global info for this file ---
     istat = NF90_INQUIRE (ncid_in, nDimensions=ndim, nVariables=nvar, nAttributes=natt, unlimitedDimId=recdim, formatNum=in_format)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
    !--- get information related to the record (time) dimension ----
     istat = NF90_INQUIRE_DIMENSION (ncid_in, recdim, name_recdim, numtime)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_INQ_VARID (ncid_in, name_recdim, varid_recdim)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

     if (in_format == NF90_FORMAT_64BIT .or. in_format == NF90_FORMAT_CLASSIC) then
         if (user_deflation > 0) then
             print *, 'Not using compression on NetCDF3 file'
             user_deflation = 0
         endif
     else
         if (user_deflation == -1) then
             print *, 'Will use input file settings for NetCDF4 deflation level'
         else
             print *, 'Will use NetCDF4 deflation level=', user_deflation
         endif
         if (user_shuffle == -1) then
             print *, 'Will use input file settings for NetCDF4 shuffle'
         else
             print *, 'Will use NetCDF4 shuffle=', user_shuffle
         endif
    endif

     if (verbose) then
         write (*,*) '    ncid=', ncid_in
         write (*,101) ndim, nvar, natt, trim(name_recdim), recdim, numtime
     101 format (4x,'ndim=',i3,4x,'nvar=',i3,4x,'natt=',i2, &
               /,4x,'recdim: name=',a,4x,'id=',i2,4x,'length=',i5)
     endif

     ! get units of record dimension
     do i=1,NF90_MAX_NAME; time_att_units(i:i) = ' '; enddo
     do i=1,NF90_MAX_NAME; time_units    (i:i) = ' '; enddo
     istat = NF90_GET_ATT (ncid_in, varid_recdim, 'units', time_att_units)
     ! extract time units from time axis units attribute
     nca = len_trim(time_att_units)
     ncu = index(time_att_units(1:nca),' ')-1  ! use first word only
     time_units = time_att_units(1:ncu)

  if (ifile == 1) then
     ! get time axis bounds attribute if it is present
     do i=1,NF90_MAX_NAME; name_time_bounds(i:i) = ' '; enddo
     istat = NF90_GET_ATT (ncid_in, varid_recdim, 'bounds', name_time_bounds)
     time_bounds_present = .false.
     if (name_time_bounds(1:1) .ne. ' ') then
         istat = NF90_INQ_VARID (ncid_in, trim(name_time_bounds), varid)
         if (istat == NF90_NOERR) time_bounds_present = .true.
     endif
     do_time_bounds = time_bounds_present
     ! or check for climatological bounds attribute
     if (.not.time_bounds_present) then
        do i=1,NF90_MAX_NAME; name_time_bounds(i:i) = ' '; enddo
        istat = NF90_GET_ATT (ncid_in, varid_recdim, 'climatology', name_time_bounds)
        climo_bounds_present = name_time_bounds(1:1) .ne. ' '
        do_time_bounds = climo_bounds_present
     endif
     ! replace null character at end of string with space (bug fix)
     nc = len_trim(name_time_bounds)
     if (iachar(name_time_bounds(nc:nc))==0) name_time_bounds(nc:nc)=' '
     if (verbose) then
         print *, 'time_bounds_present=',time_bounds_present
         print *, 'climo_bounds_present=',climo_bounds_present
         print *, 'name_time_bounds=',trim(name_time_bounds)
     endif
  endif

    !--- contiguous time average or climatological time average? ---
     if (ifile == 1) then
         if (numtime > 1) then
             call check_for_climo_avg ( ncid_in, name_recdim, time_units, numtime, &
                                        do_climo_avg, period_in, period_out,       &
                                        name_time_bounds)
             if (do_climo_avg) climo_time = 0.0
         else
             do_climo_avg = .false.
         endif
         ! is the time bounds attribute changing from "bounds" to "climatology"?
         if (time_bounds_present .and. do_climo_avg) then
             change_bounds_to_climo = .true.
         else
             change_bounds_to_climo = .false.
         endif
     endif

!print *, 'do_climo_avg,change_bounds_to_climo=',do_climo_avg,change_bounds_to_climo
!print *, 'add_cell_methods,time_bounds_present=',add_cell_methods,time_bounds_present

    !--- add CF convention for time bounds and cell methods? ----
     if (add_cell_methods .and. .not.time_bounds_present .and. .not.climo_bounds_present) then
         do_time_bounds  = .true.
         add_time_bounds = .true.
         if (verbose) then
             write (*,*) 'Adding CF convention for time averaging.'
         endif
         nvar_out = nvar+1
         ! note: will add axis called "nv", this axis can not exist
         istat = NF90_INQ_DIMID (ncid_in, 'nv', dimid)
         if (istat == NF90_NOERR) call error_handler &
                        ('Dimension "nv" already exists, can not create')
     else
         add_time_bounds = .false.
         nvar_out = nvar
     endif

    !--- time bnds (assume the same name use by all variables) ---
    !--- should work for FMS generated files, but need to make this more robust ---
     do_avg = .true.
     istat = NF90_INQ_VARID (ncid_in, trim(name_time_bounds), time_bnds_id)
     if (istat /= NF90_NOERR) then
        do_avg = .false.
     endif

    !--- allocate array space when doing first file ---
     if (ifile == 1) then
         allocate (Vars(nvar_out))
         numdims = ndim
         numvars = nvar
         numatts = natt
         Vars(:)%use_data_wt = .false.
         Vars(:)%skip        = .false.
         Vars(:)%do_avg      = .false.
         Vars(:)%do_missval  = .false.
         Vars(:)%uses_tod    = .false.
         max_shape = 1
         do i=1,NF90_MAX_NAME
           Vars(:)%name(i:i)=' '
           Vars(:)%avg_info(i:i)=' '
         enddo

!print *, trim(Vars(nvar_out)%name),': use_data_wt,do_missval=',Vars(nvar_out)%use_data_wt,Vars(nvar_out)%do_missval

         istat = NF90_INQ_DIMID(ncid_in, "time_of_day_24", dimid_tod)
         if (istat /= NF90_NOERR) dimid_tod = -1

        !--- compute variable/axis lengths ---
        !--- loop thru variables (in first file) ---
         do ivar = 1, nvar
            istat = NF90_INQUIRE_VARIABLE (ncid_in, ivar, name=Vars(ivar)%name, &
                                           xtype=Vars(ivar)%xtype, ndims=Vars(ivar)%ndim, &
                                           dimids=Vars(ivar)%axes, natts=Vars(ivar)%natt)
            if (istat /= NF90_NOERR) call error_handler (ncode=istat)

            ! Check if var uses new time_of_day axis.
            do idim = 1,Vars(ivar)%ndim
               if (Vars(ivar)%axes(idim) .eq. dimid_tod) then
                  Vars(ivar)%uses_tod = .true.
               endif
            enddo

            call init_variable_type  ( ncid_in, ivar, Vars(ivar) )
            call alloc_variable_type ( numtime, Vars(ivar) )
            if (verbose) then
               write (*,33) trim(Vars(ivar)%name),Vars(ivar)%skip,Vars(ivar)%shape
            33 format ('name=',a,2x,'skip=',l5,2x,'shape=',4i3)
            endif
         enddo

         ! add time bounds variable
         if (add_time_bounds) then
            ivar = nvar+1
            if (do_climo_avg) then
                name_time_bounds = 'climatology_bounds'
            else
                name_time_bounds = trim(name_recdim)//'_bounds'
            endif
            Vars(ivar)%name  = name_time_bounds
!print *, trim(Vars(ivar)%name),': use_data_wt,do_missval=',Vars(ivar)%use_data_wt,Vars(ivar)%do_missval
            Vars(ivar)%xtype = NF90_REAL8
            Vars(ivar)%ndim  = 1
            Vars(ivar)%axes(1:2) = (/numdims+1,recdim/)
            Vars(ivar)%natt  = 0
            Vars(ivar)%static = .false.
            Vars(ivar)%shape  = (/2,1,1,1/)
            call alloc_variable_type ( numtime, Vars(ivar) )
         endif
        !--- done variable loop (for first file) ----
        !--------------------------------------------

        !--- array for reading/writing of netcdf fields ---
         allocate (ddata(max_shape(1),max_shape(2),max_shape(3),max_shape(4)))

     else
        !--- check subsequent files to make sure that number of variables is the same ----
        !--- will check variable names below ---
         if (nvar /= numvars) call error_handler ('Input netcdf files &
                            &do not have the same number of variables')
     endif

     if (verbose) then
         write (*,*) '  time_bounds_present=',time_bounds_present
         write (*,*) '  add_time_bounds=',add_time_bounds
         write (*,*) '  do_climo_avg=',do_climo_avg
         write (*,*) '  name_time_bounds=',trim(name_time_bounds)
     endif

   !-------------------------------------------------------------
   !---------------- OUTPUT FILE INITIALIZATION -----------------

     if (ifile == 1) then
       if (in_format==NF90_FORMAT_NETCDF4) then
          cmode = NF90_NETCDF4
       elseif (in_format==NF90_FORMAT_NETCDF4_CLASSIC) then
          cmode=IOR(NF90_NETCDF4,NF90_CLASSIC_MODEL)
       elseif (in_format==NF90_FORMAT_64BIT) then
          cmode=IOR(NF90_CLOBBER,NF90_64BIT_OFFSET)
       elseif (in_format==NF90_FORMAT_CLASSIC) then
          cmode=IOR(NF90_CLOBBER,NF90_CLASSIC_MODEL)
       else
          call error_handler ('Unknown netCDF format')
       endif
       ! create file
         blksz = BSIZE  ! reset since this may have been altered by nf90_open
         istat = NF90_CREATE (trim(file_name_out), cmode, ncid_out, chunksize=blksz)
         if (istat /= NF90_NOERR) call error_handler ('creating output file', ncode=istat)

       ! copy global attributes
         do attnum = 1, numatts
            ! get name
            istat = NF90_INQ_ATTNAME (ncid_in, NF90_GLOBAL, attnum, attname)
            if (istat /= NF90_NOERR) call error_handler (ncode=istat)
            ! change global filename attribute
            if (trim(attname) == 'filename') then
                istat = NF90_PUT_ATT (ncid_out, NF90_GLOBAL, 'filename', trim(file_name_out))
                if (istat /= NF90_NOERR) call error_handler (ncode=istat)
            else if (trim(attname) == 'time_averaging') then
                cycle ! skip (create new attribute later)
            else if (trim(attname) == 'comment') then
                cycle ! skip (create new attribute later)
            else
                ! copy all others
                istat = NF90_COPY_ATT (ncid_in, NF90_GLOBAL, attname, ncid_out, NF90_GLOBAL)
                if (istat /= NF90_NOERR) call error_handler ('copying global attributes', &
                                                              ncode=istat)
            endif
         enddo
         ! add attribute
         istat = NF90_PUT_ATT (ncid_out, NF90_GLOBAL, 'comment', trim(version))
         if (istat /= NF90_NOERR) call error_handler (ncode=istat)

       ! copy dimensions
         do idim = 1, numdims
            istat = NF90_INQUIRE_DIMENSION (ncid_in, idim, name, dimlen)
            if (istat /= NF90_NOERR) call error_handler (ncode=istat)
            if (idim == recdim) then
                istat = NF90_DEF_DIM (ncid_out, name, NF90_UNLIMITED, dimid)
                if (istat /= NF90_NOERR) call error_handler (ncode=istat)
            else
                istat = NF90_DEF_DIM (ncid_out, name, dimlen, dimid)
                if (istat /= NF90_NOERR) call error_handler (ncode=istat)
            endif
            ! mapping between input and output dimension ids
            dimids_out(idim) = dimid
         enddo

       ! add axis for time bounds variable
         if (add_time_bounds) then
            ! Add the bounds vertex bounds for the time bounds
             istat = NF90_DEF_DIM (ncid_out, 'nv', 2, nv_dimid)
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             dimids_out(numdims+1) = nv_dimid

            ! Add the vertex number variable to the netCDF file
             istat = NF90_DEF_VAR(ncid_out, 'nv', nf90_int, nv_dimid, nv_varid)
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             istat = NF90_PUT_ATT(ncid_out, nv_varid, 'long_name', 'vertex number')
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             istat = NF90_PUT_ATT(ncid_out, nv_varid, 'units', 'none')
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             istat = NF90_PUT_ATT(ncid_out, nv_varid, 'cartesian_axis', 'N')
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             istat = NF90_PUT_ATT(ncid_out, nv_varid, 'axis', 'N')
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
         endif

       ! copy variables
         do ivar = 1, numvars
            if (Vars(ivar)%skip) cycle
            ndim = Vars(ivar)%ndim
            if (.not.Vars(ivar)%static) ndim = ndim+1
            ! remap input dimensions to output dimensions
            do idim = 1, ndim
              axes(idim) = dimids_out(Vars(ivar)%axes(idim))
            enddo
            ! define output variable
            if (trim(Vars(ivar)%name) == trim(name_time_bounds) .and. change_bounds_to_climo) then
                istat = NF90_DEF_VAR (ncid_out, 'climatology_bounds', Vars(ivar)%xtype, &
                                      axes(1:ndim), Vars(ivar)%varid_out)
            else
                istat = NF90_DEF_VAR (ncid_out, Vars(ivar)%name, Vars(ivar)%xtype, &
                                      axes(1:ndim), Vars(ivar)%varid_out)
            endif
            if (istat /= NF90_NOERR) call error_handler &
                           ('defining output variable '//trim(Vars(ivar)%name), ncode=istat)
            if (verbose) write (*,*) 'Defining output variable: '//trim(Vars(ivar)%name)

            ! get compression
            if (in_format == NF90_FORMAT_64BIT .or. in_format == NF90_FORMAT_CLASSIC) then
                shuffle_in = 0
                deflate_in = 0
                deflation_in = 0
            else
                istat = NF90_INQ_VAR_DEFLATE (ncid_in, Vars(ivar)%varid_out, shuffle_in, &
                    deflate_in, deflation_in)
                if (istat /= NF90_NOERR) call error_handler &
                    ('getting NetCDF4 compression: '//trim(Vars(ivar)%name), ncode=istat)
            endif
            if (deflation_in > max_input_deflation) max_input_deflation = deflation_in
            if (shuffle_in > max_input_shuffle) max_input_shuffle = shuffle_in

            ! set compression
            call apply_compression_defaults (deflation_in, shuffle_in, deflation_out, shuffle_out )
            if (deflation_out > 0) then
                istat = NF90_DEF_VAR_DEFLATE (ncid_out, Vars(ivar)%varid_out, shuffle_out, 1, &
                    deflation_out)
                if (istat /= NF90_NOERR) call error_handler &
                           ('setting NetCDF4 compression: '//trim(Vars(ivar)%name), ncode=istat)
            endif

            ! copy variable attributes
            do i=1,NF90_MAX_NAME; cell_methods(i:i) = ' '; enddo
            do i=1,NF90_MAX_NAME; name_bnds   (i:i) = ' '; enddo
            do attnum = 1, Vars(ivar)%natt
               istat = NF90_INQ_ATTNAME (ncid_in, ivar, attnum, attname)
               if (istat /= NF90_NOERR) call error_handler (ncode=istat)
               if (trim(attname) == 'time_avg_info') then
                  !istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, attname, &
                  !     trim(avg_name)//'_T1,'//trim(avg_name)//'_T2,'//trim(avg_name)//'_DT')
                  !if (istat /= NF90_NOERR) call error_handler &
                  !   ('creating attribute time_avg_info for variable '//trim(Vars(ivar)%name), &
                  !        ncode=istat)
               else if (trim(attname) == 'cell_methods' .and. .not.Vars(ivar)%static) then
                   istat = NF90_GET_ATT (ncid_in, ivar, 'cell_methods', cell_methods)
                   if (istat /= NF90_NOERR) call error_handler  &
                             ('getting cell_methods'//trim(Vars(ivar)%name), ncode=istat)
               else if (trim(Vars(ivar)%name) == trim(name_recdim) .and. &
                        trim(attname) == 'bounds' .and. do_climo_avg) then
                   ! change bounds to climatology
                   istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, 'climatology', 'climatology_bounds')
                   if (istat /= NF90_NOERR) call error_handler (ncode=istat)
               else
                   istat = NF90_COPY_ATT (ncid_in, ivar, attname, ncid_out, Vars(ivar)%varid_out)
                   if (istat /= NF90_NOERR) call error_handler &
                        ('copying attribute for variable '//trim(Vars(ivar)%name), ncode=istat)
               endif
            enddo
            !--- when necessary create time average attributes for all variables ---
            if (Vars(ivar)%static) cycle
            if (trim(Vars(ivar)%name) == trim(name_time_bounds)) cycle
            if (trim(Vars(ivar)%name) == trim(name_recdim))      cycle
            ! CF convention
            cell_methods_tavg = trim(name_recdim)//': mean'
            do_period_out = .false.
            if (do_climo_avg) do_period_out = .true.
            if (index(trim(lowercase(cell_methods)),trim(name_recdim)//': point') > 0) do_period_out = .true.
            if (do_period_out) cell_methods_tavg = trim(cell_methods_tavg)//' over '//trim(period_out)
            if (do_time_bounds) then
               if (index(trim(lowercase(cell_methods)),trim(lowercase(cell_methods_tavg))) == 0) then
                   if (len_trim(cell_methods) > 0 .and. do_period_out .and. &
                       index(trim(lowercase(cell_methods)),'within') == 0) then
                       cell_methods = trim(cell_methods)//' within '//trim(period_in)
                   else if (len_trim(cell_methods) == 0 .and. Vars(ivar)%do_avg .and. do_period_out) then
                       cell_methods = trim(name_recdim)//': mean within '//trim(period_in)
                   endif
                   cell_methods = trim(cell_methods)//' '//trim(cell_methods_tavg)
               endif
              !if (trim(lowercase(cell_methods_tavg)) /= trim(lowercase(cell_methods))) then
              !    cell_methods = trim(cell_methods)//' '//trim(cell_methods_tavg)
              !endif
               istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, &
                                     'cell_methods', trim(cell_methods))
               if (istat /= NF90_NOERR) call error_handler ('creating attribute '// &
                   'cell_methods for variable '//trim(Vars(ivar)%name), ncode=istat)
            endif
         enddo

        ! extract time units from time axis units attribute
         nca = len_trim(time_att_units)
         ncu = index(time_att_units(1:nca),' ')-1  ! use first word only

        ! more metadata needed when adding time bounds variable
         if (add_time_bounds) then
             ! add attribute to time axis
             ivar = varid_recdim
             if (do_climo_avg) then
                 istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, 'climatology', &
                                       'climatology_bounds')
             else
                 istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, 'bounds', &
                                       trim(Vars(ivar)%name)//'_bounds')
             endif
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)

             ! add time bounds variable if necessary
             ivar = nvar+1
             istat = NF90_DEF_VAR (ncid_out, Vars(ivar)%name, Vars(ivar)%xtype, &
                                      (/nv_dimid,recdim/), Vars(ivar)%varid_out)
             if (verbose) write (*,*) 'Defining output variable: '//trim(Vars(ivar)%name)
             if (istat /= NF90_NOERR) call error_handler &
                       ('defining output variable '//trim(Vars(ivar)%name), ncode=istat)

            ! set compression
            call apply_compression_defaults (max_input_deflation, max_input_shuffle, &
                deflation_out, shuffle_out )
            if (deflation_out > 0) then
                istat = NF90_DEF_VAR_DEFLATE (ncid_out, Vars(ivar)%varid_out, shuffle_out, 1, &
                    deflation_out)
                if (istat /= NF90_NOERR) call error_handler &
                           ('setting NetCDF4 compression: '//trim(Vars(ivar)%name), ncode=istat)
            endif

             ! plus attributes
             istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, 'long_name', &
                                      trim(name_recdim)//' axis boundaries')
             if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             if (ncu > 0) then
               istat = NF90_PUT_ATT (ncid_out, Vars(ivar)%varid_out, 'units', time_att_units(1:ncu))
               if (istat /= NF90_NOERR) call error_handler (ncode=istat)
             endif
         endif

        ! end of defining metadata for output file
         istat = NF90_ENDDEF (ncid_out,header_buffer_val,4,0,4)
         if (istat /= NF90_NOERR) call error_handler (ncode=istat)

     endif

   !------------ END OF OUTPUT FILE INITIALIZATION --------------
   !-------------------------------------------------------------

    !---- zero out sums -----
     do ivar = 1, nvar
        Vars(ivar)%data_sum = 0.0
        Vars(ivar)%tavg(3)  = 0.0
        if (Vars(ivar)%use_data_wt) then
            Vars(ivar)%data_wt  = 0.0
            Vars(ivar)%data_wt_max = 0.
        else
            Vars(ivar)%data_wt_max = 0.
        endif
     enddo

    !--- this will not really be used ---
    !  initialize to avoid divid by zero
     if (nvar_out == nvar+1) then
        ivar = nvar_out
        Vars(ivar)%data_sum = 0.0
        Vars(ivar)%tavg(3)  = 0.0
        Vars(ivar)%data_wt  = 0.0
        Vars(ivar)%data_wt_max = 1.
     endif

!-----------------------------------------------------------------------
!-------------------- Loop through time axis ---------------------------
!-----------------------------------------------------------------------

     do itime = 0, numtime

      ! get time coordinate value
        if (itime > 0) then
            istat =  NF90_GET_VAR (ncid_in, varid_recdim, time, start=(/itime/))
            if (istat /= NF90_NOERR) call error_handler ('getting time coord value', ncode=istat)
        endif

      !--- read time bnds_info
        istat = NF90_GET_VAR (ncid_in, time_bnds_id, tavg(1:2), (/itime/))
        if (istat /= NF90_NOERR) call error_handler &
          ('reading time bnds', ncode=istat)
        tavg(3) = tavg(2) - tavg(1)
        if (do_climo_avg .and. itime == numtime) then
          ! use midpoint of last time interval for climatological time
          if (change_bounds_to_climo) then
            climo_time = 0.5*(tavg(1)+tavg(2))
          else
            ! use last time if already climo time
            climo_time = time
          endif
        endif
!-----------------------------------------------------------------------
!-------------------- Loop through variables ---------------------------
!-----------------------------------------------------------------------

     do ivar = 1, nvar

      !--- skip certain fields ---
        if (.not.Vars(ivar)%static .and. itime == 0) cycle
        if (     Vars(ivar)%static .and. itime  > 0) cycle
        if (     Vars(ivar)%skip                   ) cycle

      !--- get variable data and check dimensions ----
      !  use name from first file to get variable id
         istat = NF90_INQ_VARID (ncid_in, Vars(ivar)%name, varid )
         if (istat /= NF90_NOERR) call error_handler ('variable '//trim(Vars(ivar)%name)// &
                                ' not found in file '//trim(file_names(ifile)), ncode=istat)
         istat = NF90_INQUIRE_VARIABLE (ncid_in, varid, ndims=ndim, dimids=axes)
         if (istat /= NF90_NOERR) call error_handler (ncode=istat)
         ! determine shape
         shape = 1
         do idim = 1, ndim
            if (axes(idim) == recdim) cycle
            istat = NF90_INQUIRE_DIMENSION (ncid_in, axes(idim), len=shape(idim))
            if (istat /= NF90_NOERR) call error_handler (ncode=istat)
         enddo
         ! check shape (must match shape from first file)
         if (count(shape /= Vars(ivar)%shape) > 0) call error_handler &
                  ('inconsistent shape between files for variable '//trim(Vars(ivar)%name))

       ! setup variable size and start indices
         start = 1
         if (.not.Vars(ivar)%static) start(Vars(ivar)%ndim+1) = itime
         ix = Vars(ivar)%shape(1);  jy = Vars(ivar)%shape(2);  kz = Vars(ivar)%shape(3)
         mt = Vars(ivar)%shape(4)

       !--- read variable ---
         select case (Vars(ivar)%ndim)
            case(0)
               istat = NF90_GET_VAR (ncid_in, varid, ddata(1,1,1,1), start(1:1))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(1)
               istat = NF90_GET_VAR (ncid_in, varid, ddata(1:ix,1,1,1), start(1:2))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(2)
               istat = NF90_GET_VAR (ncid_in, varid, ddata(1:ix,1:jy,1,1), start(1:3))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(3)
               istat = NF90_GET_VAR (ncid_in, varid, ddata(1:ix,1:jy,1:kz,1), start(1:4))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(4)
               istat = NF90_GET_VAR (ncid_in, varid, ddata(1:ix,1:jy,1:kz,1:mt), start(1:5))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         end select

       !--- keep track of time averaging info: begin, end, and delta times ---
         if (itime > 0) then
            if (Vars(ivar)%do_avg) then
               if (itime == 1)       Vars(ivar)%tavg(1) = tavg(1)   ! begin time
               if (itime == numtime) Vars(ivar)%tavg(2) = tavg(2)   ! end time
               Vars(ivar)%tavg(3)  = Vars(ivar)%tavg(3) + tavg(3)   ! total time avg length
            else
              !----- no time average info -----
              !    assume equal avg intervals
              !    this will only work if numtime > 1
               if (itime == 1)      Vars(ivar)%tavg(1) = time
               if (itime == 2) then
                   Vars(ivar)%tavg(3) = time-Vars(ivar)%tavg(1)
                   Vars(ivar)%tavg(1) = Vars(ivar)%tavg(1)-Vars(ivar)%tavg(3)
                   Vars(ivar)%tavg(3) = Vars(ivar)%tavg(3)*real(numtime)
               endif
               if (itime == numtime) Vars(ivar)%tavg(2) = time
               if (numtime == 1    ) Vars(ivar)%tavg(3) = 1.0 ! avoid divid by zero
            endif
         endif

        !---------- SUM FIELDS ----------

         if (Vars(ivar)%use_data_wt) then
           if (Vars(ivar)%do_missval) then
              where (ddata(1:ix,1:jy,1:kz,1:mt) .ne. Vars(ivar)%missval)
                Vars(ivar)%data_sum = Vars(ivar)%data_sum+tavg(3)*ddata(1:ix,1:jy,1:kz,1:mt)
                Vars(ivar)%data_wt  = Vars(ivar)%data_wt +tavg(3)
              endwhere
           else
                Vars(ivar)%data_sum = Vars(ivar)%data_sum+tavg(3)*ddata(1:ix,1:jy,1:kz,1:mt)
                Vars(ivar)%data_wt  = Vars(ivar)%data_wt +tavg(3)
           endif
           Vars(ivar)%data_wt_max = Vars(ivar)%data_wt_max+tavg(3)
         else
           if (numtime == 1) then
                Vars(ivar)%data_sum = ddata(1:ix,1:jy,1:kz,1:mt)
           else
                Vars(ivar)%data_sum = Vars(ivar)%data_sum+tavg(3)*ddata(1:ix,1:jy,1:kz,1:mt)
                Vars(ivar)%data_wt_max = Vars(ivar)%data_wt_max+tavg(3)
           endif
         endif

!-----------------------------------------------------------------------
     enddo    ! variable loop
!-----------------------------------------------------------------------
     enddo    ! time loop
!-----------------------------------------------------------------------

!--- check length of average ---
!--- must be the same for all fields ---

     ! find the first avg field
     first = .true.
     do ivar = 1, nvar
        if (Vars(ivar)%do_avg) then
           time_bounds(1:3) = Vars(ivar)%tavg(1:3)
           first = .false.
           exit
        endif
     enddo

     ! check time bounds for all variable
     do ivar = 1, nvar
       if (Vars(ivar)%static) cycle
       if (Vars(ivar)%skip  ) cycle
       if (trim(Vars(ivar)%name) == trim(name_time_bounds)) cycle
       if (trim(Vars(ivar)%name) == trim(name_recdim))      cycle
       if (first) then
          ! in case all vars are do_avg=F
          time_bounds(1:3) = Vars(ivar)%tavg(1:3)
          first = .false.
       else
          if (count(time_bounds /= Vars(ivar)%tavg) > 0) then
              if (.not.Vars(ivar)%do_avg .and. skip_tavg_errors) then
                  if (.not.suppress_warnings) write (*,98) trim(Vars(ivar)%name)
               98 format ('WARNING: average information does not agree for variable ',a)
              else
                  print *, 'time_bounds for file=',time_bounds
                  print *, 'tavg for var=',Vars(ivar)%tavg
                  call error_handler &
                    ('average information does not agree for variable '//trim(Vars(ivar)%name))
              endif
          endif
       endif
     enddo

     ! save last time value from the input file
     last_time = time

!=======================================================================
!------------------------ OUTPUT SECTION -------------------------------
!=======================================================================

  !---- write variables with no time axis ----
   if (ifile == 1) then
      do ivar = 1, nvar_out
         if (.not.Vars(ivar)%static) cycle
         if (     Vars(ivar)%skip)   cycle

         select case (Vars(ivar)%ndim)
            case(0)
               istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, Vars(ivar)%data_sum(1,1,1,1))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(1)
               istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, Vars(ivar)%data_sum(:,1,1,1))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(2)
               istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, Vars(ivar)%data_sum(:,:,1,1))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(3)
               istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, Vars(ivar)%data_sum(:,:,:,1))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
            case(4)
               istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, Vars(ivar)%data_sum(:,:,:,:))
               if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         end select
      enddo
   endif

  !---- output time coordinate value ----

    if (use_end_time) then
        time_out = last_time
    else
        time_out = 0.5*(time_bounds(1)+time_bounds(2))
        if (do_climo_avg) time_out = climo_time
    endif

    if (verbose) then
      print '(a,4f15.3)', 'time_out,time_bounds=',time_out,time_bounds
    endif

  !----- OUTPUT TIME AVERAGES FOR STATIC FIELDS ------

   do ivar = 1, nvar_out
      if (Vars(ivar)%static) cycle    ! process only non-static fields
      if (Vars(ivar)%skip)   cycle

    ! size of this variable
      ix = Vars(ivar)%shape(1);  jy = Vars(ivar)%shape(2);  kz = Vars(ivar)%shape(3)
      mt = Vars(ivar)%shape(4)

    ! compute time average
    ! divid by data weight when appropriate
    ! insert missing value where needed
      if (Vars(ivar)%use_data_wt) then
          data_wt_min = max(min(frac_valid_data,1.-small),0.)*Vars(ivar)%data_wt_max
          if (Vars(ivar)%do_missval) then
             where (Vars(ivar)%data_wt > data_wt_min)
                ddata(1:ix,1:jy,1:kz,1:mt) = Vars(ivar)%data_sum/Vars(ivar)%data_wt
             elsewhere
                ddata(1:ix,1:jy,1:kz,1:mt) = Vars(ivar)%missval
             endwhere
          else
             ! this where loop should always be true
             where (Vars(ivar)%data_wt > data_wt_min)
                ddata(1:ix,1:jy,1:kz,1:mt) = Vars(ivar)%data_sum/Vars(ivar)%data_wt
             endwhere
          endif
      else
          if (numtime == 1) then
             ddata(1:ix,1:jy,1:kz,1:mt) = Vars(ivar)%data_sum
          else
             ddata(1:ix,1:jy,1:kz,1:mt) = Vars(ivar)%data_sum/Vars(ivar)%data_wt_max
          endif
      endif

    ! correct non-float data for truncation error
      if (Vars(ivar)%xtype /= NF90_REAL8 .and. Vars(ivar)%xtype /= NF90_REAL4) then
          if (Vars(ivar)%do_missval) then
              where (ddata(1:ix,1:jy,1:kz,1:mt) /= Vars(ivar)%missval) &
                ddata(1:ix,1:jy,1:kz,1:mt) = ddata(1:ix,1:jy,1:kz,1:mt) + sign(half,ddata(1:ix,1:jy,1:kz,1:mt))
          else
                ddata(1:ix,1:jy,1:kz,1:mt) = ddata(1:ix,1:jy,1:kz,1:mt) + sign(half,ddata(1:ix,1:jy,1:kz,1:mt))
          endif
      endif

    ! starting index in output file for this time
      start = 1
      if (.not.Vars(ivar)%static) start(Vars(ivar)%ndim+1) = ifile

    !**** special processing for variable "time_bounds" ****
      if (trim(Vars(ivar)%name) == trim(name_time_bounds)) then
          istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, time_bounds(1:2), start(1:2))
          if (istat /= NF90_NOERR) call error_handler ('var= '//trim(Vars(ivar)%name), ncode=istat)
          cycle
      endif

    !**** special processing for record dimension ****
      if (trim(Vars(ivar)%name) == trim(name_recdim)) then
          istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, time_out, start(1:1))
          if (istat /= NF90_NOERR) call error_handler ('writing time coord', ncode=istat)
          cycle
      endif

    !---- write data to output file ----
      select case (Vars(ivar)%ndim)
         case(0)
            istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, ddata(1,1,1,1), start(1:1))
            if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         case(1)
            istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, ddata(1:ix,1,1,1), start(1:2))
            if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         case(2)
            istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, ddata(1:ix,1:jy,1,1), start(1:3))
            if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         case(3)
            istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, ddata(1:ix,1:jy,1:kz,1), start(1:4))
            if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
         case(4)
            istat = NF90_PUT_VAR (ncid_out, Vars(ivar)%varid_out, ddata(1:ix,1:jy,1:kz,1:mt), start(1:5))
            if (istat /= NF90_NOERR) call error_handler &
                                        ('var= '//trim(Vars(ivar)%name), ncode=istat)
      end select

   enddo ! end of output variable loop

   ! write time average info (probably do only once)
     do i = 1, 3
        istat = NF90_PUT_VAR (ncid_out, tid(i), time_bounds(i), (/ifile/))
        if (istat /= NF90_NOERR) call error_handler ('writing time avg info', ncode=istat)
     enddo

     if (add_time_bounds) then
        istat = NF90_PUT_VAR(ncid_out, nv_varid, (/0,1/))
        if (istat /= NF90_NOERR) call error_handler ('writing nv bounds info', ncode=istat)
     endif


   ! close input file
     istat = NF90_CLOSE (ncid_in)
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)

 enddo
!---------------- end of input file (ifile) loop --------------------

! close output file
  istat = NF90_CLOSE (ncid_out)
  if (istat /= NF90_NOERR) call error_handler (ncode=istat)

contains
!#######################################################################

 subroutine error_handler (string, ncode)
 character(len=*), intent(in), optional :: string
 integer         , intent(in), optional :: ncode
 character(len=80) :: errstrg
   if (present(ncode))  then
       print *, 'NETCDF ERROR in program time_average'
   else
       print *, 'ERROR in program time_average'
   endif
   if (present(string)) print *, trim(string)
   if (present(ncode))  then
       errstrg = NF90_STRERROR (ncode)
       print *, trim(errstrg)
   endif
   call abort ()
!  stop 111
!  call error_mesg ('program time_average', string, FATAL)
 end subroutine error_handler

!#######################################################################

! change string to all lower case

 function lowercase (cs)
 character(len=*), intent(in) :: cs
 character(len=len(cs))       :: lowercase
 character :: ca(len(cs))

 integer, parameter :: co=iachar('a')-iachar('A') ! case offset

    ca = transfer(cs,"x",len(cs))
    where (ca >= "A" .and. ca <= "Z") ca = achar(iachar(ca)+co)
    lowercase = transfer(ca,cs)

 end function lowercase

!#######################################################################

 subroutine init_variable_type ( ncid, varid, Var )
 integer            , intent(in)    :: ncid, varid
 type(variable_type), intent(inout) :: Var

 integer :: idim, nc, istat, i
 logical :: do_fillval
 real    :: fillval

    ! will not allow fields with more than five dimensions
    if (Var%ndim > 5) call error_handler &
                ('input field ('//trim(Var%name)//') has more than 5 axes')

    ! skip all fields stored as 1-byte data
    if (Var%xtype == NF90_CHAR .or. Var%xtype == NF90_INT1) then
        Var%skip = .true.
        write (*,100) trim(Var%name), '(data type not supported)'
    100 format (4x,'Skipping field=',a,4x,a)
        return
    endif

   !--- variable has a time axis (non-static) ---
    if (Var%axes(Var%ndim) == recdim) then ! assume last dimension is time
        Var%static = .false.
        Var%ndim = Var%ndim - 1
        ! skip certain fields that only have record dimension
        if (Var%ndim == 0) then
            nc = len_trim(Var%name)
            if (Var%name(nc-2:nc) == '_T1' .or. &
                Var%name(nc-2:nc) == '_T2' .or. &
                Var%name(nc-2:nc) == '_DT' .or. &
                Var%name(nc-6:nc) == '_NITEMS') Var%skip = .true.
            ! skip record dimension field
         !!!if (Var%axes(1) == recdim) Var%skip = .true.
        endif
        if (Var%skip) then
            if (verbose) write (*,100) trim(Var%name), ' '
            return
        endif
   !--- variable does not have a time axis (static) ---
    else
        Var%static = .true.
    endif

    ! check compilation precision for double precision variables
    ! skip for time_bounds and record dimension
    if (Var%xtype == NF90_REAL8) then
        if (trim(Var%name).ne.trim(name_time_bounds) .and. &
            trim(Var%name).ne.trim(name_recdim)) then
                if (precision(ddata) < precision(time)) then
                    if (.not.suppress_warnings) write (*,90) trim(Var%name)
                90  format ('WARNING: double precision field ',a, &
                            ' will be averaged using lower precision')
                    endif
        endif
    endif

   !--- compute variable size/shape ---
    shape(1:4) = 1
    do idim = 1, Var%ndim
       istat = NF90_INQUIRE_DIMENSION (ncid, Var%axes(idim), len=shape(idim))
       if (istat /= NF90_NOERR) call error_handler (ncode=istat)
    enddo
    Var%shape = shape

   !--- check shape of "time_bounds" ---
    if (trim(Var%name) == trim(name_time_bounds)) then
        if (Var%ndim > 1 .or. Var%shape(1) /= 2) then
            if (.not.suppress_warnings) write (*,95)
         95 format ('WARNING: wrong dimensions for time_bounds ... skipping')
            Var%skip = .true.
            return
        endif
    endif

   !--- change the name of the time bounds for climatology bounds ----
   !if (trim(Var%name) == trim(name_time_bounds)) then
   !    if (change_bounds_to_climo) then
   !        do i=1,NF90_MAX_NAME; name_time_bounds(i:i) = ' '; enddo
   !        do i=1,NF90_MAX_NAME; Var%name        (i:i) = ' '; enddo
   !        name_time_bounds = 'climatology_bounds'
   !        Var%name         = 'climatology_bounds'
   !    endif
   !endif

   !--- save the maximum size ---
    do idim = 1, Var%ndim
       max_shape(idim) = max( max_shape(idim), shape(idim) )
    enddo

   !--- get missing value if present ---
    Var%do_missval = .false.
    istat = NF90_GET_ATT (ncid, varid, 'missing_value', Var%missval)
    if (istat == NF90_NOERR) Var%do_missval = .true.
   !--- also check for fill value ---
    do_fillval = .false.
    istat = NF90_GET_ATT (ncid, varid, '_FillValue', fillval)
    if (istat == NF90_NOERR) do_fillval = .true.
   !--- use fill value when missing value not present ---
    if (.not.Var%do_missval .and. do_fillval) then
        Var%do_missval = do_fillval
        Var%missval    = fillval
    endif

   !--- get time average info for non-static variables ----
    Var%do_avg = .false.
    if (.not.Var%static) then
      Var%do_avg = .true.
    endif

 end subroutine init_variable_type

!#######################################################################

 subroutine check_for_climo_avg ( ncid, tname, tunits, ntimes,  &
                                  climo, period_in, period_out, &
                                  tbnds_name)
 integer,          intent(in)  :: ncid
 character(len=*), intent(in)  :: tname, tunits
 integer,          intent(in)  :: ntimes
 logical,          intent(out) :: climo
 character(len=*), intent(out) :: period_in, period_out
 character(len=*), intent(in)  :: tbnds_name
 real(8) :: tbnds(2,2), times(2), dt(2)
 integer :: varid(2), i, istat, nc
 logical :: tbnds_present

 ! read the time bounds for the two time periods

  if (tbnds_name(1:1) .ne. ' ') tbnds_present = .true.
  if (tbnds_present) then
      istat = NF90_INQ_VARID ( ncid, trim(tbnds_name), varid(1) )
      if (istat /= NF90_NOERR) call error_handler ('error getting varid for time bounds', ncode=istat)
      istat = NF90_GET_VAR ( ncid, varid(1), tbnds )
      if (istat /= NF90_NOERR) call error_handler ('error getting data for time bounds', ncode=istat)
  else
      tbnds = 0.0
  endif

! if the time bounds are contiguous this is not a climatological average

  climo = .false.
  if (tbnds_present) then
     if (tbnds(2,1) == tbnds(1,2)) then
        climo = .false.
     else
        climo = .true.
     endif
  endif

! get the first two time axis values

  if (ntimes >= 2) then
     istat = NF90_INQ_VARID ( ncid, trim(tname), varid(1) )
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     istat = NF90_GET_VAR ( ncid, varid(1), times )
     if (istat /= NF90_NOERR) call error_handler (ncode=istat)
     dt = times(2)-times(1)
  else
     if (.not.tbnds_present) then
         call error_handler ('cannot average when only one time level')
     endif
  endif

! convert time units to days

  if (trim(tunits) == 'days') then
  else if (trim(tunits) == 'hours') then
       tbnds = tbnds/24.
  else if (trim(tunits) == 'minutes') then
       tbnds = tbnds/1440.
  else if (trim(tunits) == 'seconds') then
       tbnds = tbnds/86400.
  else
     !call error_handler ('invalid units')
      tbnds = 0.
  endif

! determine the "within average period"
! this is the length of the input period

  if (tbnds_present) then
      dt(1) = tbnds(2,1)-tbnds(1,1)
  else
      dt(1) = times(2)-times(1)
      dt(2) = ntimes*dt(1) ! estimate total length of average
  endif

! choose days, months, years
  if (dt(1) > 27. .and. dt(1) < 300.) then
      period_in = 'months'
  else if (dt(1) < 1.01) then
      period_in = 'days'
  else if (dt(1) > 359.) then
      period_in = 'years'
  else
      period_in = ' ' ! unknown
  endif

! determine the "over average period"
! this is the time between input values when climatology
! this is the averaging time when not climatology

  if (tbnds_present) then
      if (climo) then
          dt(2) = tbnds(1,2)-tbnds(1,1)
      else
          dt(2) = ntimes*(times(2)-times(1)) ! estimate total length of average
      endif
  endif

! choose days, months, years
  if (dt(2) > 27. .and. dt(2) < 300.) then
      period_out = 'months'
  else if (dt(1) < 1.01) then
      period_out = 'days'
  else if (dt(2) > 359.) then
      period_out = 'years'
  else
      period_out = ' ' ! unknown
  endif


 end subroutine check_for_climo_avg

!#######################################################################

 subroutine alloc_variable_type ( ntime, Var )
 integer            , intent(in)    :: ntime
 type(variable_type), intent(inout) :: Var

   !--- allocate space for data (and data wt when necessary) ---
    allocate (Var%data_sum(shape(1),shape(2),shape(3),shape(4)))
    if (.not.Var%use_data_wt) then
       if (ntime > 1 .and. Var%do_missval .and. .not.Var%static) then
           allocate (Var%data_wt(shape(1),shape(2),shape(3),shape(4)))
           Var%use_data_wt = .true.
       endif
    endif

    if (verbose) then  ! global variable
        write (*,102) trim(Var%name), Var%static, &
                      Var%use_data_wt, Var%ndim, Var%shape
        if (Var%do_missval) write (*,103) Var%missval
        if (Var%do_avg)   write (*,104) Var%do_avg
    102 format (4x,'Field=',a,4x,'static=',L1,4x,'use_data_wt=',L1, &
              /,10x,'ndim=',i1,4x,'shape=',3i4)
    103 format (10x,'missing value=',g15.7)
    104 format (10x,'time averaged=',L1)
    endif

 end subroutine alloc_variable_type

!#######################################################################
! apply_compression_defaults -- applies the user-specified settings
!   deflation and shuffle to the input variable parameters.
!   If the user-specified settings are -1, use the input file settings.

subroutine apply_compression_defaults (deflation_in, shuffle_in, deflation_out, shuffle_out )
    integer, intent(in) :: deflation_in, shuffle_in
    integer, intent(out) :: deflation_out, shuffle_out

    if (user_deflation == -1) then
        deflation_out = deflation_in
    else if (user_deflation == 0) then
        deflation_out = 0
    else
        deflation_out = user_deflation
    endif

    if (user_shuffle == -1) then
        shuffle_out = shuffle_in
    else
        shuffle_out = user_shuffle
    endif

    if (verbose) then
        if (deflation_out > 0) then
            print *, 'Using deflation=', deflation_out, 'and shuffle=', shuffle_out
        else
            print *, 'Not using deflation'
        endif
    endif

end subroutine apply_compression_defaults
!#######################################################################

end program time_average
