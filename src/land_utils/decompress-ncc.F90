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
! This program reads several input netcdf files, presumably containing "compressed
! by gathering" data, and combines them into a single output file
!
!-----------------------------------------------------------------------

#include <config.h>

#define CHECK_NF_ERRSTAT(ierr) call nfu_check_err(ierr,__FILE__,__LINE__)
program decompress

  use nfu_mod
  use nfu_compress_mod
  implicit none
  include 'netcdf.inc'
  include 'version.inc'

  integer, parameter :: PATH_MAX = 1024 ! max len of the file name;
  integer, parameter :: HEADERPAD = 16384 ! Use mpp_io large headers;
  integer            :: blksz = 65536  ! blksz must be writable for nf__create

  character(PATH_MAX), allocatable :: files(:) ! names of all files on the command line
  character(PATH_MAX)              :: outfile  ! name of the output file
  integer :: nfiles    ! number of files on command line
  integer :: debug = 0 ! debug output verbosity level
  integer, allocatable :: input(:)             ! netcdf IDs of input files
  integer :: i,iret,ncid,dimid,varid,varid1,xtype,dimlen,vsize,ndims,nvars,natts
  integer :: in_format, cmode
  integer :: dimids(NF_MAX_VAR_DIMS), dimlens(NF_MAX_VAR_DIMS)
  character(NF_MAX_NAME) :: dimname,varname,attname
  logical :: is_dim, is_compressed, has_records, add_missing_value = .FALSE.
  real   , allocatable :: buffer(:)
  logical, allocatable :: mask(:)
  real :: missing
  real :: ocean_value !< the ocean fillvalue
  logical :: do_oceanValue !< A flag to turn on the ocean fillvalue

  ! get command line options and list of files
  call parse_command_line() ! modigies global data!

  call assert(nfiles>0,'at least one input file must be specified')
  if(debug>0) then
     do i = 1,nfiles
        write(*,'("input file",i3,":",a)')i, '"'//trim(files(i))//'"'
     enddo
     write(*,'("output file:",a)')'"'//trim(outfile)//'"'
  endif

  ! open all input files
  allocate(input(nfiles))
  do i = 1,nfiles
     CHECK_NF_ERRSTAT(nf_open(files(i),NF_NOWRITE,input(i)))
     CHECK_NF_ERRSTAT(nf_inq_format(input(i),in_format))
  enddo

  if (in_format==NF_FORMAT_NETCDF4) then
     cmode = NF_NETCDF4
  elseif (in_format==NF_FORMAT_NETCDF4_CLASSIC) then
     cmode=IOR(NF_NETCDF4,NF_CLASSIC_MODEL)
  elseif (in_format==NF_FORMAT_64BIT) then
     cmode=IOR(NF_CLOBBER,NF_64BIT_OFFSET)
     if(debug>0)write(*,'("output file is 64-bit netcdf")')
  elseif (in_format==NF_FORMAT_CLASSIC) then
     cmode=IOR(NF_CLOBBER,NF_CLASSIC_MODEL)
     if(debug>0)write(*,'("output file is 32-bit netcdf")')
  else
     call assert(.false.,'Unknown netCDF format')
  endif

  ! create output file
  CHECK_NF_ERRSTAT(nf__create(outfile,cmode,0,blksz,ncid))

  ! Create netcdf structure in the output NetCDF file, using last input file
  ! as a template.

  ! clone all dimensions except compressed ones; compressed are just skipped
  CHECK_NF_ERRSTAT(nf_inq_ndims(input(nfiles),ndims))
  do dimid = 1,ndims
     CHECK_NF_ERRSTAT(nfu_inq_dim(input(nfiles),dimid,dimname=dimname,dimlen=dimlen,is_unlim=has_records))
     if(nfu_inq_att(input(nfiles),dimname,'compress')==NF_NOERR) cycle
     if(has_records)&
          dimlen=NF_UNLIMITED
     if(debug>0)&
          write(*,*)'defining dimension "'//trim(dimname)//'" with length',dimlen
     CHECK_NF_ERRSTAT(nf_def_dim(ncid,dimname,dimlen,i)) ! i is just a space for id
  enddo

  ! clone all variable definitions, replacing compressed dimensions with sets
  ! of uncompressed ones
  CHECK_NF_ERRSTAT(nf_inq_nvars(input(nfiles),nvars))
  do varid = 1,nvars
     iret = nfu_inq_compressed_var(input(nfiles),varid,varname,xtype,ndims,dimids,&
          natts=natts,is_dim=is_dim, is_compressed=is_compressed)
     if(debug>0)&
          write(*,*)'defining variable "'//trim(varname)//'"'
     CHECK_NF_ERRSTAT(iret) ! just because the line is going to be too long
     if(is_dim.and.is_compressed) cycle
     do i = 1,ndims
        CHECK_NF_ERRSTAT(nf_inq_dimname(input(nfiles),dimids(i),dimname))
        CHECK_NF_ERRSTAT(nf_inq_dimid(ncid,dimname,dimids(i)))
     enddo
     CHECK_NF_ERRSTAT(nf_def_var(ncid,varname,xtype,ndims,dimids,varid1))
     do i=1,natts
        CHECK_NF_ERRSTAT(nf_inq_attname(input(nfiles),varid,i,attname))
        CHECK_NF_ERRSTAT(nf_copy_att(input(nfiles),varid,attname,ncid,varid1))
     enddo
     if(add_missing_value.and.is_compressed) then
        ! check if missing_value or _FillValue attributes are present
        if ( nf_inq_atttype(input(nfiles),varid,'missing_value',iret)/=NF_NOERR .and. &
             nf_inq_atttype(input(nfiles),varid,'_FillValue',iret)/=NF_NOERR ) then
           ! if not, define the missing value attribute
           select case(xtype)
           case(NF_DOUBLE)
              missing = NF_FILL_DOUBLE
           case(NF_FLOAT)
              missing = NF_FILL_FLOAT
           case(NF_INT)
              missing = NF_FILL_INT
           end select
           ! and add it to the output variable
           CHECK_NF_ERRSTAT(nf_put_att_double(ncid,varid1,'missing_value',xtype,1,missing))
        endif
     endif
  enddo

  ! clone all global attributes
  CHECK_NF_ERRSTAT(nf_inq_natts(input(nfiles),natts))
  do i = 1,natts
     CHECK_NF_ERRSTAT(nf_inq_attname(input(nfiles),NF_GLOBAL,i,attname))
     CHECK_NF_ERRSTAT(nf_copy_att(input(nfiles),NF_GLOBAL,attname,ncid,NF_GLOBAL))
  enddo

  ! ---- end of definition stage
  CHECK_NF_ERRSTAT(nf__enddef(ncid,HEADERPAD,4,0,4))


  ! grow unlimited dimension, if necessary
  do varid=1,nvars
     CHECK_NF_ERRSTAT(nfu_inq_compressed_var(input(nfiles),varid,name=varname,dimlens=dimlens,has_records=has_records))
     if(has_records) then
        ! just write an integer at the very end of the variable -- that extends the
        ! record dimensions as well
        CHECK_NF_ERRSTAT(nf_inq_varid(ncid,varname,varid))
        CHECK_NF_ERRSTAT(nf_put_var1_int(ncid,varid,dimlens,0))
        exit ! this loop
     endif
  enddo

  ! gather and copy data
  CHECK_NF_ERRSTAT(nf_inq_nvars(ncid,nvars))
  do varid = 1,nvars
!     CHECK_NF_ERRSTAT(nfu_inq_compressed_var(ncid,varid,name=varname,varsize=vsize))
     CHECK_NF_ERRSTAT(nfu_inq_var(ncid,varid,varname,xtype,varsize=vsize))
     if(debug>0) &
          write(*,*)'processing var "'//trim(varname)//'"'
     allocate(buffer(vsize),mask(vsize))
     mask(:) = .false.
     do_oceanValue=.false.
     ocean_value = 0.0
     if(nfu_get_att(ncid,varname,'ocean_fillvalue',ocean_value)==NF_NOERR) then
       do_oceanValue=.true.
     endif

     ! obtain the missing value
     if(nfu_get_att(ncid,varname,'missing_value',missing)==NF_NOERR) then
        ! do nothing, the value is already in the "missing" variable
     else if(nfu_get_att(ncid,varname,'_FillValue',missing)==NF_NOERR) then
        ! do nothing, the value is already in the "missing" variable
     else
        ! get fill value for the type instead of the missing value
        select case(xtype)
        case(NF_DOUBLE)
           missing = NF_FILL_DOUBLE
        case(NF_FLOAT)
           missing = NF_FILL_FLOAT
        case(NF_INT)
           missing = NF_FILL_INT
        end select
     endif
     ! fill the buffer with the missing value
     buffer=missing

     ! read the variable
     do i=1,nfiles
        CHECK_NF_ERRSTAT(nfu_get_compressed_var_r8n(input(i),varname,buffer,mask,ocean=do_oceanValue,ocean_value=ocean_value))
     enddo
     ! write the variable
     CHECK_NF_ERRSTAT(nfu_put_var_r8(ncid,varname,buffer))
     deallocate(buffer,mask)
  enddo

  CHECK_NF_ERRSTAT(nf_sync(ncid))
  CHECK_NF_ERRSTAT(nf_close(ncid))

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ---- parses command line arguments, getting options and gathering list of
! file names
! NOTE: updates global variables.
subroutine parse_command_line()
  character(PATH_MAX) :: arg, param, command_name

  integer :: nargs     ! number of command-line arguments
  logical :: do_interpret_arguments
  integer :: i, iostat

  call getarg(0,command_name)

  nargs = command_argument_count()

  if(nargs==0) then
     call usage()
     call exit(EXIT_FAILURE)
  endif

  allocate(files(nargs))  ! allocate storage for all file names
  do_interpret_arguments = .true.
  add_missing_value = .false.
  i=1        ! counter of all command-line arguments
  nfiles = 0 ! counter of input files
  do while (i<=nargs)
     call getarg(i,arg)
     if(debug>1) write(*,*)'argument ',i, trim(arg)
     if(arg(1:1)=='-'.and.do_interpret_arguments) then
        select case(trim(arg))
        case('--')
           do_interpret_arguments = .false.

        case('-D','--debug-level','-v','--verbosity-level')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           call getarg(i+1,param)
           read(param,*,iostat=iostat) debug
           call assert(iostat==0,trim(arg)//' flag must be followed by integer verbosity level')
           i=i+1

        case('-m','--add-missing-value')
           add_missing_value = .TRUE.

        case ('-h','-?','--help')
           call usage()
           call exit(EXIT_SUCCESS)

        case ('-V', '--version')
           call print_version(trim(command_name))
           call exit(EXIT_SUCCESS)

        case default
           call usage()
           call assert(.false.,'argument "'//trim(arg)//'" is illegal')
        end select
     else
        ! argument is either input or output file name, add it to the list
        nfiles = nfiles+1
        files(nfiles) = arg
     endif
     i = i+1
  enddo
  if (nfiles>0) then
     outfile = files(nfiles)
     nfiles  = nfiles-1
  endif
end subroutine


! ---- prints usage information
subroutine usage()
! this program reads several input netcdf files, presumably containing "compressed
! by gathering" data, and combines them into a single output file
  character(len=PATH_MAX) :: name
  call getarg(0,name)
  write(*,'(a)')'Usage:'
  write(*,'(a)')'  '//trim(name)//' [-v verbosity-level] [-m] in.nc [...] out.nc'
  write(*,'(a)')'Converts one or several compressed-by-gathering netcdf file into'
  write(*,'(a)')'one regular netcdf. Normally used to convert lm3 restarts into a'
  write(*,'(a)')'form suitable for visualization applications.'
  write(*,'(a)')
  write(*,'(a)')'  -v verbosity-level   Specifies level of debug output verbosity'
  write(*,'(a)')'  -m                   Forces adding a missing_value attribute to the variables'
  write(*,'(a)')'                       that do not have it'
  write(*,'(a)')'  -h, --help           display this help and exit'
  write(*,'(a)')'  -V, --version        output version information and exit'
  write(*,'(a)')'  in.nc                Input file name(s)'
  write(*,'(a)')'  out.nc               Output file name'
  write(*,'(a)')
  write(*,'(a)')'WARNING: output file is overwritten.'
end subroutine

! ===========================================================================
! ---- prints error message an exits if condition is not satisfied
subroutine assert(cond,message)
  logical     , intent(in) :: cond    ! condition to check
  character(*), intent(in) :: message ! error message to print if condition is not satisfied

  if(.not.cond) then
     write(*,*) 'ERROR :: ',trim(message)
     call exit(EXIT_FAILURE)
  endif
end subroutine

end program decompress
