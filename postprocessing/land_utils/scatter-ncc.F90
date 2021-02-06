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
! This program is distributed under the terms of the GNU General Public
! License.
!
! This program reads several input netcdf files, presumably containing
! "compressed by gathering" data, and combines them into a single output file
!
!-----------------------------------------------------------------------

#define __NF_ASRT__(ierr) call nfu_check_err(ierr,__FILE__,__LINE__)
program combine_res

  use nfu_mod
  use nfu_compress_mod
  implicit none
  include 'netcdf.inc'

  integer, parameter :: PATH_MAX = 1024 ! max len of the file name;
  integer, parameter :: HEADERPAD = 16384 ! Use mpp_io large headers;
  integer            ::  blksz = 65536  ! blksz must be writable for nf__create

  character(PATH_MAX), allocatable :: files(:) ! names of all files on the command line
  character(PATH_MAX)              :: outfile  ! name of the output file
  character(PATH_MAX)              :: infile   ! name of input file
  integer :: nfiles    ! number of files on command line
  integer :: nfiles_out ! number of output file
  integer :: npex_io    ! number of division in x-direction
  integer :: npey_io    ! number of division in y-direction
  integer :: nlon, nlat             ! global domain size
  integer :: nlon_local, nlat_local ! io_domain size
  integer :: vsize_local, recsize
  integer :: debug = 0 ! debug output verbosity level
  integer :: in_ncid
  integer, allocatable :: out_ncid(:), dimlen_list(:)
  integer :: dimids(NF_MAX_VAR_DIMS)
  integer :: i,dimid,varid,dimlen,vsize,ndims,nvars,ngatts
  integer :: n, l, ll, j, ii, jj, nn, t, npts, npts_local, ij
  integer :: dimlens(NF_MAX_DIMS)
  logical :: has_records, is_compressed
  integer :: in_format ! format of input files
  integer :: cmode     ! mode for output file creation
  character(NF_MAX_NAME) :: dimname,varname,attname
  real   , allocatable :: buffer(:), buffer_2d(:,:)
  logical, allocatable :: mask(:), mask_2d(:,:)
  integer :: nz, start(4), nread(4), nwrite(4), k, nrec, tlev

  integer :: t_id, k_id, nz_saved

  ! get command line options and list of files
  call parse_command_line() ! modigies global data!

  if(debug>0) then
     write(*,'("input file: ",a)') trim(infile)
  endif

  ! open input files and determine the creation mode of output file:
  ! if the input file is 64-bit then the output is 64-bit as well,
  ! otherwise it's 32-bit
  __NF_ASRT__(nf_open(infile,NF_NOWRITE,in_ncid))
  __NF_ASRT__(nf_inq_format(in_ncid,in_format))

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

  !--- create the file
  nfiles_out = npex_io*npey_io

  ! get the grid size
  nlon = 0
  nlat = 0
  nz = 1
  __NF_ASRT__(nf_inq_ndims(in_ncid,ndims))
  do dimid = 1,ndims
     __NF_ASRT__(nfu_inq_dim(in_ncid,dimid,dimname=dimname,dimlen=dimlen))
     if( trim(dimname) == "lon" ) then
        nlon = dimlen
     else if( trim(dimname) == "lat" ) then
        nlat = dimlen
     else if( trim(dimname) == "zfull" ) then
        nz = dimlen
     endif
  enddo

  if(nlon==0)  call assert(.false.,'file '//trim(infile)//' does not have dimension nlon')
  if(nlat==0)  call assert(.false.,'file '//trim(infile)//' does not have dimension nlat')

  ! make sure grid size is divisible by number of division.
  if(mod(nlon, npex_io) .NE. 0)  call assert(.false.,' nlon is not divisible by npex_io (-i)')
  if(mod(nlat, npey_io) .NE. 0)  call assert(.false.,' nlat is not divisible by npey_io (-i)')
  nlon_local = nlon/npex_io
  nlat_local = nlat/npey_io
  nfiles_out = npex_io*npey_io

  ! create output file
  allocate(out_ncid(nfiles_out))
  allocate(dimlen_list(nfiles_out))
  do n = 1, nfiles_out
     ! set the outfile name
     write(outfile, '(a,i4.4)') trim(infile)//".", n-1
     __NF_ASRT__(nf__create(outfile,cmode,0,blksz,out_ncid(n)))
  enddo


  ! clone all dimensions; for compressed dimensions calculate the length
  nrec = 1
  nz_saved = nz
  __NF_ASRT__(nf_inq_ndims(in_ncid,ndims))
  do dimid = 1,ndims
     __NF_ASRT__(nfu_inq_dim(in_ncid,dimid,dimname=dimname,dimlen=dimlen,is_unlim=has_records))
     if(debug>0)&
          write(*,*)'defining dimension "'//trim(dimname)//'" with length',dimlen
     if(has_records)then
        nrec = dimlen
        dimlen = NF_UNLIMITED
     endif
     dimlen_list = 0
     is_compressed = .false.
     if(nfu_inq_att(in_ncid,dimname,'compress')==NF_NOERR) then
        is_compressed = .true.
        __NF_ASRT__(nfu_inq_compressed_var(in_ncid,dimname,varsize=vsize))
        allocate(buffer(vsize),mask(vsize))
        mask(:) = .false.
        __NF_ASRT__(nfu_get_compressed_var_r8n(in_ncid,dimname,buffer,mask))
        !--- figure out dimension length for each file
        do l = 1, vsize
           if(mask(l)) then
              ij = mod((l-1), nlon*nlat)+1
              i = mod((ij-1), nlon)+1
              j = (ij-1)/nlon+1
              ii = (i-1)/nlon_local+1
              jj = (j-1)/nlat_local+1
              nn = (jj-1)*npex_io+ii
              dimlen_list(nn) = dimlen_list(nn)+1
           endif
        enddo
        ! can't have 0-length dimension, since it is (mis-)understood by netcdf as
        ! a record one.
        deallocate(buffer,mask)
     endif
     do n = 1, nfiles_out
        if(is_compressed) dimlen = max(dimlen_list(n),1)
        __NF_ASRT__(nf_def_dim(out_ncid(n),dimname,dimlen,i)) ! i is just a space for id
     enddo
  enddo

  ! clone all variable definitions
  __NF_ASRT__(nf_inq_nvars(in_ncid,nvars))
  do i = 1,nvars
     do n = 1, nfiles_out
        __NF_ASRT__(nfu_clone_var(in_ncid,i,out_ncid(n)))
        ! NOTE: since cloning of variable definition relies on dimension names,
        ! each variable tile and compressed dimensions automaticaly get the right
        ! size, as defined while creating dimensions in the output file
     enddo
  enddo

  ! clone all global attributes
  __NF_ASRT__(nf_inq_natts(in_ncid,ngatts))
  do i = 1,ngatts
     do n = 1, nfiles_out
        __NF_ASRT__(nf_inq_attname(in_ncid,NF_GLOBAL,i,attname))
        __NF_ASRT__(nf_copy_att(in_ncid,NF_GLOBAL,attname,out_ncid(n),NF_GLOBAL))
     enddo
  enddo

  ! ---- end of definition stage
  do n = 1, nfiles_out
     __NF_ASRT__(nf__enddef(out_ncid(n),HEADERPAD,4,0,4))
  enddo


  npts = nlon*nlat
  npts_local = nlon_local*nlat_local
  !--- loop through each record
  do tlev = 1, nrec


     ! write out the data
     do varid = 1,nvars
        !-- make sure number of levels is no greater than 2.
        __NF_ASRT__(nfu_inq_var(in_ncid,varid, ndims=ndims, dimids=dimids, dimlens=dimlens,recsize=recsize,has_records=has_records))

        if(.not. has_records .and. tlev>1) cycle

        is_compressed = .false.
        if(ndims>0) then
           !--- restrict that only the first dimension could be compressed dimension.
           __NF_ASRT__(nfu_inq_dim(in_ncid,dimids(1),dimname=dimname))
           if(nfu_inq_att(in_ncid,dimname,'compress')==NF_NOERR) is_compressed = .true.
        endif

        __NF_ASRT__(nfu_inq_compressed_var(in_ncid,varid,name=varname,varsize=vsize, first_dim_only=.true.))

        if(debug>0) &
             write(*,*)'processing var "'//trim(varname)//'"'
        allocate(buffer(vsize),mask(vsize))
        if(is_compressed) then
           vsize_local = vsize/nfiles_out
           allocate(buffer_2d(vsize_local,nfiles_out), mask_2d(vsize_local,nfiles_out))

           k_id = 0
           t_id = 0
           do dimid = 1, ndims
              __NF_ASRT__(nfu_inq_dim(in_ncid,dimids(dimid),dimname=dimname,dimlen=dimlen,is_unlim=has_records))
              if(has_records) then
                 t_id = dimid
              else if(trim(dimname) == "zfull") then
                 k_id = dimid
              endif
           enddo

           if(k_id==0) then
              nz = 1
           else
              nz = nz_saved
           endif
           recsize = recsize/nz

           start(:) = 1; nread(:) = 1
           nread(1) = recsize
           if(t_id > 0) start(t_id) = tlev
           do k = 1, nz
              if(k_id>0) start(k_id)=k
              buffer(:) = 0
              mask(:) = .false.
              mask_2d(:,:) = .false.
              buffer_2d(:,:) = 0
              __NF_ASRT__(nfu_get_compressed_var_r8n(in_ncid,varname,buffer,mask,start=start,count=nread))
              do l = 1, vsize
                 if(mask(l)) then
                    ij = mod((l-1), npts)+1
                    t = (l-1)/npts + 1
                    i = mod((ij-1), nlon)+1
                    j = (ij-1)/nlon+1
                    ii = (i-1)/nlon_local+1
                    jj = (j-1)/nlat_local+1
                    nn = (jj-1)*npex_io+ii
                    ii = mod((i-1),nlon_local) + 1
                    jj = mod((j-1),nlat_local) + 1
                    ll = (t-1)*npts_local + (jj-1)*nlon_local + ii
                    buffer_2d(ll,nn) = buffer(l)
                    mask_2d(ll,nn) = mask(l)
                 endif
              enddo
              do n = 1, nfiles_out
                 if (count(mask_2d(:,n))>0) then
                    start = 1; nwrite = 1
                    nwrite(1) = count(mask_2d(:,n))
                    if(t_id>0) start(t_id) = tlev
                    if(k_id>0) start(k_id) = k
                    __NF_ASRT__(nfu_put_vara_r8(out_ncid(n),varname,start, nwrite, pack(buffer_2d(:,n),mask_2d(:,n))))
                 endif
              enddo
           enddo
           deallocate(buffer_2d, mask_2d)
        else
           buffer(:) = 0
           mask(:) = .false.
           call assert((ndims .LE. 1), "ndims must be <= 1 when is compressed is false")
           if(has_records) then  ! record variable
              start(:) = 1; nread(:) = 1; nwrite(:) = 1
              start(1) = tlev
              __NF_ASRT__(nfu_get_compressed_var_r8n(in_ncid,varname,buffer,mask,start=start,count=nread))
                 do n = 1, nfiles_out
                    __NF_ASRT__(nfu_put_vara_r8(out_ncid(n),varname,start, nwrite,buffer))
                 enddo
           else
              __NF_ASRT__(nfu_get_compressed_var_r8n(in_ncid,varname,buffer,mask))
              if (count(mask)>0) then
                 do n = 1, nfiles_out
                    __NF_ASRT__(nfu_put_var_r8(out_ncid(n),varname,pack(buffer,mask)))
                 enddo
              endif
           endif
        endif
        deallocate(buffer,mask)
     enddo
  enddo
  __NF_ASRT__(nf_close(in_ncid))

  do n = 1, nfiles_out
     __NF_ASRT__(nf_close(out_ncid(n)))
  enddo

contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

! ---- parses command line arguments, getting the input file name (single file).
! NOTE: updates global variables.
subroutine parse_command_line()
  character(PATH_MAX) :: arg, param

  integer :: nargs     ! number of command-line arguments
  logical :: do_interpret_arguments
  integer :: i, iostat

  !integer, external :: iargc

  nargs = iargc()
  if(nargs==0) then
     call usage()
     call exit(1)
  endif

  npex_io = 0
  npey_io = 0

  allocate(files(nargs))  ! allocate storage for all file names
  do_interpret_arguments = .true.
  i=1        ! counter of all command-line arguments
  nfiles = 0 ! counter of input files
  do while (i<=nargs)
     call getarg(i,arg)
     if(debug>1) write(*,*)'argument ',i, trim(arg)
     if(arg(1:1)=='-'.and.do_interpret_arguments) then
        select case(trim(arg))
        case('--')
           do_interpret_arguments = .false.

        case('-D','--debug-level')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           call getarg(i+1,param)
           read(param,*,iostat=iostat) debug
           call assert(iostat==0,trim(arg)//' flag must be followed by integer verbosity level')
           i=i+1
        case('-i')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           call getarg(i+1,param)
           read(param,*,iostat=iostat) npex_io
           call assert(iostat==0,trim(arg)//' flag must be followed by integer verbosity level')
           i=i+1
        case('-j')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           call getarg(i+1,param)
           read(param,*,iostat=iostat) npey_io
           call assert(iostat==0,trim(arg)//' flag must be followed by integer verbosity level')
           i=i+1
        case ('-h','-?','--help')
           call usage()
           call exit(1)

        case default
           call usage()
           call assert(.false.,'argument "'//trim(arg)//'" is illegal')
        end select
     else
        ! argument is input or output file
        nfiles = nfiles+1
        files(nfiles) = arg
     endif
     i = i+1
  enddo
  if (nfiles==1) then
     infile = files(1)
  else if(nfiles==0) then
     call assert(.false.,'infile is not specified')
  else
     call assert(.false.,'number of files specified must be 1')
  endif

end subroutine


! ---- prints usage information
subroutine usage()
  character(len=PATH_MAX) :: name
  call getarg(0,name)
  write(*,'(a)')'Scatters one file into several distributed file.'
  write(*,'(a)')'Normally used to scatter bombined lm3 restart file'
  write(*,'(a)')'The output files name is in.nc.????'
  write(*,'(a)')
  write(*,'(a)')'Usage:'
  write(*,'(a)')'  '//trim(name)//' [-D debug-level] -i ndiv_x -j ndiv_y in.nc '
  write(*,'(a)')
  write(*,'(a)')'-D debug-level   Specifies level of debug output verbosity'
  write(*,'(a)')'-i ndiv_x        Specifies number of divisions in x-direction (Same as io_layout(1))'
  write(*,'(a)')'-j ndiv_y        Specifies number of divisions in y-direction (Same as io_layout(2))'
  write(*,'(a)')'in.nc            Input file name'
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
     call exit(1)
  endif
end subroutine

end program combine_res
