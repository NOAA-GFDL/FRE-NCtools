!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License. See the file COPYING contained in this directory
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

  type dimtype
     character(NF_MAX_NAME) :: name = ''
     integer :: len    = -1 ! length in output file
     logical :: compressed = .FALSE.
     ! the rest makes sense only for compressed dimensions
     integer :: buflen = -1 ! sum of lengths in input files, size of input buffer
  end type

  character(PATH_MAX), allocatable :: files(:) ! names of all files on the command line
  character(PATH_MAX)              :: outfile  ! name of the output file
  integer :: nfiles    ! number of files on command line
  integer :: verbosity = 0 ! debug output verbosity level
  integer, allocatable :: input(:)             ! netcdf IDs of input files
  integer :: ncid,dimid,varid,varid1,ovarid,dimlen,ndims,nvars,ngatts,xtype
  integer :: dimids(NF_MAX_DIMS), start(NF_MAX_DIMS), cnt(NF_MAX_DIMS)
  integer :: dimlens(NF_MAX_DIMS)
  logical :: has_records
  integer :: in_format ! format of input files
  integer :: cmode     ! mode for output file creation
  character(NF_MAX_NAME) :: varname,attname
  real     , allocatable :: buffer(:), obuffer(:)
  character, allocatable :: text(:)
  integer  , allocatable :: rank(:) ! re-ordering of the indices
  integer  , allocatable :: sizes(:) ! length of compressed dimension for each of the files
  integer :: nrec, recsize
  integer :: i,j,k,n,rec,ii,cdim,ifile
  type(dimtype), allocatable :: dim(:)


  ! get command line options and list of files
  call parse_command_line() ! modigies global data!

  call assert(nfiles>0,'at least one input file must be specified')
  if(verbosity>0) then
     do i = 1,nfiles
        write(*,'("input file",i4," : ",a)')i, '"'//trim(files(i))//'"'
     enddo
     write(*,'("output file : ",a)')'"'//trim(outfile)//'"'
  endif

  ! open all input files and determine the creation mode of output file:
  ! if any of the input files is 64-bit then the output is 64-bit as well,
  ! otherwise it's 32-bit
  allocate(input(nfiles))
  do i = 1,nfiles
     __NF_ASRT__(nf_open(files(i),NF_NOWRITE,input(i)))
     __NF_ASRT__(nf_inq_format(input(i),in_format))
  enddo

  if (in_format==NF_FORMAT_NETCDF4) then
     cmode = NF_NETCDF4
  elseif (in_format==NF_FORMAT_NETCDF4_CLASSIC) then
     cmode=IOR(NF_NETCDF4,NF_CLASSIC_MODEL)
  elseif (in_format==NF_FORMAT_64BIT) then
     cmode=IOR(NF_CLOBBER,NF_64BIT_OFFSET)
     if(verbosity>0)write(*,'("output file is 64-bit netcdf")')
  elseif (in_format==NF_FORMAT_CLASSIC) then
     cmode=IOR(NF_CLOBBER,NF_CLASSIC_MODEL)
     if(verbosity>0)write(*,'("output file is 32-bit netcdf")')
  else
     call assert(.false.,'Unknown netCDF format')
  endif

  ! create output file

  ! mpp_io supports environment variables to set these. For Riga, we'll simply use the defaults`

  __NF_ASRT__(nf__create(outfile,cmode,0,blksz,ncid))

  ! Create netcdf structure in the output NetCDF file, using last input file
  ! as a template.

  ! clone all dimensions; for compressed dimensions calculate the total length in all
  ! input files
  __NF_ASRT__(nf_inq_ndims(input(nfiles),ndims))
  allocate(dim(ndims))
  do dimid = 1,ndims
     associate(d => dim(dimid))
     __NF_ASRT__(nfu_inq_dim(input(nfiles),dimid,dimname=d%name))
     call inquire_dimension(input(:), d%name, &
         len=dimlen, siz=d%buflen, compressed=d%compressed, is_unlim=has_records)
     d%len = dimlen

     ! TODO: check that there are no duplicate values in the index, make it a warning

     ! TODO: check that the total size is > 0
     ! can't have 0-length dimension, since it is (mis-)understood by netcdf as
     ! a record one.

     if(has_records)then
        dimlen = NF_UNLIMITED
     else
        dimlen = max(dimlen,1)
     endif
     if(verbosity>0)&
           write(*,'(x,a,i8)')'defining dimension "'//trim(d%name)//'" with length',d%len
     __NF_ASRT__(nf_def_dim(ncid,d%name,dimlen,i)) ! i is just a dummy var for dimid, unused
     end associate
  enddo

  ! clone all variable definitions
  __NF_ASRT__(nf_inq_nvars(input(nfiles),nvars))
  do i = 1,nvars
     __NF_ASRT__(nfu_clone_var(input(nfiles),i,ncid))
     ! NOTE: since cloning of variable definition relies on dimension names,
     ! each variable tile and compressed dimensions automaticaly get the right
     ! size, as defined while creating dimensions in the output file
  enddo

  ! clone all global attributes
  __NF_ASRT__(nf_inq_natts(input(nfiles),ngatts))
  do i = 1,ngatts
     __NF_ASRT__(nf_inq_attname(input(nfiles),NF_GLOBAL,i,attname))
     __NF_ASRT__(nf_copy_att(input(nfiles),NF_GLOBAL,attname,ncid,NF_GLOBAL))
  enddo

  ! ---- end of definition stage
  __NF_ASRT__(nf__enddef(ncid,HEADERPAD,4,0,4))

  ! copy all uncompressed vars
  do varid = 1, nvars
     __NF_ASRT__(nfu_inq_var(input(nfiles),varid,name=varname,ndims=ndims,dimids=dimids,dimlens=dimlens,recsize=recsize,nrec=nrec,xtype=xtype))
     n = 0
     do k = 1,ndims
        if (dim(dimids(k))%compressed) n = n+1
     enddo
     if (n>1) then
        write(*,*)'Variable "'//trim(varname)//'" has more then one compressed dimension. Cannot handle that.'
        call exit(255)
     endif
     if (n==0) then ! no compressed dims => variable is uncompressed
        if(verbosity>0) write(*,'(2x,a)')'copy uncompressed variable "'//trim(varname)//'"'
        if (xtype==NF_CHAR) then
           ! we are not bothering with writing CHAR variables by record since (1)
           ! they are relatively small (2) they are unlikely to have record
           ! dimension anyway, and (3) there is no convenient interface (yet) in
           ! nfu utilities for by-record i/o for CHAR variables.
           allocate(text(recsize*nrec))
           __NF_ASRT__(nf_get_var_text(input(nfiles),varid,text))
           __NF_ASRT__(nfu_inq_var(ncid,varname,id=varid1))
           __NF_ASRT__(nf_put_var_text(ncid,varid1,text))
           deallocate(text)
        else
           allocate(buffer(recsize))
           do rec=1,nrec
              __NF_ASRT__(nfu_get_rec_r8(input(nfiles),varid,rec,buffer))
              __NF_ASRT__(nfu_put_rec_r8(ncid,varname,rec,buffer))
           enddo
           deallocate(buffer)
        endif
     endif
  enddo

  ! copy compressed vars
  allocate(sizes(nfiles))
  do dimid = 1,size(dim)
     associate(d=>dim(dimid))
     if (.not.d%compressed) cycle ! skip nonn-compressed dimensions

     if (d%len==0) then
        ! if compressed dimension has zero length, then the variables
        ! that depend on it have no values; by skipping them we leave
        ! the output filled with respective _FillValues
        if (verbosity>0) then
           write(*,'(x,a)')'compressed dimension "'//trim(d%name)//'" has no values, skipping it and all variables that depend on it'
           write(*,'(x,a)')'variables in output file will be filled with respective _FillValues'
        endif
        cycle
     endif
     if (verbosity>0) write(*,'(x,a)')'processing compressed dimension "'//trim(d%name)//'"'

     ! get the size of compressed dim in every file
     call inquire_dimension(input(:),d%name,sizes=sizes)
     ! allocate i/o and reordering buffers
     allocate(rank(d%buflen), buffer(d%buflen), obuffer(d%len))
     ! create re-ordering index
     call reorder_compressed_index(input,d%name,rank)

     ! process all vars that depend on this compressed dimension
     do varid = 1, nvars
        __NF_ASRT__(nfu_inq_var(input(nfiles),varid,name=varname,ndims=ndims,dimids=dimids,dimlens=dimlens))
        if (.not.any(dimids(1:ndims)==dimid)) cycle ! skip variables that do not depend on our compressed dim
        if(verbosity>0) write(*,'(2x,a)')'copy compressed variable "'//trim(varname)//'"'

        ! get the output variable ID
        __NF_ASRT__(nfu_inq_var(ncid,varname,id=ovarid))

        ! find index of the compressed dimension
        cdim = 1
        do while(dimids(cdim)/=dimid)
           cdim = cdim+1
        enddo

        ! loop over all uncompressed dimensions
        cnt(:) = 1
        do i = 1, product(dimlens(1:ndims))/dimlens(cdim)
           ! define starting indices for current slice
           ii = i-1
           do j = 1,ndims
              if (j==cdim) then
                 start(j) = 1
              else
                 start(j) = modulo(ii,dimlens(j))+1; ii = ii/dimlens(j)
              endif
           enddo

           ! read slice from the variable in all input files
           k = 1
           do ifile = 1,nfiles
              cnt(cdim) = sizes(ifile)
              __NF_ASRT__(nfu_inq_var(input(ifile),varname,id=varid1))
              __NF_ASRT__(nf_get_vara_double(input(ifile),varid1,start,cnt,buffer(k)))
              k = k+sizes(ifile)
           enddo
           ! reshuffle variable values in desired order
           do k = 1,size(obuffer)
              obuffer(k) = buffer(rank(k))
           enddo
           ! write slice to output file
           cnt(cdim) = d%len
           __NF_ASRT__(nf_put_vara_double(ncid,ovarid,start,cnt,obuffer))
        enddo
     enddo
     deallocate(rank, buffer, obuffer)
     end associate
  enddo

  i = nf_close(ncid)
contains ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

subroutine inquire_dimension(input,name,len,compressed,siz,sizes,is_unlim)
  integer, intent(in) :: input(:) ! netcdf IDs of input files
  character(*), intent(in) :: name ! name of the dimension
  integer, intent(out), optional :: len       ! size of dimension in output file -- may be
  logical, intent(out), optional :: compressed
  integer, intent(out), optional :: siz      ! size of dimension in input file(s)
           ! for un-compressed dimension it is equal to len; for compressed
           ! dimensions it is the sum of all the sizes in all files. It may be
           ! different than len, because we exclude invalid indices in output.
           ! Invalid indices come from compute domains that have no data.
  integer, intent(out), optional :: sizes(:)  ! size of dimension in each of input files
  logical, intent(out), optional :: is_unlim

  integer :: nfiles
  integer :: len_
  logical :: compressed_
  integer :: sizes_(size(input))
  integer, allocatable :: buff(:)
  integer :: i

  nfiles = size(input)
  __NF_ASRT__(nfu_inq_dim(input(nfiles), name, dimlen=len_, is_unlim=is_unlim))
  if (present(len)) len = len_
  if (present(siz)) siz = len_
  if (present(sizes)) sizes(:)=len_

  compressed_ = (nfu_inq_att(input(nfiles), name,'compress')==NF_NOERR)
  if (present(compressed)) compressed = compressed_

  ! ---- the rest is only for compressed data ---------------------------------
  if (.not.compressed_) return

  ! calculate total size of the compressed dimension across all files
  do i=1,nfiles
     __NF_ASRT__(nfu_inq_dim(input(i),name,dimlen=sizes_(i)))
  enddo
  if (present(sizes)) sizes(:) = sizes_(:)
  if (present(siz)) siz = sum(sizes_(:))
  if (present(len)) then
     len = 0
     allocate(buff(maxval(sizes_(:))))
     do i=1,nfiles
        __NF_ASRT__(nfu_get_var_int(input(i),name,buff))
        len = len + count(buff(1:sizes_(i))>=0)
     enddo
     deallocate(buff)
  endif
end subroutine inquire_dimension


subroutine reorder_compressed_index(input,name,rank)
  integer, intent(in) :: input(:) ! netcdf IDs of input files
  character(*), intent(in) :: name
  integer, intent(out) :: rank(:) ! re-ordering index

  integer, allocatable :: buff(:)
  integer :: i, k, dimlen

  allocate(buff(size(rank)))
  k = 1
  do i=1,nfiles
     __NF_ASRT__(nfu_get_var_int(input(i),name,buff(k:)))
     __NF_ASRT__(nfu_inq_dim(input(i),name,dimlen=dimlen))
     k = k+dimlen
  enddo

  ! rank dimension index for re-ordering
  call rank_ascending(buff,rank)

  ! skip leading negatives: they are artifacts of the compute domains with no tiles
  do k = 1,size(rank)
     if (buff(rank(k))>=0) exit ! from loop
  enddo
!  if(k>1) then
     do i = 1,size(rank)-k+1
        rank(i) = rank(i+k-1)
     enddo
     do i = size(rank)-k+2,size(rank)
        rank(i) = -1
     enddo
!  endif
   deallocate(buff)
end subroutine reorder_compressed_index

! ======================================================================
! ranks array x in descending order: on return, idx() contains indices
! of elements of array x in descending order of x values
subroutine rank_ascending(x,idx)
   integer, intent(in)  :: x(:)
   integer, intent(out) :: idx(:)

   integer :: i,n
   integer, allocatable :: t(:)

   n = size(x)
   do i = 1,n
      idx(i) = i
   enddo

   allocate(t((n+1)/2))
   call mergerank(x,idx,n,t)
   deallocate(t)
end subroutine


! =====================================================================
! based on:
! http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
recursive subroutine mergerank(x,a,n,t)
  integer, intent(in) :: n
  integer, intent(in) :: x(*)
  integer, dimension(n), intent(inout) :: a
  integer, dimension((n+1)/2), intent (out) :: t

  integer :: na,nb
  integer :: v

  if (n < 2) return
  if (n == 2) then
     if ( x(a(1)) > x(a(2)) ) then
        v = a(1) ; a(1) = a(2) ; a(2) = v
     endif
     return
  endif
  na=(n+1)/2
  nb=n-na

  call mergerank(x,a,na,t)
  call mergerank(x,a(na+1),nb,t)

  if (x(a(na)) > x(a(na+1))) then
     t(1:na)=a(1:na)
     call merge(x,t,na,a(na+1),nb,a,n)
  endif
end subroutine mergerank

subroutine merge(x,a,na,b,nb,c,nc)
   integer, intent(in) :: na,nb,nc ! Normal usage: NA+NB = NC
   integer, intent(in)    :: x(*)
   integer, intent(in)    :: a(na)    ! B overlays C(NA+1:NC)
   integer, intent(in)    :: b(nb)
   integer, intent(inout) :: c(nc)

   integer :: i,j,k

   i = 1; j = 1; k = 1;
   do while(i <= na .and. j <= nb)
      if (x(a(i)) <= x(b(j))) then
         c(k) = a(i) ; i = i+1
      else
         c(k) = b(j) ; j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= na)
      c(k) = a(i) ; i = i + 1 ; k = k + 1
   enddo
end subroutine merge

! ---- parses command line arguments, getting options and gathering list of
! file names
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

  allocate(files(nargs))  ! allocate storage for all file names
  do_interpret_arguments = .true.
  i=1        ! counter of all command-line arguments
  nfiles = 0 ! counter of input files
  do while (i<=nargs)
     call getarg(i,arg)
     if(verbosity>1) write(*,*)'argument ',i, trim(arg)
     if(arg(1:1)=='-'.and.do_interpret_arguments) then
        select case(trim(arg))
        case('--')
           do_interpret_arguments = .false.

        case('-v','--verbose')
           call assert(i<nargs,trim(arg)//' flag must be followed by integer verbosity level')
           verbosity = verbosity+1

        case ('-h','-?','--help')
           call usage()
           call exit(1)

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
  if(verbosity>0) &
       write(*,*) nfiles, ' input files'
end subroutine


! ---- prints usage information
subroutine usage()
  character(len=PATH_MAX) :: name
  call getarg(0,name)
  write(*,'(a)')'Combines several compressed-by-gathering netcdf files into one.'
  write(*,'(a)')'Normally used to combine lm3 restarts generated by each processor.'
  write(*,'(a)')
  write(*,'(a)')'Usage:'
  write(*,'(a)')'  '//trim(name)//' [-v verbosity-level] in.nc [...] out.nc'
  write(*,'(a)')
  write(*,'(a)')'-v verbosity-level   Specifies level of verbosity output verbosity'
  write(*,'(a)')'in.nc                Input file name(s)'
  write(*,'(a)')'out.nc               Output file name'
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
