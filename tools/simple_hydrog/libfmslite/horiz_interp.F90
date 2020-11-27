!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module horiz_interp_mod
!> \author Zhi Liang
!! \author Bruce Wyman
!!
!! \brief Performs spatial interpolation between grids.
!!
!! This module can interpolate data from any logically rectangular grid
!! to any logically rectangular grid. Four interpolation schems are used here:
!! conservative, bilinear, bicubic and inverse of square distance weighted.
!! The four interpolation schemes are implemented seperately in
!! horiz_interp_conserver_mod, horiz_interp_blinear_mod, horiz_interp_bicubic_mod
!! and horiz_interp_spherical_mod. bicubic interpolation requires the source grid
!! is regular lon/lat grid. User can choose the interpolation method in the
!! public interface horiz_interp_new through optional argument interp_method,
!! with acceptable value "conservative", "bilinear", "bicubic" and "spherical".
!! The default value is "conservative". There is an optional mask field for
!! missing input data. An optional output mask field may be used in conjunction with
!! the input mask to show where output data exists.

  use mpp_mod,                    only: mpp_error, FATAL, WARNING
  use constants_mod,              only: PI
  use horiz_interp_type_mod,      only: horiz_interp_type, assignment(=)
  use horiz_interp_type_mod,      only: BILINEAR
  use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_init, horiz_interp_bilinear
  use horiz_interp_bilinear_mod,  only: horiz_interp_bilinear_new, horiz_interp_bilinear_del

  implicit none
  private

  public horiz_interp_type, horiz_interp, horiz_interp_new, horiz_interp_del
  public horiz_interp_init, horiz_interp_end, assignment(=)

  ! <INTERFACE NAME="horiz_interp_new">
  !   <OVERVIEW>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights.
  !   </OVERVIEW>
  !   <DESCRIPTION>
  !      Allocates space and initializes a derived-type variable
  !      that contains pre-computed interpolation indices and weights
  !      for improved performance of multiple interpolations between
  !      the same grids. This routine does not need to be called if you
  !      are doing a single grid-to-grid interpolation.
  !   </DESCRIPTION>
  !   <IN NAME="lon_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
  !      Longitude (in radians) for source data grid. You can pass 1-D lon_in to
  !      represent the geographical longitude of regular lon/lat grid, or just
  !      pass geographical longitude(lon_in is 2-D). The grid location may be
  !      located at grid cell edge or center, decided by optional argument "grid_at_center".
  !   </IN>
  !   <IN NAME="lat_in" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians">
  !      Latitude (in radians) for source data grid. You can pass 1-D lat_in to
  !      represent the geographical latitude of regular lon/lat grid, or just
  !      pass geographical latitude(lat_in is 2-D). The grid location may be
  !      located at grid cell edge or center, decided by optional argument "grid_at_center".
  !   </IN>
!   <IN NAME="lon_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Longitude (in radians) for destination data grid. You can pass 1-D lon_out to
!      represent the geographical longitude of regular lon/lat grid, or just
!      pass geographical longitude(lon_out is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="lat_out" TYPE="real" DIM="dimension(:), dimension(:,:)" UNITS="radians" >
!      Latitude (in radians) for destination data grid. You can pass 1-D lat_out to
!      represent the geographical latitude of regular lon/lat grid, or just
!      pass geographical latitude(lat_out is 2-D). The grid location may be
!      located at grid cell edge or center, decided by optional argument "grid_at_center".
!   </IN>
!   <IN NAME="verbose" TYPE="integer">
!      Integer flag that controls the amount of printed output.
!      verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="interp_method" TYPE="character(len=*)" >
!      interpolation method, = "conservative", using conservation scheme,
!      = "bilinear", using bilinear interpolation, = "spherical",using spherical regrid.
!      = "bicubic", using bicubic interpolation. The default value is "convervative".
!   </IN>
!   <IN NAME = "src_modulo" >
!      Indicate the source data grid is cyclic or not.
!   </IN>
!   <IN NAME = "grid_at_center" >
!      Indicate the data is on the center of grid box or the edge of grid box.
!      When true, the data is on the center of grid box. default vaule is false.
!      This option is only available when interp_method = "bilinear" or "bicubic".
!   </IN>
!   <OUT NAME="Interp" >
!      A derived-type variable containing indices and weights used for subsequent
!      interpolations. To reinitialize this variable for a different grid-to-grid
!      interpolation you must first use the "horiz_interp_del" interface.
!   </OUT>

interface horiz_interp_new
  module procedure horiz_interp_new_1d     ! Source grid is 1d, destination grid is 1d
  module procedure horiz_interp_new_1d_src ! Source grid is 1d, destination grid is 2d
  module procedure horiz_interp_new_2d     ! Source grid is 2d, destination grid is 2d
  module procedure horiz_interp_new_1d_dst ! Source grid is 2d, destination grid is 1d
end interface
! </INTERFACE>

! <INTERFACE NAME="horiz_interp">
!
!   <OVERVIEW>
!     Subroutine for performing the horizontal interpolation between two grids.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Subroutine for performing the horizontal interpolation between
!     two grids. There are two forms of this interface.
!     Form A requires first calling horiz_interp_new, while Form B
!     requires no initialization.
!   </DESCRIPTION>

!   <IN NAME="Interp" >
!     Derived-type variable containing interpolation indices and weights.
!     Returned by a previous call to horiz_interp_new.
!   </IN>
!   <IN NAME="data_in">
!      Input data on source grid.
!   </IN>
!   <IN NAME="verbose">
!      flag for the amount of print output.
!               verbose = 0, no output; = 1, min,max,means; = 2, still more
!   </IN>
!   <IN NAME="mask_in">
!      Input mask, must be the same size as the input data. The real value of
!      mask_in must be in the range (0.,1.). Set mask_in=0.0 for data points
!      that should not be used or have missing data. It is Not needed for
!      spherical regrid.
!   </IN>
!   <IN NAME="missing_value" >
!      Use the missing_value to indicate missing data.
!   </IN>
!   <IN NAME="missing_permit">
!      numbers of points allowed to miss for the bilinear interpolation. The value
!      should be between 0 and 3.
!   </IN>
!   <IN NAME="lon_in, lat_in" >
!      longitude and latitude (in radians) of source grid. More explanation can
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <IN NAME="lon_out, lat_out" >
!      longitude and latitude (in radians) of destination grid. More explanation can
!      be found in the documentation of horiz_interp_new.
!   </IN>
!   <OUT NAME="data_out">
!      Output data on destination grid.
!   </OUT>
!   <OUT NAME="mask_out">
!      Output mask that specifies whether data was computed.
!   </OUT>

!   <ERROR MSG="size of input array incorrect" STATUS="FATAL">
!      The input data array does not match the size of the input grid edges
!      specified. If you are using the initialization interface make sure you
!      have the correct grid size.
!   </ERROR>
!   <ERROR MSG="size of output array incorrect" STATUS="FATAL">
!      The output data array does not match the size of the input grid
!      edges specified. If you are using the initialization interface make
!      sure you have the correct grid size.
!   </ERROR>

interface horiz_interp
  module procedure horiz_interp_base_2d
  module procedure horiz_interp_base_3d
  module procedure horiz_interp_solo_1d
  module procedure horiz_interp_solo_1d_src
  module procedure horiz_interp_solo_2d
  module procedure horiz_interp_solo_1d_dst
  module procedure horiz_interp_solo_old
end interface
! </INTERFACE>

logical :: module_is_initialized = .FALSE.

contains

  !> Initializes the module
  !!
  !! Initialize the horiz_interp module, this is essentillay a no-op.  This routine
  !! is to be compatiable with libFMS.
  subroutine horiz_interp_init
    module_is_initialized = .true.
  end subroutine horiz_interp_init

!#######################################################################
!  <SUBROUTINE NAME="horiz_interp_new_1d" INTERFACE="horiz_interp_new">
!  <IN NAME="lon_in" TYPE="real" DIM="(:),(:,:)" UNITS="radians"></IN>
!  <IN NAME="lat_in" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lon_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="lat_out" TYPE="real" DIM="(:),(:,:)"></IN>
!  <IN NAME="verbose" TYPE="integer, optional"></IN>
!  <IN NAME="interp_method" TYPE="character(len=*),optional"></IN>
!  <IN NAME="src_modulo" TYPE="logical, optional" > </IN>
!  <OUT NAME="Interp" TYPE="type(horiz_interp_type)"></OUT>

!<PUBLICROUTINE INTERFACE="horiz_interp_new">
  subroutine horiz_interp_new_1d(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
      & num_nbrs, max_dist, src_modulo, grid_at_center, mask_in, mask_out)
!</PUBLICROUTINE>
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in), dimension(:) :: lon_in , lat_in
    real, intent(in), dimension(:) :: lon_out, lat_out
    integer, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: interp_method
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    logical, intent(in), optional :: grid_at_center
    real, intent(in), dimension(:,:), optional :: mask_in  ! dummy
    real, intent(inout),dimension(:,:), optional :: mask_out ! dummy

    real, dimension(:,:), allocatable :: lon_src, lat_src, lon_dst, lat_dst
    real, dimension(:), allocatable :: lon_src_1d, lat_src_1d, lon_dst_1d, lat_dst_1d
    integer :: i, j, nlon_in, nlat_in, nlon_out, nlat_out
    logical :: center
    character(len=40) :: method

    call horiz_interp_init

    ! Default method
    method = 'conservative'
    if (present(interp_method)) method = interp_method

    select case (trim(method))
    case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if (present(grid_at_center)) center = grid_at_center
      if (center) then
        nlon_out = size(lon_out(:))
        nlat_out = size(lat_out(:))
        allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
        do i=1, nlon_out
          lon_dst(i,:) = lon_out(i)
        end do
        do j=1, nlat_out
          lat_dst(:,j) = lat_out(j)
        end do

        call horiz_interp_bilinear_new ( Interp, lon_in, lat_in, lon_dst, lat_dst,&
            & verbose, src_modulo)
          deallocate(lon_dst, lat_dst)
      else
        nlon_in  = size(lon_in(:))-1
          nlat_in  = size(lat_in(:))-1
          nlon_out = size(lon_out(:))-1
          nlat_out = size(lat_out(:))-1
          allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
          allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
          do i=1, nlon_in
            lon_src_1d(i) = (lon_in(i)+lon_in(i+1))*0.5
          end do
          do j=1, nlat_in
             lat_src_1d(j) = (lat_in(j)+lat_in(j+1))*0.5
          end do
          do i=1, nlon_out
             lon_dst(i,:) = (lon_out(i)+lon_out(i+1))*0.5
          enddo
          do j=1, nlat_out
             lat_dst(:,j) = (lat_out(j)+lat_out(j+1)) * 0.5
          enddo
          call horiz_interp_bilinear_new(Interp, lon_src_1d, lat_src_1d, lon_dst, lat_dst,&
              & verbose, src_modulo)
          deallocate(lon_src_1d, lat_src_1d, lon_dst, lat_dst)
      end if
    case default
      call mpp_error(FATAL, 'horiz_interp_mod::horiz_interp_new_1d: Unknown interp_method. Should'//&
          & 'be conservative, bilinear, bicubic, spherical')
    end select

    Interp%I_am_initialized = .true.
  end subroutine horiz_interp_new_1d
!  </SUBROUTINE>

  subroutine horiz_interp_new_1d_src(Interp, lon_in, lat_in, lon_out, lat_out, verbose,&
      & interp_method, num_nbrs, max_dist, src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out)
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in), dimension(:) :: lon_in , lat_in
    real, intent(in), dimension(:,:) :: lon_out, lat_out
    integer, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: interp_method
    integer, intent(in), optional :: num_nbrs  ! minimum number of neighbors
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    logical, intent(in), optional :: grid_at_center
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out),dimension(:,:), optional :: mask_out
    logical, intent(in), optional :: is_latlon_out

    real, dimension(:,:), allocatable :: lon_src, lat_src
    real, dimension(:), allocatable :: lon_src_1d, lat_src_1d
    integer :: i, j, nlon_in, nlat_in
    character(len=40) :: method
    logical :: center
    logical :: dst_is_latlon

    call horiz_interp_init

    method = 'bilinear'
    if (present(interp_method)) method = interp_method

    select case (trim(method))
   case ("bilinear")
      Interp%interp_method = BILINEAR
      center = .false.
      if (present(grid_at_center)) center = grid_at_center
      if (center) then
        call horiz_interp_bilinear_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo)
      else
        nlon_in  = size(lon_in(:))-1
        nlat_in  = size(lat_in(:))-1
        allocate(lon_src_1d(nlon_in), lat_src_1d(nlat_in))
        do i=1, nlon_in
          lon_src_1d(i) = (lon_in(i) + lon_in(i+1)) * 0.5
        end do
        do j=1, nlat_in
          lat_src_1d(j) = (lat_in(j) + lat_in(j+1)) * 0.5
        end do
        call horiz_interp_bilinear_new(Interp, lon_src_1d, lat_src_1d, lon_out, lat_out,&
            & verbose, src_modulo )
        deallocate(lon_src_1d,lat_src_1d)
      end if
   case default
      call mpp_error(FATAL,'horiz_interp_mod::horiz_interp_new_1d: Unknown interp_method. Should '//&
          & 'be bilinear')
   end select

   Interp%I_am_initialized = .true.
  end subroutine horiz_interp_new_1d_src

  subroutine horiz_interp_new_2d(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
      & num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in, is_latlon_out  )
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in), dimension(:,:) :: lon_in , lat_in
    real, intent(in), dimension(:,:) :: lon_out, lat_out
    integer, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: interp_method
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out),dimension(:,:), optional :: mask_out
    logical, intent(in), optional :: is_latlon_in, is_latlon_out

    logical :: src_is_latlon, dst_is_latlon
    character(len=40) :: method

    call horiz_interp_init

    method = 'bilinear'
    if (present(interp_method)) method = interp_method

    select case (trim(method))
    case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, src_modulo)
    case default
      call mpp_error(FATAL,'horiz_interp_mod::horiz_interp_new_2d: Unknown interp_method. Should '//&
          & 'be bilinear')
    end select

    Interp%I_am_initialized = .true.
  end subroutine horiz_interp_new_2d

  subroutine horiz_interp_new_1d_dst(Interp, lon_in, lat_in, lon_out, lat_out,&
      & verbose, interp_method, num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in)
    type(horiz_interp_type), intent(inout) :: Interp
    real, intent(in), dimension(:,:) :: lon_in , lat_in
    real, intent(in), dimension(:) :: lon_out, lat_out
    integer, intent(in), optional :: verbose
    character(len=*), intent(in), optional :: interp_method
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out),dimension(:,:), optional :: mask_out
    logical, intent(in), optional :: is_latlon_in

    character(len=40) :: method
    integer :: i, j, nlon_out, nlat_out
    real, dimension(:,:), allocatable :: lon_dst, lat_dst
    logical :: src_is_latlon

    call horiz_interp_init

    method = 'bilinear'
    if (present(interp_method)) method = interp_method

    nlon_out = size(lon_out(:))
    nlat_out = size(lat_out(:))
    allocate(lon_dst(nlon_out,nlat_out), lat_dst(nlon_out,nlat_out))
    do i=1, nlon_out
      lon_dst(i,:) = lon_out(i)
    end do
    do j=1, nlat_out
      lat_dst(:,j) = lat_out(j)
    end do

    select case (trim(method))
    case ("bilinear")
      Interp%interp_method = BILINEAR
      call horiz_interp_bilinear_new(Interp, lon_in, lat_in, lon_dst, lat_dst,&
          & verbose, src_modulo)
    case default
      call mpp_error(FATAL,'when source grid are 2d, interp_method should be bilinear')
    end select

    deallocate(lon_dst,lat_dst)

    Interp%I_am_initialized = .true.
  end subroutine horiz_interp_new_1d_dst

  subroutine horiz_interp_base_2d(Interp, data_in, data_out, verbose, mask_in, mask_out,&
      & missing_value, missing_permit, err_msg, new_missing_handle )
    type (horiz_interp_type), intent(in) :: Interp
    real, intent(in), dimension(:,:) :: data_in
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    character(len=*), intent(out), optional :: err_msg
    logical, intent(in), optional :: new_missing_handle

    if (present(err_msg)) err_msg = ''
    if (.not.Interp%I_am_initialized) then
      call mpp_error(FATAL, 'horiz_interp_mod::horiz_interp_base_2d: The horiz_interp_type '//&
          & 'variable is not initialized')
    end if

    select case(Interp%interp_method)
    case(BILINEAR)
    case default
      call mpp_error(FATAL,'interp_method should be "bilinear".')
    end select

    return
  end subroutine horiz_interp_base_2d
  ! </SUBROUTINE>

  !-----------------------------------------------------------------------
  !   overload of interface horiz_interp_base_2d
  !   uses 3d arrays for data and mask
  !   this allows for multiple interpolations with one call
  !-----------------------------------------------------------------------
  subroutine horiz_interp_base_3d(Interp, data_in, data_out, verbose, mask_in, mask_out,&
      & missing_value, missing_permit, err_msg)
    type (horiz_interp_type), intent(in) :: Interp
    real, intent(in), dimension(:,:,:) :: data_in
    real, intent(out), dimension(:,:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:,:), optional :: mask_in
    real, intent(out), dimension(:,:,:), optional :: mask_out
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    character(len=*), intent(out), optional :: err_msg

    integer :: n

    if (present(err_msg)) err_msg = ''
    if (.not.Interp%I_am_initialized) then
      call mpp_error('horiz_interp_mod::horiz_interp_base_3d',&
          & 'The horiz_interp_type variable is not initialized',&
          & FATAL)
      return
    end if

    do n=1, size(data_in,3)
      if (present(mask_in))then
        if (present(mask_out)) then
          call horiz_interp_base_2d(Interp, data_in(:,:,n), data_out(:,:,n), verbose,&
              & mask_in(:,:,n), mask_out(:,:,n), missing_value, missing_permit)
        else
          call horiz_interp_base_2d(Interp, data_in(:,:,n), data_out(:,:,n), verbose,&
              & mask_in(:,:,n), missing_value=missing_value, missing_permit=missing_permit)
        end if
      else
        if (present(mask_out)) then
          call horiz_interp_base_2d(Interp, data_in(:,:,n), data_out(:,:,n), verbose,&
              & mask_out=mask_out(:,:,n), missing_value=missing_value, missing_permit=missing_permit)
        else
          call horiz_interp_base_2d(Interp, data_in(:,:,n), data_out(:,:,n), verbose,&
              & missing_value=missing_value, missing_permit=missing_permit)
        end if
      end if
    end do
    return
  end subroutine horiz_interp_base_3d

  ! interpolates from a rectangular grid to rectangular grid.
  ! interp_method can be the value conservative, bilinear or spherical.
  ! horiz_interp_new don't need to be called before calling this routine.
  !<PUBLICROUTINE INTERFACE="horiz_interp">
  subroutine horiz_interp_solo_1d(data_in, lon_in, lat_in, lon_out, lat_out, data_out, verbose,&
      & mask_in, mask_out, interp_method, missing_value, missing_permit, num_nbrs, max_dist,&
      & src_modulo, grid_at_center)
  !</PUBLICROUTINE>
    real, intent(in), dimension(:,:) :: data_in
    real, intent(in), dimension(:) :: lon_in , lat_in
    real, intent(in), dimension(:) :: lon_out, lat_out
    real, intent(out),dimension(:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    character(len=*), intent(in), optional :: interp_method
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    logical, intent(in), optional :: grid_at_center

    type (horiz_interp_type) :: Interp

    call horiz_interp_init
    call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
        & num_nbrs, max_dist, src_modulo, grid_at_center)
    call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit)
    call horiz_interp_del(Interp)
  end subroutine horiz_interp_solo_1d

  ! interpolates from a uniformly spaced grid to any output grid.
  ! interp_method can be the value "onservative","bilinear" or "spherical".
  ! horiz_interp_new don't need to be called before calling this routine.
  subroutine horiz_interp_solo_1d_src ( data_in, lon_in, lat_in, lon_out, lat_out, data_out,&
      & verbose, mask_in, mask_out, interp_method, missing_value, missing_permit, num_nbrs,&
      & max_dist, src_modulo, grid_at_center)
    real, intent(in), dimension(:,:) :: data_in
    real, intent(in), dimension(:) :: lon_in , lat_in
    real, intent(in), dimension(:,:) :: lon_out, lat_out
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    character(len=*), intent(in), optional :: interp_method
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo
    logical, intent(in), optional :: grid_at_center

    type(horiz_interp_type) :: Interp
    logical :: dst_is_latlon
    character(len=128) :: method

    call horiz_interp_init
    method = 'conservative'
    if (present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    if (trim(method) == 'conservative') dst_is_latlon = is_lat_lon(lon_out, lat_out)

    if (dst_is_latlon) then
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, grid_at_center, is_latlon_out=dst_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit)
    else
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, grid_at_center, mask_in, mask_out, is_latlon_out=dst_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, missing_value=missing_value, missing_permit=missing_permit )
    end if
    call horiz_interp_del ( Interp )
  end subroutine horiz_interp_solo_1d_src

  ! interpolates from any grid to any grid. interp_method should be "spherical"
  ! horiz_interp_new don't need to be called before calling this routine.
  subroutine horiz_interp_solo_2d(data_in, lon_in, lat_in, lon_out, lat_out, data_out, verbose,&
      & mask_in, mask_out, interp_method, missing_value, missing_permit, num_nbrs, max_dist, src_modulo)
    real, intent(in), dimension(:,:) :: data_in
    real, intent(in), dimension(:,:) :: lon_in , lat_in
    real, intent(in), dimension(:,:) :: lon_out, lat_out
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    character(len=*), intent(in), optional :: interp_method
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo

    type(horiz_interp_type) :: Interp
    logical :: dst_is_latlon, src_is_latlon
    character(len=128) :: method

    call horiz_interp_init

    method = 'conservative'
    if (present(interp_method)) method = interp_method
    dst_is_latlon = .true.
    src_is_latlon = .true.
    if (trim(method) == 'conservative') then
      dst_is_latlon = is_lat_lon(lon_out, lat_out)
      src_is_latlon = is_lat_lon(lon_in, lat_in)
    end if

    if (dst_is_latlon .and. src_is_latlon) then
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, is_latlon_in=dst_is_latlon, is_latlon_out=dst_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit)
    else
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in=dst_is_latlon, is_latlon_out=dst_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, missing_value=missing_value, missing_permit=missing_permit)
    end if
    call horiz_interp_del ( Interp )
  end subroutine horiz_interp_solo_2d

!   interpolates from any grid to rectangular longitude/latitude grid.
!   interp_method should be "spherical".
!   horiz_interp_new don't need to be called before calling this routine.
  subroutine horiz_interp_solo_1d_dst(data_in, lon_in, lat_in, lon_out, lat_out, data_out, verbose,&
      & mask_in, mask_out,interp_method,missing_value, missing_permit, num_nbrs, max_dist, src_modulo)
    real, intent(in), dimension(:,:) :: data_in
    real, intent(in), dimension(:,:) :: lon_in , lat_in
    real, intent(in), dimension(:) :: lon_out, lat_out
    real, intent(out), dimension(:,:) :: data_out
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out
    character(len=*), intent(in), optional :: interp_method
    real, intent(in), optional :: missing_value
    integer, intent(in), optional :: missing_permit
    integer, intent(in), optional :: num_nbrs
    real, intent(in), optional :: max_dist
    logical, intent(in), optional :: src_modulo

    type(horiz_interp_type) :: Interp
    logical :: src_is_latlon
    character(len=128) :: method

    call horiz_interp_init

    method = 'conservative'
    if (present(interp_method)) method = interp_method
    src_is_latlon = .true.
    if (trim(method) == 'conservative') src_is_latlon = is_lat_lon(lon_in, lat_in)

    if (src_is_latlon) then
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, is_latlon_in=src_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, mask_in, mask_out, missing_value, missing_permit)
    else
      call horiz_interp_new(Interp, lon_in, lat_in, lon_out, lat_out, verbose, interp_method,&
          & num_nbrs, max_dist, src_modulo, mask_in, mask_out, is_latlon_in=src_is_latlon)
      call horiz_interp(Interp, data_in, data_out, verbose, missing_value=missing_value,&
          & missing_permit=missing_permit )
    end if
    call horiz_interp_del ( Interp )
  end subroutine horiz_interp_solo_1d_dst

  !> Overloaded version of interface horiz_interp_solo_2
  subroutine horiz_interp_solo_old(data_in, wb, sb, dx, dy, lon_out, lat_out, data_out, verbose,&
      & mask_in, mask_out)
    real, intent(in), dimension(:,:) :: data_in !< Global input data stored from west to east
        !! (first dimension), south to north (second dimension).
    real, intent(in) :: wb !< Longitude (in radians) that corresponds to western-most boundary of
        !! grid box i=1 in array data_in.
    real, intent(in) :: sb !< Latitude (in radians) that corresponds to southern-most boundary of
        !! grid box j=1 in array data_in.
    real, intent(in) :: dx !< Grid spacing (in radians) for the longitude axis (first dimension) for
        !! the input data.
    real, intent(in) :: dy !< Grid spacing (in radians) for the latitude axis (second dimension) for
        !! the input data.
    real, intent(in), dimension(:) :: lon_out !< The longitude edges (in radians) for output data
        !! grid boxes.  The values are for adjacent grid boxes and must increase in value. If there
        !! are MLON grid boxes there must be MLON+1 edge values.
    real, intent(in), dimension(:) :: lat_out !< The latitude edges (in radians) for output data
        !! grid boxes.  The values are for adjacent grid boxes and may increase or decrease in
        !! value. If there are NLAT grid boxes there must be NLAT+1 edge values.
    real, intent(out), dimension(:,:) :: data_out !< Output data on the output grid defined by grid
        !! box edges: blon_out and blat_out.
    integer, intent(in), optional :: verbose
    real, intent(in), dimension(:,:), optional :: mask_in
    real, intent(out), dimension(:,:), optional :: mask_out

    real, dimension(size(data_in,1)+1) :: blon_in
    real, dimension(size(data_in,2)+1) :: blat_in
    integer :: i, j, nlon_in, nlat_in
    real :: tpi

    call horiz_interp_init

    tpi = 2.*PI
    nlon_in = size(data_in,1)
    nlat_in = size(data_in,2)

    do i=1, nlon_in+1
      blon_in(i) = wb + float(i-1)*dx
    end do
    if (abs(blon_in(nlon_in+1)-blon_in(1)-tpi) < epsilon(blon_in))&
        & blon_in(nlon_in+1)=blon_in(1)+tpi

    do j=2, nlat_in
      blat_in(j) = sb + float(j-1)*dy
    end do
    blat_in(1) = -0.5*pi
    blat_in(nlat_in+1) =  0.5*pi

    call horiz_interp_solo_1d (data_in, blon_in, blat_in, lon_out, lat_out, data_out, verbose,&
        & mask_in, mask_out)
  end subroutine horiz_interp_solo_old

  !>     Deallocates memory used by "horiz_interp_type" variables.
  !!       Must be called before reinitializing with horiz_interp_new.
  subroutine horiz_interp_del(Interp)
    type(horiz_interp_type), intent(inout) :: Interp !< A derived-type variable returned by previous
        !! call to horiz_interp_new. The input variable must have allocated arrays. The returned
        !! variable will contain deallocated arrays.

    select case(Interp % interp_method)
    case (BILINEAR)
      call horiz_interp_bilinear_del(Interp )
    end select

    Interp%I_am_initialized = .false.
  end subroutine horiz_interp_del

  !> Dummy routine.
  subroutine horiz_interp_end
    return
  end subroutine horiz_interp_end

  function is_lat_lon(lon, lat)
    real, dimension(:,:), intent(in) :: lon, lat
    logical :: is_lat_lon
    integer :: i, j, nlon, nlat, num

    is_lat_lon = .true.
    nlon = size(lon,1)
    nlat = size(lon,2)
    LOOP_LAT: do j=1, nlat
      do i=2, nlon
        if (lat(i,j).NE.lat(1,j)) then
          is_lat_lon = .false.
          exit LOOP_LAT
        end if
      end do
    end do LOOP_LAT

    if (is_lat_lon) then
      LOOP_LON: do i=1, nlon
        do j=2, nlat
          if (lon(i,j).NE.lon(i,1)) then
            is_lat_lon = .false.
            exit LOOP_LON
          end if
        end do
      end do LOOP_LON
    end if

    num = 0
    if (is_lat_lon) num = 1
    if (num == 1) then
      is_lat_lon = .true.
    else
      is_lat_lon = .false.
    end if

    return
  end function is_lat_lon
end module horiz_interp_mod
