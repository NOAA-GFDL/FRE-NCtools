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
module mpp_mod
  use ISO_FORTRAN_ENV

  implicit none
  private

  integer, public, parameter :: NOTE=0, WARNING=1, FATAL=2

  public :: mpp_error

  interface mpp_error
    module procedure mpp_error_basic
    module procedure mpp_error_mesg
    module procedure mpp_error_noargs
    module procedure mpp_error_is
    module procedure mpp_error_rs
    module procedure mpp_error_ia
    module procedure mpp_error_ra
    module procedure mpp_error_ia_ia
    module procedure mpp_error_ia_ra
    module procedure mpp_error_ra_ia
    module procedure mpp_error_ra_ra
    module procedure mpp_error_ia_is
    module procedure mpp_error_ia_rs
    module procedure mpp_error_ra_is
    module procedure mpp_error_ra_rs
    module procedure mpp_error_is_ia
    module procedure mpp_error_is_ra
    module procedure mpp_error_rs_ia
    module procedure mpp_error_rs_ra
    module procedure mpp_error_is_is
    module procedure mpp_error_is_rs
    module procedure mpp_error_rs_is
    module procedure mpp_error_rs_rs
  end interface

  interface array_to_char
    module procedure iarray_to_char
    module procedure rarray_to_char
  end interface

contains

  function iarray_to_char(iarray) result(string)
    integer, intent(in) :: iarray(:)
    character(len=256) :: string
    character(len=32)  :: chtmp
    integer :: i, len_tmp, len_string

    string = ''
    do i=1,size(iarray)
      write(chtmp,'(i16)') iarray(i)
      chtmp = adjustl(chtmp)
      len_tmp = len_trim(chtmp)
      len_string  = len_trim(string)
      string(len_string+1:len_string+len_tmp) = trim(chtmp)
      string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
    end do
    len_string = len_trim(string)
    string(len_string:len_string) = ' ' ! remove trailing comma
  end function iarray_to_char

  function rarray_to_char(rarray) result(string)
    real, intent(in) :: rarray(:)
    character(len=256) :: string
    character(len=32)  :: chtmp
    integer :: i, len_tmp, len_string

    string = ''
    do i=1,size(rarray)
      write(chtmp,'(G16.9)') rarray(i)
      chtmp = adjustl(chtmp)
      len_tmp = len_trim(chtmp)
      len_string  = len_trim(string)
      string(len_string+1:len_string+len_tmp) = trim(chtmp)
      string(len_string+len_tmp+1:len_string+len_tmp+1) = ','
    end do
    len_string = len_trim(string)
    string(len_string:len_string) = ' ' ! remove trailing comma
  end function rarray_to_char

  !> basic error handler
  !!
  !! A very basic error handler. Uses `ABORT` and `FLUSH` calls
  subroutine mpp_error_basic(errortype, errormsg)
    integer, intent(in) :: errortype
    character(len=*), intent(in), optional :: errormsg

    character(len=512) :: text
    logical :: opened
    integer :: istat, errunit, outunit

    select case(errortype)
    case(NOTE)
      text = 'NOTE'         !just FYI
    case(WARNING)
      text = 'WARNING'      !probable error
    case(FATAL)
      text = 'FATAL'        !fatal error
    case default
      text = 'WARNING: non-existent errortype (must be NOTE|WARNING|FATAL)'
    end select

    if (PRESENT(errormsg)) text = trim(text)//': '//trim(errormsg)

    errunit = ERROR_UNIT
    outunit = OUTPUT_UNIT

    select case(errortype)
    case(NOTE)
      write( outunit,'(a)' )trim(text)
    case default
      write( errunit,'(/a/)' )trim(text)
      write( outunit,'(/a/)' )trim(text)
      if (errortype.EQ.FATAL) then
        call FLUSH(outunit)
        call ABORT() !automatically calls traceback on Cray systems
      end if
    end select

    return
  end subroutine mpp_error_basic

  !> overloads to mpp_error_basic, support for error_mesg routine in FMS
  subroutine mpp_error_mesg(routine, errormsg, errortype)
    character(len=*), intent(in) :: routine, errormsg
    integer, intent(in) :: errortype

    call mpp_error(errortype, trim(routine)//': '//trim(errormsg))
    return
  end subroutine mpp_error_mesg

  subroutine mpp_error_noargs()
    call mpp_error(FATAL)
  end subroutine mpp_error_noargs

  subroutine mpp_error_Is(errortype, errormsg1, value, errormsg2)
    integer, intent(in) :: errortype
    integer, intent(in) :: value
    character(len=*), intent(in) :: errormsg1
    character(len=*), intent(in), optional :: errormsg2
    call mpp_error( errortype, errormsg1, (/value/), errormsg2)
  end subroutine mpp_error_Is

  subroutine mpp_error_Rs(errortype, errormsg1, value, errormsg2)
    integer, intent(in) :: errortype
    REAL, intent(in) :: value
    character(len=*), intent(in) :: errormsg1
    character(len=*), intent(in), optional :: errormsg2
    call mpp_error( errortype, errormsg1, (/value/), errormsg2)
  end subroutine mpp_error_Rs

  subroutine mpp_error_Ia(errortype, errormsg1, array, errormsg2)
    integer, intent(in) :: errortype
    integer, dimension(:), intent(in) :: array
    character(len=*), intent(in) :: errormsg1
    character(len=*), intent(in), optional :: errormsg2
    character(len=512) :: string
    string = errormsg1//trim(array_to_char(array))
    if (present(errormsg2)) string = trim(string)//errormsg2
    call mpp_error_basic( errortype, trim(string))
  end subroutine mpp_error_Ia

  subroutine mpp_error_Ra(errortype, errormsg1, array, errormsg2)
    integer, intent(in) :: errortype
    REAL, dimension(:), intent(in) :: array
    character(len=*), intent(in) :: errormsg1
    character(len=*), intent(in), optional :: errormsg2
    character(len=512) :: string
    string = errormsg1//trim(array_to_char(array))
    if (present(errormsg2)) string = trim(string)//errormsg2
    call mpp_error_basic( errortype, trim(string))
  end subroutine mpp_error_Ra

#define _SUBNAME_ mpp_error_ia_ia
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include "mpp_error_a_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ia_ra
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include "mpp_error_a_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ra_ia
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include "mpp_error_a_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ra_ra
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include "mpp_error_a_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ia_is
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include "mpp_error_a_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ia_rs
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include "mpp_error_a_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ra_is
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include "mpp_error_a_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_ra_rs
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include "mpp_error_a_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_is_ia
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include "mpp_error_s_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_is_ra
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include "mpp_error_s_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_rs_ia
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include "mpp_error_s_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_rs_ra
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include "mpp_error_s_a.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_is_is
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ integer
#include "mpp_error_s_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_is_rs
#define _ARRAY1TYPE_ integer
#define _ARRAY2TYPE_ real
#include "mpp_error_s_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_rs_is
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ integer
#include "mpp_error_s_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_

#define _SUBNAME_ mpp_error_rs_rs
#define _ARRAY1TYPE_ real
#define _ARRAY2TYPE_ real
#include "mpp_error_s_s.h"
#undef _SUBNAME_
#undef _ARRAY1TYPE_
#undef _ARRAY2TYPE_
end module mpp_mod
