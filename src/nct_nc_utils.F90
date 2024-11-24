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

#include "config.h"

module nct_nc_utils

   use iso_fortran_env, only: error_unit
   use netcdf

   implicit none
contains
   subroutine nct_hst_att_cat (out_id, hist_string)
      integer, intent(in) :: out_id
      character(len=*), intent(in) :: hist_string
      integer :: ncstatus, numatts, attnum
      integer :: atttype, attlen
      character(len=*), parameter :: history_att_name = 'history'
      character(len=NF90_MAX_NAME) :: attname
      character(len=:), allocatable :: new_history_string, old_history_string

      ncstatus = NF90_INQUIRE(out_id, nattributes=numatts)
      do attnum=1, numatts
         ncstatus = nf90_inq_attname(out_id, NF90_GLOBAL, attnum, attname)
         if (trim(attname) == trim(history_att_name)) then
            exit
         end if
      end do
      if (attnum > numatts) then
         new_history_string = date_and_time_str()//": "//trim(hist_string)
      else
         ncstatus = NF90_INQUIRE_ATTRIBUTE(out_id, NF90_GLOBAL, history_att_name, len=attlen, xtype=atttype)
         if (atttype /= NF90_CHAR) then
            write (error_unit,*) 'History attribute is not a character type'
            return
         end if
         allocate(character(len=attlen) :: old_history_string)
         ncstatus = NF90_GET_ATT(out_id, NF90_GLOBAL, history_att_name, old_history_string)
         new_history_string = date_and_time_str()//": "//trim(hist_string)//new_line('A')//old_history_string
      end if
      ! Add/prepend history attribute
      ncstatus = NF90_PUT_ATT(out_id, NF90_GLOBAL, history_att_name, new_history_string)
      ! Add FRE NCtools version information to the file
      ncstatus = NF90_PUT_ATT(out_id, NF90_GLOBAL, PACKAGE,&
           & PACKAGE_NAME//" version "//PACKAGE_VERSION//" (git: "//GIT_REVISION//")")
   end subroutine nct_hst_att_cat

   character(len=24) function date_and_time_str()
      character(len=3), dimension(12), parameter :: month_name =&
           & ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
      integer, dimension(8) :: time

      call date_and_time(VALUES=time)
      write (date_and_time_str, FMT='(A3," ",A3," ",I0.2," ",I0.2,":",I0.2,":",I0.2," ",I0.4)') &
           & day_of_week_str(time(1), time(2), time(3)), &
           & trim(adjustl(month_name(time(2)))), &
           & time(3), &
           & time(5), time(6), time(7), time(1)
   end function date_and_time_str

   character(len=3) function day_of_week_str(year, month, day)
      integer, intent(in) :: year, month, day
      integer :: day_of_week
      integer, dimension(12) :: month_key = [1, 4, 4, 0, 2, 5, 0, 3, 6, 1, 4, 6]
      integer, dimension(12) :: month_key_leap = [0, 3, 4, 0, 2, 5, 0, 3, 6, 1, 4, 6]
      character(len=3), dimension(7) :: day_of_week_string = ["Sat", "Sun", "Mon", "Tue", "Wed", "Thu", "Fri"]

      ! Add to the year without the centery, 1/4 of the year without
      ! the century and the day of the month
      day_of_week = mod(year, 100) + mod(year, 100)/4 + day

      ! add the month key to the day of the week
      if (is_leap_year(year)) then
         day_of_week = day_of_week + month_key_leap(month)
      else
         day_of_week = day_of_week + month_key(month)
      end if

      ! Adjust based on the century
      if (year < 1800) then
         day_of_week = day_of_week + 4
      else if (year < 1900) then
         day_of_week = day_of_week + 2
      else if (year >= 2000 .and. year <= 2099) then
         day_of_week = day_of_week - 1
      end if

      ! divide by 7 and keep the remainder
      day_of_week = mod(day_of_week, 7)
      day_of_week_str = day_of_week_string(day_of_week+1)
   end function day_of_week_str

   logical function is_leap_year(year)
      integer, intent(in) :: year
      if (mod(year, 4) == 0) then
         if (mod(year, 100) /= 0 .or. mod(year, 400) == 0) then
            is_leap_year = .true.
         else
            is_leap_year = .false.
         end if
      else
         is_leap_year = .false.
      end if
   end function is_leap_year

   function itoa(i)
      integer, intent(in) :: i
      character(len=:), allocatable :: itoa
      character(len=256) :: itoa_tmp
      write (itoa_tmp, '(I0)') i
      itoa = trim(adjustl(itoa_tmp))
   end function itoa

   function rtoa(r)
      real, intent(in) :: r
      character(len=:), allocatable :: rtoa
      character(len=256) :: rtoa_tmp
      write (rtoa_tmp, '(F0.0)') r
      rtoa = trim(adjustl(rtoa_tmp))
   end function rtoa

   character(len=1) function ltoa(l)
      logical, intent(in) :: l
      if (l) then
         ltoa = 'T'
      else
         ltoa = 'F'
      end if
   end function ltoa
end module nct_nc_utils
