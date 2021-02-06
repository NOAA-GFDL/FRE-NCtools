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

program nc_null

use netcdf
use iso_c_binding

implicit none

character (len=200) :: filename !< The name of the file to be checked
character (len=10) :: varname !< The variable to be checked
character (len=6) :: attribute_name="bounds" !< The attirbute to be checked
character (len=20) :: attribute_value !< THe value of the attribute
integer :: ncid !< The netcdf file ID
integer :: varid !< The netcdf variable ID
integer :: att_len !< The length of the attribute string
integer :: i !< For looping
logical :: null_found = .false. !< True if a null character is found
integer :: arglen, istatus

!> Get the file name from the command line
call get_command_argument(1 , filename, arglen, istatus)
if (istatus .ne. 0) then
     !> Print error message if no file is given
     write (6,*) "Please enter a FILENAME argument on the command line"
     write (6,*) "nc-null-check FILENAME"
     stop
endif
!> Tell the user the name of the file being processed
write (6,*) "Checking "//trim(filename)
write (6,*) " "
!> Check the lat variable
varname = "lat"
!> Open file
  call check( nf90_open(FILENAME, nf90_nowrite, ncid) )
!> The the varid
  call check( nf90_inq_varid(ncid, varname, varid) )
!> Get the value of the desired attribute
  call check( nf90_get_att(ncid, varid, attribute_name, attribute_value) )
!> Get the length of the attribute value
  call check( nf90_inquire_attribute(ncid, varid, attribute_name, len = att_len) )

!> Loop through each character in the attribute value string to check for a null character
do i = 1,len(attribute_value)
     if (attribute_value(i:i) == c_null_char) then
          !> If a null character is found, inform the user
          write (6,*) "The variable "//trim(varname)//" in "//trim(filename)//" has a null character "//&
          "in attribute "//trim(attribute_name)//"::"//trim(attribute_value)
          null_found = .true.
     endif
enddo
!> If not null character is found, inform the user
if (.not. null_found) write (6,*) "No null character found"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
contains
!> Check for netcdf errors
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check

end program nc_null

