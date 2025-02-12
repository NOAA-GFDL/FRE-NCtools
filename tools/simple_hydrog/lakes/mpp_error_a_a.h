/***********************************************************************
 *                   GNU Lesser General Public License
 *
 * This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
 *
 * FRE-NCtools is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FRE-NCTools.  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/
subroutine _SUBNAME_(errortype, errormsg1, array1, errormsg2, array2, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_, dimension(:), intent(in) :: array1
  _ARRAY2TYPE_, dimension(:), intent(in) :: array2
  character(len=*),      intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3
  character(len=512) :: string

  string = errormsg1//trim(array_to_char(array1)) 
  string = trim(string)//errormsg2//trim(array_to_char(array2)) 
  if(present(errormsg3)) string = trim(string)//errormsg3
  call mpp_error_basic( errortype, trim(string))

end subroutine _SUBNAME_
