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
!! \brief Defines useful constants for Earth. Constants are defined as real
!! parameters. Constants are accessed through the "use" statement.
!!
!! \author Bruce Wyman <Bruce.Wyman@noaa.gov>
module constants_mod

implicit none
private

real, public, parameter :: PI = 3.14159265358979323846 !< Ratio of circle circumference to diameter [N/A]
end module constants_mod
