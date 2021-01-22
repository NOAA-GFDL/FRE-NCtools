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
!-----------------------------------------------------------------------
! Copyright 2011 NOAA Geophysical Fluid Dynamics Lab, Princeton, NJ
! This program is distributed under the terms of the GNU General Public
! License.
!-----------------------------------------------------------------------

module plev_constants_mod
implicit none
private

real, public, parameter :: GRAV   = 9.80
real, public, parameter :: RDGAS  = 287.04
real, public, parameter :: RVGAS  = 461.50
real, public, parameter :: TFREEZE  = 273.16

end module plev_constants_mod
