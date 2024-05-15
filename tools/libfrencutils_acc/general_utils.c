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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "general_utils.h"

/*******************************************************************************
  double maxval_double(int size, double *data)
  get the maximum value of double array
*******************************************************************************/
double maxval_double_acc(int size, const double *data)
{
  int n;
  double maxval;

  maxval = data[0];
  for(n=1; n<size; n++){
    if( data[n] > maxval ) maxval = data[n];
  }

  return maxval;

} /* maxval_double */


/*******************************************************************************
  double minval_double(int size, double *data)
  get the minimum value of double array
*******************************************************************************/
double minval_double_acc(int size, const double *data)
{
  int n;
  double minval;

  minval = data[0];
  for(n=1; n<size; n++){
    if( data[n] < minval ) minval = data[n];
  }

  return minval;

} /* minval_double */

/*******************************************************************************
  double avgval_double(int size, double *data)
  get the average value of double array
*******************************************************************************/
double avgval_double_acc(int size, const double *data)
{
  int n;
  double avgval;

  avgval = 0;
  for(n=0; n<size; n++) avgval += data[n];
  avgval /= size;

  return avgval;

} /* avgval_double */
