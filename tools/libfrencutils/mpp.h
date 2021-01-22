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
 * License along with FRE-NCTools (LICENSE.md).  If not, see
 * <http://www.gnu.org/licenses/>.
 **********************************************************************/

/*********************************************************************
                    mpp.h
  This header contains subroutine for parallel programming.
  only MPI parallel is implemented. 
  Contact: Zhi.Liang@noaa.gov
 ********************************************************************/
#ifndef MPP_H_
#define MPP_H_

void mpp_init(int *argc, char ***argv);          /* start parallel programming, create communicator */
void mpp_end();           /* end of parallel programming, abort the program */
int mpp_pe();      /* return processor ID */
int mpp_root_pe(); /* return root pe of current pelist */
int mpp_npes();    /* return number of processor used */
int* mpp_get_pelist();
void mpp_send_double(const double* data, int size, int to_pe); /* send data */
void mpp_send_int(const int* data, int size, int to_pe); /* send data */
void mpp_recv_double(double* data, int size, int from_pe); /* recv data */
void mpp_recv_int(int* data, int size, int from_pe); /* recv data */
void mpp_error(char *str);
void mpp_sum_int(int count, int *data);
void mpp_sum_double(int count, double *data);
void mpp_min_double(int count, double *data);
void mpp_max_double(int count, double *data);
void print_mem_usage(const char* text);
void print_time(const char* text, double t);
void mpp_sync_self();
void mpp_sync();
#endif
