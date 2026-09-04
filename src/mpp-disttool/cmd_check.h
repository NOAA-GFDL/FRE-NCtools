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
#ifndef MPPDISTTOOL_CMD_CHECK_H
#define MPPDISTTOOL_CMD_CHECK_H

#include <stdbool.h>

/* Returns true (exit 0) if any dimension variable in the file has a
   "compress" attribute.  Matches is-compressed.c logic exactly. */
bool check_is_compressed(const char *filename);

/* Returns true if the file is a valid iceberg restart (has
   file_format_major_version >= 0 and file_format_minor_version >= 1,
   and has NumFilesInSet > 0). */
bool check_is_iceberg(const char *filename);

/* Returns true if the iceberg file has icebergs (UNLIMITED dim > 0). */
bool check_iceberg_has_data(const char *filename);

/* Entry point for "mppdisttool check [--compressed|--iceberg] [--quiet] FILE"
   Returns 0 (compressed/iceberg), -1 (not compressed), 255 (not iceberg). */
int cmd_check(int argc, char **argv);

#endif /* MPPDISTTOOL_CMD_CHECK_H */
