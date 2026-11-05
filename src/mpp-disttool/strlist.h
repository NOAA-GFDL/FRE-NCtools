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
#ifndef MPPDISTTOOL_STRLIST_H
#define MPPDISTTOOL_STRLIST_H 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "xmalloc.h"

int newstringlist(char ***list, int *n, int size);
int clearstringlist(char **list, int n);
void getstringlist(char *optarg, char ***list, int *nitems);
void freestringlist(char ***list, int nitems);
int instringlist(char **list, char *str, int nitems);
int addstringtolist(char **list, char *string, int nitems);
int appendstringtolist(char ***list, char *string, int *nitems);
void printstrlist(char **list, int n, FILE *f);
int strlistu(char **list1, char **list2, char **listunion, int n1, int n2, int nu);
int strlistsd(char **list1, char **list2, char **listdiff, int n1, int n2, int nsd);
int copystrlist(char **listsrc, char **listdst, int nsrc, int ndst);
int getnumstrlist(char **list, int nlist);

#endif /* MPPDISTTOOL_STRLIST_H */
