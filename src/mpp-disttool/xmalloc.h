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
#ifndef MPPDISTTOOL_XMALLOC_H
#define MPPDISTTOOL_XMALLOC_H 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define EXIT_FATAL  2
#define EXIT_FAILED 3

#define XCALLOC(type, num) \
    ((type*) xcalloc((size_t)(num), (size_t)sizeof(type)))

#define XMALLOC(type, num) \
    ((type*) xmalloc((size_t)((num) * sizeof(type))))

#define XREALLOC(type, p, num) \
    ((type*) xrealloc((p), (size_t)((num) * sizeof(type))))

#define XFREE(stale) \
    do { if (stale) { free(stale); stale = NULL; } } while (0)

extern void *xcalloc(size_t num, size_t size);
extern void *xmalloc(size_t num);
extern void *xrealloc(void *p, size_t num);

#endif /* MPPDISTTOOL_XMALLOC_H */
