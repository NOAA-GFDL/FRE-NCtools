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
#include "xmalloc.h"

void *
xmalloc(size_t num)
{
    void *p = malloc(num);
    if (!p) {
        fprintf(stderr, "Memory exhausted\n");
        exit(EXIT_FATAL);
    }
    return p;
}

void *
xrealloc(void *p, size_t num)
{
    void *new;
    if (!p)
        return xmalloc(num);
    new = realloc(p, num);
    if (!new) {
        fprintf(stderr, "Memory exhausted\n");
        exit(EXIT_FATAL);
    }
    return new;
}

void *
xcalloc(size_t num, size_t size)
{
    void *p = xmalloc(num * size);
    memset(p, 0, num * size);
    return p;
}
