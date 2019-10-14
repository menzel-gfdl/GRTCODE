/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <string.h>
#include "extern.h"
#include "verbosity.h"


static int verbosity = RS_NONE;
static char error_buffer[4096];

EXTERN void rs_set_verbosity(int const level)
{
    verbosity = level;
}


EXTERN int rs_get_verbosity()
{
    return verbosity;
}


EXTERN void reset_error_buffer()
{
    memset(error_buffer, '\0', 4096);
}


EXTERN void append_to_error_buffer(char const * const mesg)
{
    char b[4096];
    snprintf(b, 4096, "%s", error_buffer);
    snprintf(error_buffer, 4096, "%s%s", b, mesg);
}


EXTERN void copy_error_buffer(char * const buffer, int const buffer_size)
{
    snprintf(buffer, buffer_size, "%s\n", error_buffer);
}
