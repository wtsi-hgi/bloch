/*
 * blog.c - logging
 *
 * Copyright (c) 2013 Genome Research Ltd. 
 * Author: Joshua C. Randall <jcrandall@alum.mit.edu>
 *
 * This file is part of BLOCH. 
 *
 * BLOCH is free software: you can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 3 of the License, or (at your option) any later 
 * version. 
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT 
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  
 * 
 * You should have received a copy of the GNU General Public License along with 
 * this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "config.h"

#include <stdarg.h>
#include <stdio.h>

/* gnulib headers */
#include "progname.h"

#include "blog.h"

/*
 * blog_init
 *
 * Initialise logging system
 */
void blog_init()
{
  /* init logging globals */
  verbosity = 0;
  debug_flag = 0;
}


/*
 * blog_init
 *
 * Set verbosity level
 */
void blog_set_verbosity(unsigned int v)
{
  verbosity = v;
}


/*
 * blog_increase_verbosity
 *
 * Increase the logging verbosity by 1
 */
void blog_increase_verbosity()
{
  verbosity++;
}


/*
 * blog
 *
 * If LEVEL is less than or equal to the global verbosity level or if debug_flag 
 * is set, prefixes messages with program name and log level, writes MSGFMT to 
 * stderr, passing remaining arguments to fprintf for replacement into MSGFMT, 
 * and finishing with a newline.
 *
 */
void blog(unsigned int level, const char *msgfmt, ...)
{
  if(verbosity >= level || debug_flag) {
    va_list argp;
    fprintf(stderr, "%s(%u): ", program_name, level);
    va_start(argp, msgfmt);
    vfprintf(stderr, msgfmt, argp);
    va_end(argp);
    fprintf(stderr, "\n");
  }
}

