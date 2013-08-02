/*
 * blog.h - logging functions
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
#ifndef BLOCH_FILES_H
#define BLOCH_FILES_H

/* verbosity level (0-3; increasing with each -v option): 0 is silent, 3 is maximum verbosity */
unsigned int verbosity;

/* debug flag: if non-zero, print debugging messages to stderr */
unsigned int debug_flag;

void blog_init();
void blog_set_verbosity(unsigned int v);
void blog_increase_verbosity();

void blog(unsigned int level, const char *msgfmt, ...);
void debug(const char *msgfmt, ...);

#endif
