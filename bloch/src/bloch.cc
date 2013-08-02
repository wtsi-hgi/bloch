/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * BLOCH - main entry point
 *
 * Copyright (c) 2013 Genome Research Ltd. 
 * 
 * Authors: 
 *          Michelle Parker <mp18@sanger.ac.uk>
 *          Joshua C. Randall <jcrandall@alum.mit.edu>
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

#include <getopt.h>

#include <htslib/vcf.h>

#include <lemon/list_graph.h>
#include <lemon/lgf_writer.h>



int main(int argc, char** argv) 
{
  
  char *vcf_in_file;
  char *vcf_out_file;

  /* setup progname */
  set_program_name (argv[0]);

  blog_init();
  
  blog(9, "BLOCH main: started");
  
  
  /* get command-line options */
  while (1)
    {
      int c;
      int option_index = 0;
      static struct option bloch_options[] =
	{
 	  {"verbose",	        optional_argument,	0,	 0 },
 	  {"verbose",           no_argument,		0,	'v'},
	  {"debug",		no_argument,		0,	'd'},
	  {"help",		no_argument,		0,	'h'},
	  {"vcf_in",	        required_argument,	0,	'i'},
	  {"vcf_out",	        required_argument,	0,	'o'},
	  {0, 0, 0, 0}
	};
      
      c = getopt_long(argc, argv, "vdhi:o:", bloch_options, &option_index);
      
      if (c < 0)
	break;
      
      switch (c)
	{
	case 0:
	  if (!strcmp(bloch_options[option_index].name, "verbose")) {
	    if (optarg)
	      blog_set_verbosity(atoi(optarg));
	    else 
	      blog_increase_verbosity();
	  }
	  break;
	case 'v': 
	  verbosity++;
	  break;
	case 'd':
	  debug_flag = 1;
	  break;
	case 'h':
	  print_help();
	  exit(0);
	  break;
	case 'i':
	  vcf_in_file = xstrdup(optarg);
	  break;
	case 'o':
	  vcf_out_file = xstrdup(optarg);
	  break;
	case '?':
	  /* getopt_long will have already printed an error */
	  print_usage();
	  break;
	default:
	  error(0, 0, "unhandled option [-%c]", c);
	  print_usage();
	}
    }
  
  if (verbosity > 0) 
    {
      error(0, 0, "verbosity set to %u", verbosity);
    }
  
  lemon::ListGraph g;
  g.addNode();
  
  digraphWriter(g, std::cout).
    //    arcMap("capacity", cap).
    //    node("source", src).
    //    node("target", trg).
    run();
}


