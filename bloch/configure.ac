# Process this file with autoconf to produce a configure script.
AC_INIT([bloch], [0.1], [jcrandall@alum.mit.edu])
AC_PREREQ([2.69])

# Check to make sure that the src dir actually exists (sanity check)
AC_CONFIG_SRCDIR([src/bloch.cc])

# Configure aux dirs
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_MACRO_DIRS([m4 gl/m4])

# Enable extensions to C or Posix (which gnulib will provide on other systems -- this must be called before compiler is run)
AC_USE_SYSTEM_EXTENSIONS

# Initialize libtool
LT_INIT([])
LT_PREREQ([2.4.2])

# Initialize automake
AC_PROG_MAKE_SET
AM_INIT_AUTOMAKE([foreign 1.13 -Wall -Werror dist-bzip2])
#AM_GNU_GETTEXT_VERSION([0.18.2])

# Check for C++ compiler
AC_PROG_CXX

# Check for C compiler
AC_PROG_CC
gl_EARLY

# Configure automake to use modulename_CFLAGS, modulename_LDFLAGS, etc.
AM_PROG_CC_C_O

# Test compiler functionality
AC_LANG_C

# Make sure we have install program
AC_PROG_INSTALL

# Checks for stdlib.h, stdarg.h , string.h and float.h, defines STDC_HEADERS on success
AC_HEADER_STDC

# Set VERSION
# TODO: get this from git?
VERSION="0.1"
AC_SUBST(VERSION)

# Bring in config headers
AC_CONFIG_HEADERS([config.h])

# Initialise and check gnulib modules (note: gl_EARLY must also be called before this and immediately after AC_PROG_CC)
gl_INIT

# Check for pkg-config
PKG_PROG_PKG_CONFIG

# Checks for zlib and adds -lz to LIBS and defined HAVE_LIBZ
AC_ARG_VAR([ZLIB_CFLAGS],[C compiler flags for ZLIB])
AC_ARG_VAR([ZLIB_LIBS],[linker flags for ZLIB])
PKG_CHECK_MODULES([ZLIB], [zlib >= 1.2.8])
AC_MSG_CHECKING([for zlib])
AC_CHECK_LIB([z], [zlibVersion], [], [AC_MSG_FAILURE([zlib is required but check for zlibVersion function failed! (is ZLIB_LIBS set correctly?)])], [${ZLIB_LIBS}])

# Check for htslib (which requires zlib)
AC_ARG_VAR([HTSLIB_CFLAGS],[C compiler flags for HTSLIB])
AC_ARG_VAR([HTSLIB_LIBS],[linker flags for HTSLIB])
PKG_CHECK_MODULES([HTSLIB],[htslib >= 0.0.0.1])
AC_MSG_CHECKING([for htslib])
AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_FAILURE([htslib is required but check for hts_open function failed! (is HTSLIB_LIBS set correctly?)])], [${LIBS} ${HTSLIB_LIBS} ${ZLIB_LIBS} ${LTLIBMULTITHREAD}])


# Check for LEMON (no C bindings, only C++, so cannot use AC_CHECK_LIB)
AC_ARG_VAR([LEMON_CFLAGS],[C compiler flags for LEMON])
AC_ARG_VAR([LEMON_LIBS],[linker flags for LEMON])
PKG_CHECK_MODULES([LEMON],[lemon >= 1.2.3])

#AC_MSG_CHECKING([for lemon])
#AX_CXX_CHECK_LIB([emon], [lemon::ListDigraph], [], [AC_MSG_FAILURE([lemon is required but check for lemon::ListDigraph failed! (is LEMON_LIBS set correctly?)])])
#, [${LIBS} ${HTSLIB_LIBS} ${ZLIB_LIBS} ${LTLIBMULTITHREAD}])
#AC_LANG(C++)
#SAVED_LIBS=$LIBS
#LIBS="${LIBS} ${LEMON_LIBS} ${LTLIBMULTITHREAD}"
#AC_LINK_IFELSE(
#	[AC_LANG_PROGRAM([#include <lemon/list_graph.h>],
#			 [lemon::ListDigraph dummy])],
#	[],
#	[AC_MSG_FAILURE([lemon is required but check for lemon::ListDigraph failed! (is LEMON_LIBS set correctly?)])])
#LIBS=$SAVED_LIBS
#AC_LANG_POP([C++])


# Setup GetText for internationalisation
#AM_GNU_GETTEXT([external])

# Which files to configure 
AC_CONFIG_FILES([
 Makefile 
 src/Makefile 
 gl/Makefile
 test/Makefile
 test/vars
])
# po/Makefile.in


# for TAP testing
AC_REQUIRE_AUX_FILE([tap-driver.sh])
AC_PROG_AWK 
AC_ARG_VAR([HTSCMD],[absolute path to htscmd binary (from HTSLIB), used in testing])
AC_PATH_PROG([HTSCMD], [htscmd])
if test -z "$HTSCMD"
then
	AC_MSG_WARN([htscmd not found, make check will fail])
fi

AC_ARG_VAR([DIFF],[absolute path to diff binary, used in testing])
AC_PATH_PROG([DIFF], [diff])
if test -z "$DIFF"
then
	AC_MSG_WARN([diff not found, make check will fail])
fi

# Generate all config_files
AC_OUTPUT
