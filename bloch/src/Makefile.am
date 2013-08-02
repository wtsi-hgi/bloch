AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\"

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = bloch
bloch_SOURCES = bloch.c blog.c 
bloch_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
bloch_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
bloch_LDADD = $(top_srcdir)/gl/libbloch.la 

# include these headers in the dist but not in the install
noinst_HEADERS = bloh.h	

