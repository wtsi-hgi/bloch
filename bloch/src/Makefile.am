AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\" $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS) $(LEMON_CFLAGS)

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = bloch
bloch_SOURCES = bloch.cc blog.c 
bloch_LDADD = $(top_srcdir)/gl/libbloch.la $(ZLIB_LIBS) $(HTSLIB_LIBS) $(LEMON_LIBS) 
bloch_LDFLAGS = -static

# include these headers in the dist but not in the install
noinst_HEADERS = blog.h	

