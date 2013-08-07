BLOCH
=====

Genotype imputation software written in C/C++ that implements the methods described 
in [Browning and Browning (2007)](http://dx.doi.org/10.1086/521987). 


Prerequisities
--------------

 * [LEMON](http://lemon.cs.elte.hu/) >= 1.2.3
 * [zlib](http://zlib.net/) >= 1.2.8
 * [htslib](https://github.com/samtools/htslib/)


Preparing source repository for compilation
-------------------------------------------

First, install the build dependencies:
 * [LEMON](http://lemon.cs.elte.hu/)
 * [htslib](https://github.com/samtools/htslib/)
 * [autoconf](https://www.gnu.org/software/autoconf/)
 * [automake](https://www.gnu.org/software/automake/)
 * [git](http://git-scm.com/)
  
Then clone the source via git and run ./autogen.sh
```bash
git clone https://github.com/wtsi-hgi/bloch.git
cd bloch/bloch
./autogen.sh
```
The autogen.sh script will bootstrap the checked-out sources to prepare them for compilation (in particular it will clone the current version of gnulib into the repository and generate the ./configure script using autoreconf. 

When it has completed successfully, you should now see the message:
```
./bootstrap: done.  Now you can run './configure'.
```


Compiling from source
---------------------
First, obtain source either following the directions above to prepare a checkout of the source repository for compilation or by downloading a distribution tarball. 

If all libraries are installed in standard system locations, you should be able to configure and compile by running:
```bash
./configure 
make
make install
```

However, if you have installed libraries into non-standard locations and/or if you wish to install BLOCH into a non-standard location, you should run ./configure with other options.

For example, to install BLOCH into a subdirectory of your home directory called 'local', you could run:
```bash
./configure --prefix=$(echo ~/local)
make
make install
```

Or, in addition to installing into ~/local, to let pkg-config know to search for prerequisite libraries (such as LEMON) from there as well, you could run:
```bash
./configure PKG_CONFIG_PATH=$(echo ~/local/lib/pkgconfig/):$PKG_CONFIG_PATH --prefix=$(echo ~/local)
make 
make install
```

