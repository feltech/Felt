dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(Felt, 1.0)
AC_PROG_CC
AC_PROG_CXX
AC_OPENMP
AX_CXX_COMPILE_STDCXX_11


AC_CANONICAL_SYSTEM

AM_CFLAGS = -g3 -O0
AM_CXXFLAGS = -g3 -O0


dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl Check if Libtool is present
dnl Libtool is used for building share libraries 
AC_PROG_LIBTOOL

AC_CONFIG_FILES(
	Makefile tests/Makefile libFelt/Makefile include/Makefile python/Makefile
)
#AC_CONFIG_FILES(libFelt/Makefile)                
AC_OUTPUT