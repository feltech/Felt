#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For example /usr/bin
#bin_PROGRAMS=exampleProgram

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
noinst_PROGRAMS=test

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for example.

# Sources for the a.out 
test_SOURCES = test_Grid.cpp test_ArrayMappedGrid.cpp test_Surface.cpp test_Poly.cpp test_all.cpp 

# Linker options for a.out
test_LDFLAGS = $(top_srcdir)/libFelt/libFelt.la

test_LDADD = -lboost_unit_test_framework -lgomp

# Compiler options for a.out
test_CPPFLAGS = -I$(top_srcdir)/include $(OPENMP_CFLAGS)
