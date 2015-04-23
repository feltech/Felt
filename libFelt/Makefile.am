#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
noinst_LTLIBRARIES = libFelt.la

#######################################
# Build information for each library

# Sources for libFelt
libFelt_la_SOURCES = gui/RktFileInterface.cpp gui/RktSFMLInterface.cpp \
	gui/RktSFMLRenderer.cpp gui/RktSFMLSystem.cpp

# Linker options libTestProgram
libFelt_la_LDFLAGS =

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libFelt_la_CPPFLAGS = -I$(top_srcdir)/include \
	-I$(top_srcdir)/lib/rocket/include $(OPENMP_CFLAGS)