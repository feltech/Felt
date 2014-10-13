/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : David Feltell
 Version     :
 Copyright   : All rights reserved
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/libFelt/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include "Felt.hpp"

int main(void) {
  print_hello();
  return 0;
}
