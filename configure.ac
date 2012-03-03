#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([janus],[0.1])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_SRCDIR([src/janus.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.
AC_ARG_WITH([gsl-include],
            [AS_HELP_STRING([--with-gsl-include],[location of the gsl headers])],
            [GSL_CXXFLAGS="-I$/withval"],
            [GSL_CXXFLAGS='-I/usr/include/gsl'])
AC_SUBST([GSL_CXXFLAGS])

AC_ARG_WITH([trilinos-include],
           [AS_HELP_STRING([--with-trilinos-include],[location of the trilinos
            hearders])],
           [TRILINOS_CXXFLAGS="-I$/withval"],
           [TRILINOS_CXXFLAGS='-I/usr/include/trilinos'])
AC_SUBST([TRILINOS_CXXFLAGS])

AC_ARG_WITH([lapack-include],
           [AS_HELP_STRING([--with-lapack-include],[location of the lapack
            hearders])],
           [LAPACK_CXXFLAGS="-I$/withval"],
           [LAPACK_CXXFLAGS=' '])
AC_SUBST([LAPACK_CXXFLAGS])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_ARG_WITH([gsl-lib],
            [AS_HELP_STRING([--with-gsl-lib],[location of the gsl library])],
            [GSL_LIBS="-L$withval -lgsl"],
            [GSL_LIBS="-L/usr/lib -lgsl"])
AC_SUBST([GSL_LIBS])

AC_ARG_WITH([trilinos-lib],
            [AS_HELP_STRING([--with-trilinos-lib],[location of the trilinos 
             library])],
            [TRILINOS_LIBS="-L$withval -ltrilinos_teuchos"],
            [TRILINOS_LIBS="-L/usr/lib -ltrilinos_teuchos"])
AC_SUBST([TRILINOS_LIBS])

AC_ARG_WITH([lapack-lib],
            [AS_HELP_STRING([--with-lapack-lib],[location of the lapack 
             library])],
            [LAPACK_LIBS="-L$withval -llapack"],
            [LAPACK_LIBS="-L/usr/lib -llapack"])
AC_SUBST([LAPACK_LIBS])

AC_OUTPUT()

