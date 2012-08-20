#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([janus],[0.2])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])
AC_CONFIG_SRCDIR([src/janus.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_ARG_WITH([gsl-include],
            [AS_HELP_STRING([--with-gsl-include],[location of the gsl headers])],
            [GSL_CXXFLAGS="-I$withval"],
            [GSL_CXXFLAGS="-I/usr/include/gsl"])
AC_SUBST([GSL_CXXFLAGS])

AC_ARG_WITH([trilinos-include],
           [AS_HELP_STRING([--with-trilinos-include],[location of the trilinos headers])],
           [TRILINOS_CXXFLAGS="-I$withval"],
           [TRILINOS_CXXFLAGS="-I/usr/include/trilinos"])
AC_SUBST([TRILINOS_CXXFLAGS])

AC_ARG_WITH([lapack-include],
           [AS_HELP_STRING([--with-lapack-include],[location of the lapack headers])],
           [LAPACK_CXXFLAGS="-I$withval"],
           [LAPACK_CXXFLAGS=" "])
AC_SUBST([LAPACK_CXXFLAGS])

AC_ARG_WITH([blas-include],
           [AS_HELP_STRING([--with-blas-include],[location of the blas headers])],
           [BLAS_CXXFLAGS="-I$withval"],
           [BLAS_CXXFLAGS=" "])
AC_SUBST([BLAS_CXXFLAGS])

AC_ARG_WITH([mpi-include],
           [AS_HELP_STRING([--with-mpi-include],[location of the mpi headers])],
           [MPI_CXXFLAGS="-I$withval"],
           [MPI_CXXFLAGS="-I/usr/include/mpi"])
AC_SUBST([MPI_CXXFLAGS])
AC_ARG_WITH([gsl-lib],
            [AS_HELP_STRING([--with-gsl-lib],[location of the gsl library])],
            [GSL_LIBS="-L$withval -lgsl"],
            [GSL_LIBS="-L/usr/lib -lgsl"])
AC_SUBST([GSL_LIBS])

AC_ARG_WITH([trilinos-lib],
            [AS_HELP_STRING([--with-trilinos-lib],[location of the trilinos library])],
            [TRILINOS_LIBS="-L$withval -ltrilinos_teuchos -ltrilinos_aztecoo -ltrilinos_epetra -ltrilinos_ifpack -ltrilinos_ml"],
            [TRILINOS_LIBS="-L/usr/lib -ltrilinos_teuchos -ltrilinos_aztecoo -ltrilinos_epetra -ltrilinos_ifpack -ltrilinos_ml"])
AC_SUBST([TRILINOS_LIBS])

AC_ARG_WITH([lapack-lib],
            [AS_HELP_STRING([--with-lapack-lib],[location of the lapack library])],
            [LAPACK_LIBS="-L$withval -llapack"],
            [LAPACK_LIBS="-L/usr/lib -llapack"])
AC_SUBST([LAPACK_LIBS])

AC_ARG_WITH([blas-lib],
            [AS_HELP_STRING([--with-blas-lib],[location of the blas library])],
            [BLAS_LIBS="-L$withval -lblas"],
            [BLAS_LIBS="-L/usr/lib -lblas"])
AC_SUBST([BLAS_LIBS])

AC_ARG_WITH([mpi-lib],
            [AS_HELP_STRING([--with-mpi-lib],[location of the mpi library])],
            [MPI_LIBS="-L$withval -lmpi -lmpi_cxx"],
            [MPI_LIBS="-L/usr/lib -lmpi -lmpi_cxx"])
AC_SUBST([MPI_LIBS])


AC_ARG_WITH([agmg],
            [AS_HELP_STRING([--with-agmg],[enable agmg])],
            [with_agmg=yes],
            [with_agmg=no])

AS_IF([test "x$with_agmg" != xno],
      [AC_ARG_WITH([agmg-lib],
            [AS_HELP_STRING([--agmg-path],[location of agmg])],
            [AGMG_LIBS="$withval"],
            [AGMG_LIBS="/home/bruno/Documents/vendors/AGMG/SRC/dagmg.o /home/bruno/Documents/vendors/AGMG/SRC/dagmg_mumps.o"])],
      [AC_ARG_WITH([agmg-lib],
            [AS_HELP_STRING([--agmg-path],[location of agmg])],
            [AGMG_LIBS=" "],
            [AGMG_LIBS=" "])])

AC_SUBST([AGMG_LIBS])

AS_IF([test "x$with_agmg" != xno],
      [AC_ARG_WITH([agmg-lib],[],
            [AGMG_CXXFLAGS="-D AGMG"],
            [AGMG_CXXFLAGS="-D AGMG"])],
      [AC_ARG_WITH([agmg-lib],
            [AS_HELP_STRING([--agmg-path],[location of agmg])],
            [AGMG_CXXFLAGS=" "],
            [AGMG_CXXFLAGS=" "])])

AC_SUBST([AGMG_CXXFLAGS])
 
# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files. 

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.

AC_OUTPUT()

