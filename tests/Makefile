VPATH = ../src
LAPACK = -llapack
TRILINOS = -ltrilinos_aztecoo -ltrilinos_epetra -ltrilinos_teuchos -ltrilinos_ml -ltrilinos_ifpack
CXXFLAGS = -Wall -g -I/usr/include/trilinos -I/usr/include/gsl -I/usr/include/mpi -I../src
AGMG = /home/bruno/Documents/vendors/AGMG/SRC/dagmg.o /home/bruno/Documents/vendors/AGMG/SRC/dagmg_mumps.o
LIBFLAGS = $(AGMG) -L/usr/lib $(TRILINOS) $(LAPACK) -lgsl -lmpi -lmpi_cxx -lblas
OBJECTS = cell_test dof_handler_test edge_test glc_test ls_test mip_test\
				 	parameters_test_1	parameters_test_2 pwld_test	triangulation_pol_test\
				 	transport_operator_test_1 transport_operator_test_2 triangulation_rec_test 

all : $(OBJECTS)

cell_test : cell_test.o BLD.o CELL.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

dof_handler_test : dof_handler_test.o CELL.o DOF_HANDLER.o EDGE.o PARAMETERS.o\
QUADRATURE.o TRIANGULATION.o GLC.o BLD.o PWLD.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

edge_test : edge_test.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

glc_test : glc_test.o GLC.o QUADRATURE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

ls_test : ls_test.o LS.o QUADRATURE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

mip_test : mip_test.o DOF_HANDLER.o GLC.o MIP.o PARAMETERS.o QUADRATURE.o\
TRIANGULATION.o MIP.o BLD.o EDGE.o PWLD.o CELL.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

parameters_test_1 : parameters_test_1.o PARAMETERS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

parameters_test_2 : parameters_test_2.o PARAMETERS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

pwld_test : pwld_test.o PWLD.o CELL.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

triangulation_pol_test : triangulation_pol_test.o TRIANGULATION.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

triangulation_rec_test : triangulation_rec_test.o TRIANGULATION.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

transport_operator_test_1 : transport_operator_test_1.o CELL.o EDGE.o LS.o MIP.o\
PARAMETERS.o QUADRATURE.o TRANSPORT_OPERATOR.o TRIANGULATION.o DOF_HANDLER.o BLD.o\
PWLD.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

transport_operator_test_2 : transport_operator_test_2.o CELL.o EDGE.o GLC.o LS.o MIP.o\
PARAMETERS.o QUADRATURE.o TRANSPORT_OPERATOR.o TRIANGULATION.o DOF_HANDLER.o BLD.o\
PWLD.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

.PHONY : clean check all

check : testsuite
	./testsuite

testsuite : testsuite.at $(OBJECTS)
	autom4te -l autotest testsuite.at -o testsuite

clean :
	-rm *.o
	-rm *.log
	-rm $(OBJECTS)
	-rm testsuite
	-rm -r *.dir
