code_dir = $(abs_top_srcdir)/src
VPATH = $(code_dir)
OBJECTS = cell_test cross_sections_fp_test_1 cross_sections_fp_test_2\
				 	cross_sections_regular_test_1 cross_sections_regular_test_2\
				 	cross_sections_cepxs_test dof_handler_test edge_test glc_test\
					ls_test mip_test parameters_test_1 parameters_test_2 pwld_test\
				 	transport_operator_test_1	transport_operator_test_2\
					triangulation_rec_test triangulation_pol_test

CXXFLAGS = -Wall $(AGMG_CXXFLAGS) $(GSL_CXXFLAGS) $(TRILINOS_CXXFLAGS) $(LAPACK_CXXFLAGS) $(BLAS_CXXFLAGS) $(MPI_CXXFLAGS)
LIBFLAGS = $(AGMG_LIBS) $(GSL_LIBS) $(TRILINOS_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(MPI_LIBS)

cell_test : cell_test.o BLD.o CELL.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cell_test.o : cell_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

cross_sections_fp_test_1 : cross_sections_fp_test_1.o CROSS_SECTIONS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cross_sections_fp_test_1.o : cross_sections_fp_test_1.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

cross_sections_fp_test_2 : cross_sections_fp_test_2.o CROSS_SECTIONS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cross_sections_fp_test_2.o : cross_sections_fp_test_2.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

cross_sections_regular_test_1 : cross_sections_regular_test_1.o CROSS_SECTIONS.o\
 	EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cross_sections_regular_test_1.o : cross_sections_regular_test_1.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

cross_sections_regular_test_2 : cross_sections_regular_test_2.o CROSS_SECTIONS.o\
 	EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cross_sections_regular_test_2.o : cross_sections_regular_test_2.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

cross_sections_cepxs_test : cross_sections_cepxs_test.o CROSS_SECTIONS.o\
 	EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

cross_sections_cepxs_test.o : cross_sections_cepxs_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

dof_handler_test : dof_handler_test.o CELL.o DOF_HANDLER.o EDGE.o PARAMETERS.o\
	QUADRATURE.o TRIANGULATION.o GLC.o BLD.o PWLD.o EXCEPTION.o CROSS_SECTIONS.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

dof_handler_test.o : dof_handler_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

edge_test : edge_test.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

edge_test.o : edge_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

glc_test : glc_test.o GLC.o QUADRATURE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

glc_test.o : glc_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $< 

ls_test : ls_test.o LS.o QUADRATURE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

ls_test.o : ls_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

mip_test : mip_test.o DOF_HANDLER.o GLC.o MIP.o PARAMETERS.o QUADRATURE.o\
	TRIANGULATION.o MIP.o BLD.o EDGE.o PWLD.o CELL.o EXCEPTION.o CROSS_SECTIONS.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

mip_test.o : mip_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

parameters_test_1 : parameters_test_1.o PARAMETERS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

parameters_test_1.o : parameters_test_1.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

parameters_test_2 : parameters_test_2.o PARAMETERS.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

parameters_test_2.o : parameters_test_2.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

pwld_test : pwld_test.o PWLD.o CELL.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

pwld_test.o : pwld_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

triangulation_pol_test : triangulation_pol_test.o TRIANGULATION.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

triangulation_pol_test.o : triangulation_pol_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

triangulation_rec_test : triangulation_rec_test.o TRIANGULATION.o EDGE.o EXCEPTION.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

triangulation_rec_test.o : triangulation_rec_test.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

transport_operator_test_1 : transport_operator_test_1.o CELL.o EDGE.o LS.o MIP.o\
	PARAMETERS.o QUADRATURE.o TRANSPORT_OPERATOR.o TRIANGULATION.o DOF_HANDLER.o BLD.o\
	PWLD.o EXCEPTION.o CROSS_SECTIONS.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

transport_operator_test_1.o : transport_operator_test_1.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<

transport_operator_test_2 : transport_operator_test_2.o CELL.o EDGE.o LS.o MIP.o\
	PARAMETERS.o QUADRATURE.o TRANSPORT_OPERATOR.o TRIANGULATION.o DOF_HANDLER.o BLD.o\
	PWLD.o EXCEPTION.o CROSS_SECTIONS.o GLC.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBFLAGS)

transport_operator_test_2.o : transport_operator_test_2.cc
	$(CXX) $(CXXFLAGS) -I/$(code_dir) -c $<


EXTRA_DIST = testsuite.at $(TESTSUITE)
TESTSUITE = $(srcdir)/testsuite
		      
check-local: atconfig $(TESTSUITE) $(OBJECTS)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(abs_builddir)' $(TESTSUITEFLAGS)
														     
installcheck-local: atconfig $(TESTSUITE) $(OBJECTS)
	$(SHELL) '$(TESTSUITE)' AUTOTEST_PATH='$(bindir)' $(TESTSUITEFLAGS)
																									      
clean-local:
	test ! -f '$(TESTSUITE)' || $(SHELL) '$(TESTSUITE)' --clean
	-rm testsuite
	-rm *.o
	-rm cell_test
	-rm cross_sections_fp_test_1
	-rm cross_sections_fp_test_2
	-rm cross_sections_regular_test_1
	-rm cross_sections_regular_test_2
	-rm cross_sections_cepxs_test
	-rm dof_handler_test
	-rm edge_test
	-rm glc_test
	-rm ls_test
	-rm mip_test
	-rm parameters_test_1
	-rm parameters_test_2
	-rm pwld_test
	-rm transport_operator_test_1
	-rm transport_operator_test_2
	-rm triangulation_pol_test
	-rm triangulation_rec_test

$(TESTSUITE) : $(srcdir)/testsuite.at 
	autom4te -l autotest '$(srcdir)/testsuite.at' -o testsuite
