#include <cassert>
#include <cmath>
#include <map>
#include <string>
#include <vector>
#include "gsl_math.h"
#include "mpi.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "CROSS_SECTIONS.hh"
#include "DOF_HANDLER.hh"
#include "ERROR_ESTIMATOR.hh"
#include "GLC.hh"
#include "MIP.hh"
#include "PARAMETERS.hh"
#include "QUADRATURE.hh"
#include "TRIANGULATION.hh"

#include "CELL.hh"
#include "EDGE.hh"

using namespace std;

double Compute_convergence(Epetra_MultiVector const &flux, 
    Epetra_MultiVector const &old_flux, const unsigned int n)
{
  double conv(0.);
  double* l2_norm_num = new double[n];
  double* l2_norm_denom = new double[n];
  Epetra_MultiVector tmp(flux);
  tmp.Update(-1.,old_flux,0.);
  // Compute the L2 norm of each vector of the MultiVector
  tmp.Norm2(l2_norm_num);
  flux.Norm2(l2_norm_denom);
  double norm_num(0.);
  double norm_denom(0.);
  for (unsigned int j=0; j<n; ++j)
  {
    norm_num += pow(l2_norm_num[j],2);
    norm_denom += pow(l2_norm_denom[j],2);
  }
  conv = sqrt(norm_num/norm_denom);
  delete [] l2_norm_num;
  delete [] l2_norm_denom;      

  return conv;
}
        
int main(int argc,char** argv)
{
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  d_vector solution(192);
  solution[0] = 3.65207;
  solution[1] = 1.63393;  
  solution[2] = 0.946519; 
  solution[3] = 2.29929;  
  solution[4] = 0.946519; 
  solution[5] = 0.25911;  
  solution[6] = 0.946519; 
  solution[7] = 2.29929;  
  solution[8] = 0.946519; 
  solution[9] = 1.63393;  
  solution[10] = 3.65207;   
  solution[11] = 2.29929;   
  solution[12] = 3.65207;   
  solution[13] = 5.67021;   
  solution[14] = 3.65207;   
  solution[15] = 2.29929;   
  solution[16] = 1.25823;   
  solution[17] = 2.19684;   
  solution[18] = 7.66854;   
  solution[19] = 5.94698;   
  solution[20] = 3.60261;   
  solution[21] = 7.76686;   
  solution[22] = 10.5232;   
  solution[23] = 7.76686;   
  solution[24] = 5.78508;   
  solution[25] = 7.66854;   
  solution[26] = 2.19684;   
  solution[27] = 1.25823;   
  solution[28] = 3.60261;   
  solution[29] = 5.94698;   
  solution[30] = 7.67654;   
  solution[31] = 3.08287;   
  solution[32] = 3.08287;   
  solution[33] = 7.21063;   
  solution[34] = 3.08287;   
  solution[35] = 3.08287;   
  solution[36] = 7.67654;   
  solution[37] = 7.21063;   
  solution[38] = 7.67654;   
  solution[39] = 10.8725;   
  solution[40] = 10.8725;   
  solution[41] = 7.21063;   
  solution[42] = 10.8725;   
  solution[43] = 10.8725;   
  solution[44] = 7.67654;   
  solution[45] = 7.21063;   
  solution[46] = 7.66854;   
  solution[47] = 2.19684;   
  solution[48] = 1.25823;   
  solution[49] = 3.60261;   
  solution[50] = 5.94698;   
  solution[51] = 3.65207;   
  solution[52] = 1.63393;   
  solution[53] = 0.946519;  
  solution[54] = 2.29929;   
  solution[55] = 0.946519;  
  solution[56] = 0.25911;   
  solution[57] = 0.946519;  
  solution[58] = 2.29929;  
  solution[59] = 0.946519;  
  solution[60] = 1.63393;  
  solution[61] = 3.65207;   
  solution[62] = 2.29929;  
  solution[63] = 3.65207;   
  solution[64] = 5.67021;  
  solution[65] = 3.65207;   
  solution[66] = 2.29929;   
  solution[67] = 1.25823;   
  solution[68] = 2.19684;   
  solution[69] = 7.66854;   
  solution[70] = 5.94698;   
  solution[71] = 3.60261;   
  solution[72] = 7.76686;   
  solution[73] = 10.5232;   
  solution[74] = 7.76686;   
  solution[75] = 5.78508;   
  solution[76] = 3.08287;   
  solution[77] = 3.08287;   
  solution[78] = 7.67654;   
  solution[79] = 7.21063;   
  solution[80] = 7.67654;   
  solution[81] = 10.8725;   
  solution[82] = 10.8725;   
  solution[83] = 7.21063;   
  solution[84] = 10.8725;   
  solution[85] = 10.8725;   
  solution[86] = 7.67654;   
  solution[87] = 7.21063;   
  solution[88] = 7.67654;   
  solution[89] = 3.08287;   
  solution[90] = 3.08287;   
  solution[91] = 7.21063;   
  solution[92] = 10.8178;   
  solution[93] = 10.8178;   
  solution[94] = 10.8178;   
  solution[95] = 10.8178;   
  solution[96] = 10.8178;   
  solution[97] = 10.8178;   
  solution[98] = 10.8178;   
  solution[99] = 10.8178;   
  solution[100] = 10.8725;   
  solution[101] = 10.8725;   
  solution[102] = 7.67654;   
  solution[103] = 7.21063;  
  solution[104] = 7.67654;  
  solution[105] = 3.08287;  
  solution[106] = 3.08287;  
  solution[107] = 7.21063;  
  solution[108] = 3.08287;  
  solution[109] = 3.08287;  
  solution[110] = 7.67654;  
  solution[111] = 7.21063;  
  solution[112] = 7.67654;  
  solution[113] = 10.8725;  
  solution[114] = 10.8725;  
  solution[115] = 7.21063;  
  solution[116] = 1.25823;  
  solution[117] = 2.19684;  
  solution[118] = 7.66854;  
  solution[119] = 5.94698;  
  solution[120] = 3.60261;  
  solution[121] = 7.76686;  
  solution[122] = 10.5232;  
  solution[123] = 7.76686;  
  solution[124] = 5.78508;  
  solution[125] = 7.66854;  
  solution[126] = 2.19684;  
  solution[127] = 1.25823;  
  solution[128] = 3.60261;  
  solution[129] = 5.94698;  
  solution[130] = 3.65207;  
  solution[131] = 1.63393;  
  solution[132] = 0.946519; 
  solution[133] = 2.29929;  
  solution[134] = 0.946519; 
  solution[135] = 0.25911;  
  solution[136] = 0.946519; 
  solution[137] = 2.29929;  
  solution[138] = 0.946519; 
  solution[139] = 1.63393;  
  solution[140] = 3.65207;  
  solution[141] = 2.29929;  
  solution[142] = 3.65207;  
  solution[143] = 5.67021;  
  solution[144] = 3.65207;  
  solution[145] = 2.29929;  
  solution[146] = 7.67654;  
  solution[147] = 10.8725;  
  solution[148] = 10.8725;  
  solution[149] = 7.21063;  
  solution[150] = 10.8725;  
  solution[151] = 10.8725;  
  solution[152] = 7.67654;  
  solution[153] = 7.21063;  
  solution[154] = 7.67654;  
  solution[155] = 3.08287;  
  solution[156] = 3.08287;  
  solution[157] = 7.21063;  
  solution[158] = 3.08287;  
  solution[159] = 3.08287;  
  solution[160] = 7.67654;  
  solution[161] = 7.21063;  
  solution[162] = 7.76686;  
  solution[163] = 10.5232;  
  solution[164] = 7.76686;  
  solution[165] = 5.78508;  
  solution[166] = 7.66854;  
  solution[167] = 2.19684;  
  solution[168] = 1.25823;  
  solution[169] = 3.60261;  
  solution[170] = 5.94698;  
  solution[171] = 3.65207;  
  solution[172] = 1.63393; 
  solution[173] = 0.946519; 
  solution[174] = 2.29929;  
  solution[175] = 0.946519; 
  solution[176] = 0.25911;  
  solution[177] = 0.946519; 
  solution[178] = 2.29929;  
  solution[179] = 0.946519; 
  solution[180] = 1.63393;  
  solution[181] = 3.65207;  
  solution[182] = 2.29929;  
  solution[183] = 3.65207;  
  solution[184] = 5.67021;  
  solution[185] = 3.65207;  
  solution[186] = 2.29929;  
  solution[187] = 1.25823;  
  solution[188] = 2.19684;  
  solution[189] = 7.66854;  
  solution[190] = 5.94698;  
  solution[191] = 3.60261;  

  string cross_sections_inp("cross_sections_mip_diffusion.inp");
  string geometry_inp("geometry_mip_diffusion.inp");
  string parameters_inp("parameters_mip_diffusion.inp");

  CROSS_SECTIONS cross_sections(&cross_sections_inp);
  TRIANGULATION triangulation(&geometry_inp);
  PARAMETERS parameters(&parameters_inp);
   
  triangulation.Read_geometry();
  triangulation.Build_edges();

  parameters.Read_diffusion_parameters(triangulation.Get_n_sources());

  cross_sections.Read_regular_cross_sections(triangulation.Get_n_materials(),
      parameters.Get_permutation_type(),false);
  cross_sections.Apply_ang_lvls_and_tc(parameters.Get_multigrid(),
      parameters.Get_transport_correction(),parameters.Get_optimal_tc(),
      triangulation.Get_n_materials(),parameters.Get_sn_order());
  vector<d_vector> a=cross_sections.Get_sigma_t(0);

  DOF_HANDLER *dof_handler = 
    new DOF_HANDLER(&triangulation,parameters,cross_sections);

  Epetra_Map* scalar_flux_map = new Epetra_Map(dof_handler->Get_n_dof(),0,comm);
  Epetra_MultiVector flux_moments(*scalar_flux_map,1);
  Epetra_MultiVector rhs(*scalar_flux_map,1);

  unsigned int group_iter(0);
  const unsigned int max_group_it(parameters.Get_max_group_it());
  const unsigned int n_groups(cross_sections.Get_n_groups());
  const unsigned int n_refinements(parameters.Get_n_refinements());
  const double group_tol(parameters.Get_group_tolerance());
  Epetra_MultiVector* group_flux = new Epetra_MultiVector(*scalar_flux_map,
      cross_sections.Get_n_groups());
  for (unsigned int r=0; r<=n_refinements; ++r)
  {
    double group_conv(10.*group_tol);
    Epetra_MultiVector old_group_flux(*group_flux);
    Epetra_MultiVector scalar_flux(*scalar_flux_map,1);

    MIP mip(&comm,&parameters,dof_handler);
    // Loop over the groups
    while (group_conv>group_tol)
    {
      // Loop over the supergroups.
      for (unsigned int g=0; g<n_groups; ++g)
      {
        // Set the current group
        mip.Set_group(g);
        mip.Solve_diffusion(n_groups,scalar_flux,*group_flux);

        old_group_flux[g] = (*group_flux)[g];
        (*group_flux)[g] = scalar_flux[0];
      }
      // Compute the convergence over all the groups
      group_conv = Compute_convergence(*group_flux,old_group_flux,n_groups);
      ++group_iter;
      if (group_iter>=max_group_it)
        break;
    }
    // Refine the mesh
    ui_set cells_to_refine;
    ui_set adjacent_cells;
    vector<vector<unsigned int> > projection;
    map<unsigned int,vector<vector<vector<double> > > > edge_to_refine;

    // Compute the error estimate and flag the cells
    ERROR_ESTIMATOR::Compute_refinement(cross_sections.Get_n_groups(),
        dof_handler,&parameters,group_flux,cells_to_refine,adjacent_cells,
        edge_to_refine);
    // Build the refined grid
    triangulation.Refine_mesh(cells_to_refine,adjacent_cells,edge_to_refine,
        projection);
    // Project the solution 
    unsigned int scalar_flux_size = projection.size();
    Epetra_Map* new_scalar_flux_map = new Epetra_Map(scalar_flux_size,0,comm);
    Epetra_MultiVector* new_group_flux = new Epetra_MultiVector(*new_scalar_flux_map,
        cross_sections.Get_n_groups());

    // Project group_flux on the new grid
    for (unsigned int g=0; g<cross_sections.Get_n_groups(); ++g)
    {
      for (unsigned int i=0; i<scalar_flux_size; ++i)
      {
        unsigned int n_old_vertices(projection[i].size());
        for (unsigned int j=0; j<n_old_vertices; ++j)
          (*new_group_flux)[g][i] += (*group_flux)[g][projection[i][j]];
        (*new_group_flux)[g][i] /= (double)(n_old_vertices);
      }
    }

    // Copy new_scalar_flux_map and new_group_flux to scalar_flux_map and
    // group_flux
    delete group_flux;
    delete scalar_flux_map;
    scalar_flux_map = new_scalar_flux_map;
    group_flux = new_group_flux;
    new_scalar_flux_map = NULL;
    new_group_flux = NULL;
    // Free the dof_handler
    delete dof_handler;
    // Build the new edges
    triangulation.Build_edges();
    // Create the new dof_handler
    dof_handler = new DOF_HANDLER(&triangulation,parameters,cross_sections);
  }

  for (unsigned int i=0; i<192; ++i)
    assert(fabs(solution[i]-(*group_flux)[0][i])<1e-4);

  MPI_Finalize();

  return 0;
}
