/*
Copyright (c) 2012, Bruno Turcksin.

This file is part of Janus.

Janus is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
he Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Janus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Janus.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ERROR_ESTIMATOR.hh"

namespace ERROR_ESTIMATOR
{
  void Compute_refinement(const unsigned int n_groups,DOF_HANDLER *dof_handler,
      PARAMETERS const *const parameters,Epetra_MultiVector const* const group_flux,
      ui_set &cells_to_refine,ui_set &adjacent_cells,
      map<unsigned int,vector<vector<d_vector> > > &edge_to_refine)
  {
    const unsigned int size_error_estimate(dof_handler->Get_n_cells()*n_groups);
    const unsigned int n_cells(dof_handler->Get_n_cells());
    double max_error(0.);
    const double refinement_threshold(parameters->Get_refinement_threshold());
    d_vector error_estimate(size_error_estimate,0.);

    // Loop over the edges to compute the error estimate
    vector<EDGE>::iterator edge(dof_handler->Get_edges_begin());
    vector<EDGE>::iterator edge_end(dof_handler->Get_edges_end());
    for (; edge<edge_end; ++edge)
    {  
      if (edge->Is_interior()==true)
      {
        const unsigned int edge_lid_0(edge->Get_lid(0));
        const unsigned int edge_lid_1(edge->Get_lid(1));
        CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
        CELL* next_cell(dof_handler->Get_cell(edge->Get_cell_index(1)));
        const unsigned int cell_first_dof(cell->Get_first_dof());
        const unsigned int cell_last_dof(cell->Get_last_dof());
        const unsigned int next_cell_first_dof(next_cell->Get_first_dof());
        const unsigned int next_cell_last_dof(next_cell->Get_last_dof());
        unsigned int first_dof_m(cell_first_dof+edge_lid_0);
        unsigned int second_dof_m(first_dof_m+1);
        unsigned int first_dof_p(next_cell_first_dof+edge_lid_1);
        unsigned int second_dof_p(first_dof_p+1);

        if (second_dof_m==cell_last_dof)
          second_dof_m = cell_first_dof;
        if (second_dof_p==next_cell_last_dof)
          second_dof_p = next_cell_first_dof;

        for (unsigned int g=0; g<n_groups; ++g)
        {
          double error(edge->Get_length()/2.);

          error *= pow((*group_flux)[g][first_dof_p]-(*group_flux)[g][second_dof_m],
              2)+pow((*group_flux)[g][second_dof_p]-(*group_flux)[g][first_dof_m],2);

          error_estimate[cell->Get_id()+g*n_cells] += error;
          error_estimate[next_cell->Get_id()+g*n_cells] += error;
        }
      }
      else
      {
        const unsigned int edge_lid_0(edge->Get_lid(0));
        CELL* cell(dof_handler->Get_cell(edge->Get_cell_index(0)));
        const unsigned int cell_first_dof(cell->Get_first_dof());
        const unsigned int cell_last_dof(cell->Get_last_dof());
        unsigned int first_dof_m(cell_first_dof+edge_lid_0);
        unsigned int second_dof_m(first_dof_m+1);

        if (second_dof_m==cell_last_dof)
          second_dof_m = cell_first_dof;

        for (unsigned int g=0; g<n_groups; ++g)
        {
          double error(edge->Get_length()/2.);

          double inc_flux_norm(0.);
          if (edge->Get_edge_type()==bottom_boundary)
          {
            if (edge->Get_bc_type()==isotropic)
              inc_flux_norm = parameters->Get_inc_bottom(g);
          }
          if (edge->Get_edge_type()==right_boundary) 
          {
            if (edge->Get_bc_type()==isotropic)
              inc_flux_norm = parameters->Get_inc_right(g);
          }
          if (edge->Get_edge_type()==top_boundary) 
          {
            if (edge->Get_bc_type()==isotropic)
              inc_flux_norm = parameters->Get_inc_top(g);
          }
          if (edge->Get_edge_type()==left_boundary) 
          {
            if (edge->Get_bc_type()==isotropic)
              inc_flux_norm = parameters->Get_inc_left(g);
          }
          inc_flux_norm /= parameters->Get_weight_sum();


          error *= pow(inc_flux_norm-(*group_flux)[g][second_dof_m],2)+
            pow(inc_flux_norm-(*group_flux)[g][first_dof_m],2);

          error_estimate[cell->Get_id()+g*n_cells] += error;
        }
      }
    }

    vector<CELL*>::iterator cell(dof_handler->Get_mesh_begin());
    vector<CELL*>::iterator cell_end(dof_handler->Get_mesh_end());
    for (; cell<cell_end; ++cell)
    {
      const unsigned int i_min((*cell)->Get_first_dof());
      const unsigned int i_max((*cell)->Get_last_dof());
      const double n_i(i_max-i_min);
      for (unsigned int g=0; g<n_groups; ++g)
      {
        double average_flux(0.);
        for (unsigned int i=i_min; i<i_max; ++i)
          average_flux += fabs((*group_flux)[g][i]/n_i);
        error_estimate[(*cell)->Get_id()+g*n_cells] /= pow(average_flux,2);
      }
    }

    // Find the largest error and multiply it by the refinement threshold
    max_error = refinement_threshold*
      (*max_element(error_estimate.begin(),error_estimate.end()));

    // Loop over the cell to build cells_to_refine, adjacents_cells, and
    // edge_to_refine
    cell = dof_handler->Get_mesh_begin();
    TRIANGULATION const *const triangulation(dof_handler->Get_triangulation());
    for (; cell<cell_end; ++cell)
    {
      for (unsigned int g=0; g<n_groups; ++g)
      {
        if (error_estimate[(*cell)->Get_id()+g*n_cells]>max_error)
        {
          cells_to_refine.insert((*cell)->Get_id());
          vector<EDGE*>::iterator cell_edge((*cell)->Get_cell_edges_begin());
          vector<EDGE*>::iterator cell_edge_end((*cell)->Get_cell_edges_end());
          for (; cell_edge<cell_edge_end; ++cell_edge)
          {
            d_vector const* const cell_vertex_0((*cell_edge)->Get_v0());
            d_vector const* const cell_vertex_1((*cell_edge)->Get_v1());
            unsigned int index(2);
            if ((*cell)->Get_id()!=(*cell_edge)->Get_cell_index(0))
              ++index;

            if (((*cell_edge)->Is_interior()) && ((*cell_vertex_0)[index]==0.) &&
                ((*cell_vertex_1)[index]==0.))
            {
              CELL* upwind_cell(NULL);
              if ((*cell_edge)->Get_cell_index(0)==(*cell)->Get_id())
                upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(1));
              else
                upwind_cell = dof_handler->Get_cell((*cell_edge)->Get_cell_index(0));

              d_vector vertex(2,0.);
              vector<d_vector> vertices;
              for (unsigned int i=0; i<upwind_cell->Get_n_edges(); ++i)
              {
                const double current_vertex[2] = {triangulation->Get_x_coord(
                    upwind_cell->Get_id(),i),triangulation->Get_y_coord(
                    upwind_cell->Get_id(),i)};
                if (((current_vertex[0]==(*cell_edge)->Get_v0_x()) &&
                      (current_vertex[1]==(*cell_edge)->Get_v0_y())) ||
                    ((current_vertex[0]==(*cell_edge)->Get_v1_x()) &&
                     (current_vertex[1]==(*cell_edge)->Get_v1_y())))
                {
                  vertex[0] = current_vertex[0];
                  vertex[1] = current_vertex[1];
                  vertices.push_back(vertex);
                }
              } 

              ui_set::iterator list_it(adjacent_cells.find(upwind_cell->Get_id()));
              if (list_it==adjacent_cells.end())
              {
                vector<vector<d_vector> > refined_edge;
                refined_edge.push_back(vertices);

                edge_to_refine.insert(pair<unsigned int,vector<vector<d_vector> > >
                    (upwind_cell->Get_id(),refined_edge));
                adjacent_cells.insert(upwind_cell->Get_id());
              }
              else
                edge_to_refine[upwind_cell->Get_id()].push_back(vertices);
            }
          }
          break;
        }
      }
    }
  }
}
