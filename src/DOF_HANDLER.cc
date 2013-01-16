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

#include "DOF_HANDLER.hh"

DOF_HANDLER::DOF_HANDLER(TRIANGULATION* triang,PARAMETERS const &param,
    CROSS_SECTIONS const &cross_sections) :
  n_sf_per_dir(0),
  n_dof(0),
  max_dof_per_cell(0),
  triangulation(triang)
{
  unsigned int start_dof(0);
  unsigned int start_reflective_dof(0);
  const unsigned int n_cells(triangulation->Get_n_cells());
  for (unsigned int i=0; i<n_cells; ++i)
  {
    const unsigned int n_vertices(triangulation->Get_n_vertices(i));
    d_vector x(n_vertices,0.);
    d_vector y(n_vertices,0.);
    for (unsigned int j=0; j<n_vertices; ++j)
    {
      x[j] = triangulation->Get_x_coord(i,j);
      y[j] = triangulation->Get_y_coord(i,j);
    }                                   
    FINITE_ELEMENT* fe;
    if (param.Get_fe_type()==bld)
      fe = new BLD(x,y);
    else
      fe = new PWLD(x,y);
    vector<EDGE*> cell_edges(n_vertices);
    bool cell_saf(false);
    for (unsigned int j=0; j<n_vertices; ++j)
    {
      unsigned int gid(triangulation->Get_edge_gid(i,j));
      EDGE* edge(triangulation->Get_edge(gid));
      if (edge->Get_cell_index(0)==i)
        edge->Set_lid(0,j);
      else
        edge->Set_lid(1,j);

      // Set the type of boundary condition on the edge
      if (edge->Get_edge_type()==bottom_boundary)
        edge->Set_bc_type(param.Get_bottom_bc_type());
      if (edge->Get_edge_type()==right_boundary)
        edge->Set_bc_type(param.Get_right_bc_type());
      if (edge->Get_edge_type()==top_boundary)
        edge->Set_bc_type(param.Get_top_bc_type());
      if (edge->Get_edge_type()==left_boundary)
        edge->Set_bc_type(param.Get_left_bc_type());

      // If the edge is on the reflective boundary, store the number of dof
      // associated to the cell
      if (edge->Is_interior()==false)
      {
        if (edge->Get_bc_type()==reflective)
        {
          n_sf_per_dir += fe->Get_dof_per_cell();
          cell_saf = true;
        }
      }

      cell_edges[j] = edge;
    }
    if (cell_saf==true)
    {
      saf_map_dof[i] = start_dof;
      saf_map_reflective_dof[i] = start_reflective_dof;
      start_reflective_dof += fe->Get_dof_per_cell();
    }
    // Create the cells and append them in mesh
    if (cross_sections.Sigma_e_exist()==false)
    {
      CELL* cell = new CELL(i,n_vertices,start_dof,start_dof+fe->Get_dof_per_cell(),
          param.Get_src(triangulation->Get_src_id(i)),
          cross_sections.Get_sigma_t(triangulation->Get_mat_id(i)),
          cross_sections.Get_sigma_s(triangulation->Get_mat_id(i)),cell_edges,fe);
      mesh.push_back(cell);
    }
    else
    {
      CELL* cell = new CELL(i,n_vertices,start_dof,start_dof+fe->Get_dof_per_cell(),
          param.Get_src(triangulation->Get_src_id(i)),
          cross_sections.Get_sigma_t(triangulation->Get_mat_id(i)),
          cross_sections.Get_sigma_e(triangulation->Get_mat_id(i)),
          cross_sections.Get_sigma_s(triangulation->Get_mat_id(i)),cell_edges,fe);
      mesh.push_back(cell);
    }
    start_dof += fe->Get_dof_per_cell();
    if (fe->Get_dof_per_cell()>max_dof_per_cell)
      max_dof_per_cell = fe->Get_dof_per_cell();
  }
  // Build the finite elements of each cell
  for (unsigned int i=0; i<n_cells; ++i)
  {
    FINITE_ELEMENT* fe(mesh[i]->Get_mutable_fe());
    n_dof += fe->Get_dof_per_cell();
    fe->Build_fe_1d();
    fe->Build_fe_2d();
    fe->Build_upwind_matrices(mesh[i],mesh);
  }

  Compute_exterior_normal();
}

DOF_HANDLER::~DOF_HANDLER()
{
  for (unsigned int i=0; i<mesh.size(); ++i)
  {
    delete mesh[i];
    mesh[i] = NULL;
  }
}

void DOF_HANDLER::Compute_sweep_ordering(vector<QUADRATURE*> &quad)
{
  const unsigned int n_quad(quad.size());
  sweep_order.resize(n_quad);
  most_n_bottom.resize(n_quad);
  most_n_right.resize(n_quad);
  most_n_top.resize(n_quad);
  most_n_left.resize(n_quad);

  if (triangulation->Get_cell_type()==rectangle)
  {
    const int n_x_cells(triangulation->Get_n_x_cells());
    const int n_y_cells(triangulation->Get_n_y_cells());
    for (unsigned int q=0; q<n_quad; ++q)
    {
      const unsigned int n_dir(quad[q]->Get_n_dir());
      for (unsigned int idir=0; idir<n_dir; ++idir)
      {
        // Sweep ordering for a direction
        ui_vector idir_sweep_order;
        // Pointer to omega for direction idir
        Teuchos::SerialDenseVector<int,double> const* const omega=
          quad[q]->Get_omega(idir);
        if ((*omega)[0]>0)
        {
          if ((*omega)[1]>0)
          {
            for (int y=0; y<n_y_cells; ++y)
              for (int x=0; x<n_x_cells; ++x)
                idir_sweep_order.push_back(x+y*n_x_cells);
          }
          else
          {
            for (int y=n_y_cells-1; y>-1; --y)
              for (int x=0; x<n_x_cells; ++x)
                idir_sweep_order.push_back(x+y*n_x_cells);
          }
        }
        else
        {
          if ((*omega)[1]>0)
          {
            for (int y=0; y<n_y_cells; ++y)
              for (int x=n_x_cells-1; x>-1; --x)
                idir_sweep_order.push_back(x+y*n_x_cells);
          }
          else
          {
            for (int y=n_y_cells-1; y>-1; --y)
              for (int x=n_x_cells-1; x>-1; --x)
                idir_sweep_order.push_back(x+y*n_x_cells);
          }
        }
        sweep_order[q].push_back(idir_sweep_order);
      }
      // Compute the normal directions
      Compute_most_normal_direction(quad[q],q);
    }
  }
  else
  {
    for (unsigned int q=0; q<n_quad; ++q)
    {
      const unsigned int n_dir(quad[q]->Get_n_dir());
      for (unsigned int idir=0; idir<n_dir; ++idir)
      {
        // Sweep ordering for a direction
        ui_vector idir_sweep_order;
        // Cells already in idir_sweep_order
        ui_set cells_done;
        // Edges ready to be put on idir_sweep_order
        ui_list incoming_edges;
        Teuchos::SerialDenseVector<int,double> omega(quad[q]->Get_omega_2d(idir));
        // Look for the edges with a known incoming flux on the boundary
        vector<EDGE>::iterator edge(triangulation->Get_edges_begin());
        vector<EDGE>::iterator edge_end(triangulation->Get_edges_end());
        for (; edge<edge_end; ++edge)
        {
          if (((edge->Get_edge_type()==left_boundary) && (omega(0)>0.)) ||
              ((edge->Get_edge_type()==right_boundary) && (omega(0)<0.)) ||
              ((edge->Get_edge_type()==bottom_boundary) && (omega(1)>0.)) ||
              ((edge->Get_edge_type()==top_boundary) && (omega(1)<0.)))
            incoming_edges.push_back(edge->Get_gid());
        }
        while (incoming_edges.size()!=0)
        {
          unsigned int index(0);
          EDGE* test_edge(triangulation->Get_edge(*(incoming_edges.begin())));
          CELL* cell;
          if (cells_done.count(test_edge->Get_cell_index(0))!=0)
          {
            cell = mesh[test_edge->Get_cell_index(1)];
            index = 1;
          }
          else
            cell = mesh[test_edge->Get_cell_index(0)];
          unsigned int pos(0);
          const unsigned int length(cell->Get_n_edges()-1);
          vector<bool> outgoing_tests(length,false);
          vector<bool> edge_parallel(length,false);
          vector<bool> incoming_tests(length,false);
          ui_vector edges_map(length,0.);
          vector<EDGE*>::iterator cell_edge(cell->Get_cell_edges_begin());
          vector<EDGE*>::iterator cell_edge_end(cell->Get_cell_edges_end());
          for (; cell_edge<cell_edge_end; ++cell_edge)
          {
            if ((*cell_edge)->Get_gid()!=test_edge->Get_gid())
            {
              // Check if the edge is outgoing
              unsigned int index_2(0);
              if ((*cell_edge)->Get_cell_index(index_2)!=
                  test_edge->Get_cell_index(index))
                index_2 = 1;
              if (omega.dot(*((*cell_edge)->Get_exterior_normal(index_2)))>=0.)
              {
                outgoing_tests[pos] = true;
                if (omega.dot(*((*cell_edge)->Get_exterior_normal(index_2)))==0.)
                  edge_parallel[pos] = true;
              }
              else
                if (count(incoming_edges.begin(),incoming_edges.end(),
                      (*cell_edge)->Get_gid())==1)
                  incoming_tests[pos] = true;
              edges_map[pos] = (*cell_edge)->Get_gid();
              ++pos;
            }
          }
          // The cell is going to be accepted if all the edges are outgoing, all
          // the edges but one are outgoing and the one which is not outgoing is
          // ready to be incoming, all the edges but two are outgoing and the
          // two which are not outgoing are ready to be incoming, etc.
          bool test(true);
          for (unsigned int i=0; i<length; ++i)
          {
            if ((outgoing_tests[i]==false) && (incoming_tests[i]==false))
            {
              test = false;
              break;
            }
          }
          if (test==true)
          {
            // The polygon is accepted
            idir_sweep_order.push_back(cell->Get_id());
            // Remove the edge from the list
            incoming_edges.pop_front();
            // If the other edges of the cell are outgoing, they become incoming
            // except if they are on the on the boundary. If they are incoming,
            // they are removed
            for (unsigned int i=0; i<length; ++i)
            {
              if (outgoing_tests[i]==true)
              {
                if (((triangulation->Get_edge(edges_map[i]))->Is_interior()==true) &&
                    edge_parallel[i]==false)
                  incoming_edges.push_back(edges_map[i]);
              }
              else
                incoming_edges.remove(edges_map[i]);
            }
            cells_done.insert(cell->Get_id());
          }
          else
          {
            // The polygon is rejected. The edge is put at the end of the list.
            unsigned int tmp(*(incoming_edges.begin()));
            incoming_edges.pop_front();
            incoming_edges.push_back(tmp);
          }
        }
        sweep_order[q].push_back(idir_sweep_order);
      }
      // Compute the normal directions
      Compute_most_normal_direction(quad[q],q);
    }
  }
}

void DOF_HANDLER::Compute_most_normal_direction(QUADRATURE* quad,unsigned int lvl)
{
  const unsigned int n_dir(quad->Get_n_dir());
  double bottom(0.);
  double right(0.);
  double top(0.);
  double left(0.);

  // Find the most normal direction
  for (unsigned int i=0; i<n_dir; ++i)
  {
    Teuchos::SerialDenseVector<int,double> const* omega(quad->Get_omega(i));
    if ((*omega)(1)>bottom)
      bottom = (*omega)(1);
    if ((*omega)(0)<right)
      right = (*omega)(0);
    if ((*omega)(1)<top)
      top = (*omega)(1);
    if ((*omega)(0)>left)
      left = (*omega)(0);
  }

  // Store the direction in a set
  for (unsigned int i=0; i<n_dir; ++i)
  {
    Teuchos::SerialDenseVector<int,double> const* omega(quad->Get_omega(i));
    if ((*omega)(1)==bottom)
      most_n_bottom[lvl].insert(i);
    if ((*omega)(0)==right)
      most_n_right[lvl].insert(i);
    if ((*omega)(1)==top)
      most_n_top[lvl].insert(i);
    if ((*omega)(0)==left)
      most_n_left[lvl].insert(i);
  }
}

void DOF_HANDLER::Compute_exterior_normal()
{
  vector<EDGE>::iterator edge(triangulation->Get_edges_begin());
  vector<EDGE>::iterator edge_end(triangulation->Get_edges_end());
  for (; edge<edge_end; ++edge)
  {
    Teuchos::SerialDenseVector<int,double> edge_vector(2);
    Teuchos::SerialDenseVector<int,double> n_edge(2);
    n_edge(0) = edge->Get_v1_y()-edge->Get_v0_y();
    n_edge(1) = -edge->Get_v1_x()+edge->Get_v0_x();
    n_edge *= 1./n_edge.normFrobenius();
    // Compute the sign, assume that the polygon is convex
    vector<EDGE*>::iterator cell_edge(
        mesh[edge->Get_cell_index(0)]->Get_cell_edges_begin());
    vector<EDGE*>::iterator cell_edge_end(
        mesh[edge->Get_cell_index(0)]->Get_cell_edges_end());
    for (; cell_edge<cell_edge_end; ++cell_edge)
    {
      if ((*cell_edge)->Get_gid()!=edge->Get_gid())
      {
        if ((((*cell_edge)->Get_v0_x()==edge->Get_v1_x()) &&
            ((*cell_edge)->Get_v0_y()==edge->Get_v1_y())) ||
            (((*cell_edge)->Get_v0_x()==edge->Get_v0_x()) &&
            ((*cell_edge)->Get_v0_y()==edge->Get_v0_y())))
        {
          edge_vector(0) = (*cell_edge)->Get_v1_x()-(*cell_edge)->Get_v0_x();
          edge_vector(1) = (*cell_edge)->Get_v1_y()-(*cell_edge)->Get_v0_y();
          break;
        }
        else
        {
          if ((((*cell_edge)->Get_v1_x()==edge->Get_v1_x()) &&
              ((*cell_edge)->Get_v1_y()==edge->Get_v1_y())) ||
              (((*cell_edge)->Get_v1_x()==edge->Get_v0_x()) &&
              ((*cell_edge)->Get_v1_y()==edge->Get_v0_y())))
          {
            edge_vector(0) = (*cell_edge)->Get_v0_x()-(*cell_edge)->Get_v1_x();
            edge_vector(1) = (*cell_edge)->Get_v0_y()-(*cell_edge)->Get_v1_y();
            break;
          }
        }
      }
    }
    if (n_edge.dot(edge_vector)>0.)
    {
      edge->Set_exterior_normal(1,n_edge);
      n_edge *= -1.;
      edge->Set_exterior_normal(0,n_edge);
    }
    else
    {
      edge->Set_exterior_normal(0,n_edge);
      n_edge *= -1.;
      edge->Set_exterior_normal(1,n_edge);
    }
  }
}
