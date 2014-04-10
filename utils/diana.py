#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2012-03-16 17:53:00.340629
#----------------------------------------------------------------------------#

import MESH_GENERATOR
import TRIANGLE

input_filename = 'test.1'
output_filename = 'geometry'
mesh_type = 'review'
alpha = 0.3 # Used only for random mesh and z-mesh
x = [0.,1.] # Used only for random mesh
y = [0.,1.] # Used only for random mesh
side = 0.05 # Used only for hexagon mesh
n_layer = 30 # Used only for hexagon mesh
n_refinement = 20 # Used only for amr mesh
mat_id = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2] # Used only for hexagon mesh
src_id = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0] # Used only for hexagon mesh

triangle_mesh = ['polygon','triangle','quadrilateral']
if mesh_type in triangle_mesh :
  mesh = TRIANGLE.TRIANGLE(input_filename,output_filename)
  mesh.Read_triangle_output_files()
  if mesh_type=='polygon' :
    mesh.Generate_polygonal_mesh()
  elif mesh_type=='triangle' :
    mesh.Output_remaining_cells()
    mesh.Prepend_n_cells()
    mesh.Create_apollo_file()
  elif mesh_type=='quadrilateral' :
    mesh.Generate_quadrilateral_mesh()
else :
  mesh = MESH_GENERATOR.MESH_GENERATOR(output_filename)
  if mesh_type=='random' :
    mesh.Generate_random_mesh(x,y,alpha)
  elif mesh_type=='z-mesh' :
    mesh.Generate_z_mesh(alpha)
  elif mesh_type=='hexagon' :
    mesh.Generate_hexagon_mesh(side,n_layer,mat_id,src_id)
  elif mesh_type=='amr' :
    mesh.Generate_amr_mesh(n_refinement)
  elif mesh_type=='review' :
    mesh.Generate_review_mesh(100.,100.,10)
