#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2012-03-16 17:53:00.340629
#----------------------------------------------------------------------------#

import TRIANGLE

input_filename = 'test.1'
output_filename = 'output'
mesh = 'quadrilateral'

triangle = TRIANGLE.TRIANGLE(input_filename,output_filename)
triangle.Read_triangle_output_files()
if mesh=='polygon' :
  triangle.Generate_polygonal_mesh()
elif mesh=='triangle' :
  triangle.Output_remaining_cells()
  triangle.Create_apollo_file()
elif mesh=='quadrilateral' :
  triangle.Generate_quadrilateral_mesh()
