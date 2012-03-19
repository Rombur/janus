#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2012-03-16 17:53:00.340629
#----------------------------------------------------------------------------#

import TRIANGLE

input_filename = 'test.1'
output_filename = 'output'

triangle = TRIANGLE.TRIANGLE(input_filename,output_filename)
triangle.Read_triangle_output_files()
triangle.Generate_polygonal_mesh()
