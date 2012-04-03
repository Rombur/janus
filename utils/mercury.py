#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2012-04-02 15:33:48.890257
#----------------------------------------------------------------------------#

import CONVERT_INPUT

input_filename = "input.txt"
output_filename = "parameters.inp"

converter = CONVERT_INPUT.CONVERT_INPUT(input_filename,output_filename)
converter.Read_input()
converter.Convert_input_file()
