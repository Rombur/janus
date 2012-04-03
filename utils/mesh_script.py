#! /usr/bin/env python
#----------------------------------------------------------------------------#
# Python code
# Author: Bruno Turcksin
# Date: 2012-03-19 18:05:32.329462
#----------------------------------------------------------------------------#

import numpy as np
import random

x = [0.,1.,2.]
y = [0.,1.,2.]
mat_id = [0,1,2,3]
src_id = [3,2,1,0]
# mesh_type can be Z-mesh or random
mesh_type = 'Z-mesh'
# if a Z-mesh is used, another parameter needs to be defined
alpha = 0.4
output_filename = 'output'

nx = len(x)-1
ny = len(y)-1
n_cells = nx*ny
output = open(output_filename+'.txt','w')

if mesh_type=='random' :
  x_min = x[0]
  x_max = x[-1]
  y_min = y[0]
  y_max = y[-1]
  x_full = [[] for i in xrange(ny+1)]
  y_full = [[[] for j in xrange(nx+1)] for i in xrange(ny+1)]
  for i in xrange(ny+1) :
    x_full[i] = x[:]
  for i in xrange(ny+1) :
    for j in xrange(nx+1) :
      y_full[i][j] = y[i]

  x_plot = x_full[:][:]
  y_plot = y_full[:][:]
  for i in xrange(ny+1) :
    for j in xrange(1,nx) :
      x_plot[i][j] = x_full[i][j]+0.25*random.random()*(x_full[i][j]-x_full[i][j-1])
  for i in xrange(1,ny) :
    for j in xrange(nx+1) :
      y_plot[i][j] = y_full[i][j]+0.25*random.random()*(y_full[i][j]-y_full[i-1][j])
elif mesh_type=='Z-mesh' :
  x = np.linspace(0.,1.,21)
  nx = 20
  ny = 20
  mat_id = [0 for i in xrange(nx*ny)]
  src_id = [0 for i in xrange(nx*ny)]
  x_plot = [[] for i in xrange(ny+1)]
  y_plot= [[[] for j in xrange(nx+1)] for i in xrange(ny+1)]
  for i in xrange(ny+1) :
    x_plot[i] = x[:]
  for i in xrange(ny+1) :
    for j in xrange(nx+1) :
      if j/20.<=0.2 :
        if i<11 :
          y_plot[i][j] = i/10.*alpha
        elif i<20 :
          y_plot[i][j] = alpha+(i-10.)/10.*(1.-alpha)
        else :
          y_plot[i][j] = 1.
      elif j/20.<=0.35 :
        if i<11 :
          y_plot[i][j] = i/10.*(alpha+((j-4.)/3.)*(1.-2.*alpha))
        elif i<20 :
          y_plot[i][j] = alpha+(i-10.)/10.*((1.-alpha)+(j-4.)/3.*(2.*alpha-1.)+\
              (j-4.)/3.*(1.-alpha))
        else :
          y_plot[i][j] = 1.
      elif j/20.<=0.65 :
        if i<11 :
          y_plot[i][j] = (1.-alpha)*i/10.-(j-7.)/6.*(1.-2.*alpha)*i/10.
        elif i<20 :
          y_plot[i][j] = (i-10.)/10.*alpha+(1.-alpha)+(j-7.)/6.*(2.*alpha-1+\
              (i-10.)/10.*(1.-2.*alpha))
        else :
          y_plot[i][j] = 1.
      elif j/20.<=0.8 :
        if i<11 :
          y_plot[i][j] = i/10.*(alpha+((j-13.)/3.)*(1.-2.*alpha))
        elif i<20 :
          y_plot[i][j] = alpha*+(i-10.)/10.*((1.-alpha)+(j-13.)/3.*(2.*alpha-1.)+\
              (j-13.)/3.*(1.-alpha))
        else :
          y_plot[i][j] = 1.
      else :
        if i<11 :
          y_plot[i][j] = i/10.*(1.-alpha)
        elif i<20 :
          y_plot[i][j] = (1.-alpha)+i/10.*alpha
        else :
          y_plot[i][j] = 1.
  
for i in xrange(ny) :
  for j in xrange(nx) :
    vertex_1 = [x_plot[i][j],y_plot[i][j]]
    vertex_2 = [x_plot[i][j+1],y_plot[i][j+1]]
    vertex_3 = [x_plot[i+1][j+1],y_plot[i+1][j+1]]
    vertex_4 = [x_plot[i+1][j],y_plot[i+1][j]]
    output.write('4 '+str(vertex_1[0])+' '+str(vertex_1[1])+' ')
    output.write(str(vertex_2[0])+' '+str(vertex_2[1])+' ')
    output.write(str(vertex_3[0])+' '+str(vertex_3[1])+' ')
    output.write(str(vertex_4[0])+' '+str(vertex_4[1])+' ')
    output.write(str(mat_id[j+i*nx])+' '+str(src_id[j+i*nx])+'\n')

output.close()
