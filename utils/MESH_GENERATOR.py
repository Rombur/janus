# Python code
# Author: Bruno Turcksin
# Date: 2012-05-16 16:43:53.803215

#----------------------------------------------------------------------------#
## Class MESH_GENERATOR                                                     ##
#----------------------------------------------------------------------------#

"""Create meshes which are not generated from triangle."""

import numpy as np

class MESH_GENERATOR(object) :
  """Create meshes that are not generated from triangle: random, z-mesh and
  hexagon."""

  def __init__(self,output_filename) :

    super(MESH_GENERATOR,self).__init__()

# Output filename
    self.output_filename = output_filename
    self.output_file = open(output_filename+'.inp','w')
    self.output_silo = open('silo_'+output_filename+'.txt','w')

    self.mesh = []

#----------------------------------------------------------------------------#

  def Write_cell(self,cell) :
    """Write the generated polygonal cell in the output file."""
                     
    output = str(cell[0])+' '
    for i in xrange(len(cell[1])) :
      output += str(cell[1][i][0])+' '+str(cell[1][i][1])+' '
    output += str(cell[-2])+' '+str(cell[-1])+'\n'

    self.output_file.write(output)

#----------------------------------------------------------------------------#

  def Generate_random_mesh(self,x,y,alpha) :
    """Generate a randomized mesh."""

    mat_id = [0]
    src_id = [0]
    nx = len(x)-1
    ny = len(y)-1
    n_cells = nx*ny
    self.output_file.write('polygon\n')
    self.output_file.write(str(n_cells)+'\n')

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
        x_plot[i][j] = x_full[i][j]+alpha*random.random()*\
            (x_full[i][j]-x_full[i][j-1])
    for i in xrange(1,ny) :
      for j in xrange(nx+1) :
        y_plot[i][j] = y_full[i][j]+alpha*random.random()*\
            (y_full[i][j]-y_full[i-1][j])

    for i in xrange(ny) :
      for j in xrange(nx) :
        vertex_1 = [x_plot[i][j],y_plot[i][j]]
        vertex_2 = [x_plot[i][j+1],y_plot[i][j+1]]
        vertex_3 = [x_plot[i+1][j+1],y_plot[i+1][j+1]]
        vertex_4 = [x_plot[i+1][j],y_plot[i+1][j]]
        vertices = []
        vertices.append(vertex_1)
        vertices.append(vertex_2)
        vertices.append(vertex_3)
        vertices.append(vertex_4)
        cell = [4,vertices,mat_id,src_id]
        self.mesh.append(cell)
        self.Write_cell(cell)

# Close output file
    self.output_file.close()

# Create a file readable by apollo
    self.Create_apollo_file()

#----------------------------------------------------------------------------#

  def Generate_z_mesh(self,alpha):
    """Generate a z-mesh given alpha factor [0,1]."""

    x = np.linspace(0.,1.,21)
    nx = 20
    ny = 20
    self.output_file('polygon\n')
    self.output_file('400\n')
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
            y_plot[i][j] = alpha +(j-4.)/3.*(1.-2.*alpha)+(i-10.)/10.*\
                (1.-alpha-(1.-2.*alpha)*(j-4.)/3.)
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
            y_plot[i][j] = alpha +(j-13.)/3.*(1.-2.*alpha)+(i-10.)/10.*\
                (1.-alpha-(1.-2.*alpha)*(j-13.)/3.)
          else :
            y_plot[i][j] = 1.
        else :
          if i<11 :
            y_plot[i][j] = i/10.*(1.-alpha)
          elif i<20 :
            y_plot[i][j] = (1.-alpha)+(i-10.)/10.*alpha
          else :
            y_plot[i][j] = 1.

    for i in xrange(ny) :
      for j in xrange(nx) :
        vertex_1 = [x_plot[i][j],y_plot[i][j]]
        vertex_2 = [x_plot[i][j+1],y_plot[i][j+1]]
        vertex_3 = [x_plot[i+1][j+1],y_plot[i+1][j+1]]
        vertex_4 = [x_plot[i+1][j],y_plot[i+1][j]]
        vertices = []
        vertices.append(vertex_1)
        vertices.append(vertex_2)
        vertices.append(vertex_3)
        vertices.append(vertex_4)
        cell = [4,vertices,mat_id,src_id]
        self.mesh.append(cell)
        self.Write_cell(cell)

# Close output file
    self.output_file.close()

# Create a file readable by apollo
    self.Create_apollo_file()

#----------------------------------------------------------------------------#

  def Generate_hexagon_mesh(self,side,n_layer,mat_id,src_id) :
    """Generate a mesh made of hexagonal cells of side length "side" and made 
    of "n_layer" layers."""

# Build the hexagonal cells
    self.Build_hexagonal_cells(side,n_layer,mat_id,src_id)
# Build the triangle cells on the bottom and the top of the domain
    self.Build_boundary_triangles(side,n_layer)
# Build the cells in the corner of the domain
    self.Build_cells_in_corners(side,n_layer)

# Close output file
    self.output_file.close()

# Prepend the number of cells now that the mesh is built
    self.Prepend_n_cells()

# Create a file readable by apollo
    self.Create_apollo_file()

#----------------------------------------------------------------------------#

  def Build_hexagonal_cells(self,side,n_layer,mat_id,src_id) :
    """Build the hexagonal cells of the mesh."""
    
    vertex_1 = [0.,-side]
    vertex_2 = [np.sqrt(3.)/2.*side,-0.5*side]
    vertex_3 = [np.sqrt(3.)/2.*side,0.5*side]
    vertex_4 = [0.,side]
    vertex_5 = [-np.sqrt(3.)/2.*side,0.5*side]
    vertex_6 = [-np.sqrt(3.)/2.*side,-0.5*side]
    y_offset = side
    for i in xrange(2*n_layer+1) :
      if i<n_layer+1 :
        j_max = n_layer+1+i
        x_offset = (n_layer+1-i)*np.sqrt(3.)/2.*side
      else :
        j_max = 2*n_layer+1-(i-n_layer)
        x_offset = (i-(n_layer-1))*np.sqrt(3.)/2.*side
      for j in xrange(0,j_max) :
        layer = self.Compute_layer(i,j,n_layer,j_max)
        vertices = []
        vertices.append([vertex_1[0]+x_offset,vertex_1[1]+y_offset])
        vertices.append([vertex_2[0]+x_offset,vertex_2[1]+y_offset])
        vertices.append([vertex_3[0]+x_offset,vertex_3[1]+y_offset])
        vertices.append([vertex_4[0]+x_offset,vertex_4[1]+y_offset])
        vertices.append([vertex_5[0]+x_offset,vertex_5[1]+y_offset])
        vertices.append([vertex_6[0]+x_offset,vertex_6[1]+y_offset])
        cell = [6,vertices,mat_id[layer],src_id[layer]]
        self.mesh.append(cell)
        self.Write_cell(cell)
        x_offset += np.sqrt(3.)*side
      y_offset += 1.5*side
    self.y_max = y_offset-0.5*side

#----------------------------------------------------------------------------#

  def Compute_layer(self,i,j,n_layer,j_max):

    if i==0 :
      layer = n_layer
    else :
      if i>n_layer :
        i = 2*n_layer-i
      if j<=i :
        layer = n_layer-j
      elif j>n_layer-1 :
        layer = n_layer-i+j-n_layer
      else :
        layer = n_layer-i

    return layer  

#----------------------------------------------------------------------------#

  def Build_boundary_triangles(self,side,n_layer) :
    """Build the triangular cells on the top and the bottom of the domain."""

    x_offset = (n_layer+1)*np.sqrt(3.)/2.*side
    for i in xrange(n_layer) :
      delta_x = np.sqrt(3.)*side
      vertices = []
      vertices.append([x_offset,0])
      vertices.append([x_offset+delta_x,0])
      vertices.append([x_offset+delta_x/2,0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
      vertices = []
      vertices.append([x_offset+delta_x,self.y_max])
      vertices.append([x_offset,self.y_max])
      vertices.append([x_offset+delta_x/2,self.y_max-0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
      x_offset += delta_x

#----------------------------------------------------------------------------#

  def Build_cells_in_corners(self,side,n_layer) :
    """Build the cells in the four corners of the domain."""

    x_max = (2.*n_layer+1.)*np.sqrt(3.)*side
    x_coord = [[] for i in xrange(n_layer+2)]
    x = 0.
    for i in xrange(n_layer+2) :
      x_coord[0].append(x)
      x += np.sqrt(3.)/2.*side
    for j in xrange(1,n_layer+2) :
      x_coord[j] = x_coord[0][0:n_layer+2-j]

    y_offset = 0.
    for i in xrange(n_layer) :
      for j in xrange(n_layer-i) :
# Build the rectangular cells in the bottom left corner
        vertices = []
        vertices.append([x_coord[i][j],y_offset])
        vertices.append([x_coord[i][j+1],y_offset])
        vertices.append([x_coord[i][j+1],y_offset+0.5*side])
        vertices.append([x_coord[i][j],y_offset+0.5*side])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the rectangular cells in the bottom right corner
        vertices = []
        vertices.append([x_max-x_coord[i][j+1],y_offset])
        vertices.append([x_max-x_coord[i][j],y_offset])
        vertices.append([x_max-x_coord[i][j],y_offset+0.5*side])
        vertices.append([x_max-x_coord[i][j+1],y_offset+0.5*side])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the triangular cells in the bottom left corner
      vertices = []
      vertices.append([x_coord[i][-2],y_offset])
      vertices.append([x_coord[i][-1],y_offset])      
      vertices.append([x_coord[i+1][-1],y_offset+0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
# Build the triangular cells in the bottom right corner
      vertices = []
      vertices.append([x_max-x_coord[i][-1],y_offset])
      vertices.append([x_max-x_coord[i][-2],y_offset])      
      vertices.append([x_max-x_coord[i+1][-1],y_offset+0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
      y_offset += 0.5*side
      for j in xrange(n_layer-i) :
# Build the rectangular cells in the bottom left corner
        vertices = []
        vertices.append([x_coord[i][j],y_offset])
        vertices.append([x_coord[i][j+1],y_offset])
        vertices.append([x_coord[i][j+1],y_offset+side])
        vertices.append([x_coord[i][j],y_offset+side])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the rectangular cells in the bottom right corner
        vertices = []
        vertices.append([x_max-x_coord[i][j+1],y_offset])
        vertices.append([x_max-x_coord[i][j],y_offset])
        vertices.append([x_max-x_coord[i][j],y_offset+side])
        vertices.append([x_max-x_coord[i][j+1],y_offset+side])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
      y_offset += side
# Build the triangular cell in the bottom left corner
    vertices = []
    vertices.append([x_coord[-2][-2],y_offset])
    vertices.append([x_coord[-2][-1],y_offset])      
    vertices.append([x_coord[-1][-1],y_offset+0.5*side])
    cell = [3,vertices,0,0]
    self.mesh.append(cell)
    self.Write_cell(cell)
# Build the triangular cell in the bottom right corner
    vertices = []
    vertices.append([x_max-x_coord[-2][-1],y_offset])
    vertices.append([x_max-x_coord[-2][-2],y_offset])      
    vertices.append([x_max-x_coord[-1][-1],y_offset+0.5*side])
    cell = [3,vertices,0,0]
    self.mesh.append(cell)
    self.Write_cell(cell)

    y_offset = self.y_max
    for i in xrange(n_layer) :
      for j in xrange(n_layer-i) :
# Build the rectangular cells in the top left corner
        vertices = []
        vertices.append([x_coord[i][j],y_offset-0.5*side])
        vertices.append([x_coord[i][j+1],y_offset-0.5*side])
        vertices.append([x_coord[i][j+1],y_offset])
        vertices.append([x_coord[i][j],y_offset])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the rectangular cells in the top right corner
        vertices = []
        vertices.append([x_max-x_coord[i][j+1],y_offset-0.5*side])
        vertices.append([x_max-x_coord[i][j],y_offset-0.5*side])
        vertices.append([x_max-x_coord[i][j],y_offset])
        vertices.append([x_max-x_coord[i][j+1],y_offset])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the triangular cells in the top left corner
      vertices = []
      vertices.append([x_coord[i][-2],y_offset])
      vertices.append([x_coord[i][-1],y_offset])      
      vertices.append([x_coord[i][-2],y_offset-0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
# Build the triangular cells in the top right corner
      vertices = []
      vertices.append([x_max-x_coord[i][-2],y_offset])
      vertices.append([x_max-x_coord[i][-1],y_offset])      
      vertices.append([x_max-x_coord[i][-2],y_offset-0.5*side])
      cell = [3,vertices,0,0]
      self.mesh.append(cell)
      self.Write_cell(cell)
      y_offset -= 0.5*side
      for j in xrange(n_layer-i) :
# Build the rectangular cells in the top left corner
        vertices = []
        vertices.append([x_coord[i][j],y_offset-side])
        vertices.append([x_coord[i][j+1],y_offset-side])
        vertices.append([x_coord[i][j+1],y_offset])
        vertices.append([x_coord[i][j],y_offset])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
# Build the rectangular cells in the top right corner
        vertices = []
        vertices.append([x_max-x_coord[i][j+1],y_offset-side])
        vertices.append([x_max-x_coord[i][j],y_offset-side])
        vertices.append([x_max-x_coord[i][j],y_offset])
        vertices.append([x_max-x_coord[i][j+1],y_offset])
        cell = [4,vertices,0,0]
        self.mesh.append(cell)
        self.Write_cell(cell)
      y_offset -= side
# Build the triangular cell in the top left corner
    vertices = []
    vertices.append([x_coord[-2][-2],y_offset])
    vertices.append([x_coord[-2][-1],y_offset])      
    vertices.append([x_coord[-2][-2],y_offset-0.5*side])
    cell = [3,vertices,0,0]
    self.mesh.append(cell)
    self.Write_cell(cell)
# Build the triangular cell in the top right corner
    vertices = []
    vertices.append([x_max-x_coord[-2][-2],y_offset])
    vertices.append([x_max-x_coord[-2][-1],y_offset])      
    vertices.append([x_max-x_coord[-2][-2],y_offset-0.5*side])
    cell = [3,vertices,0,0]
    self.mesh.append(cell)
    self.Write_cell(cell)

#----------------------------------------------------------------------------#

  def Create_apollo_file(self) :
    """Create the file readable by apollo which will create the file used by
    Visit."""

    n_cells = len(self.mesh)
    n_dof = 0
    offset = np.zeros(n_cells+1)
    repartition = np.zeros((8,1))

    index = 1
    for cell in self.mesh :
      n_dof += len(cell[1])
      offset[index] = n_dof
      index += 1

    self.output_silo.write(str(int(n_cells))+' ')
    self.output_silo.write(str(int(n_dof))+' ')
    self.output_silo.write(str(int(2*n_dof))+'\n')

    for off in offset :
      self.output_silo.write(str(int(off))+' ')
    self.output_silo.write('\n')

# Write vertices coordinates
    for cell in self.mesh :
      for vertex in cell[1] :
        self.output_silo.write(str(vertex[0])+' '+str(vertex[1])+'\n')
      repartition[len(cell[1])-3] += 1

# Write the material id
    for cell in self.mesh :
      for vertex in cell[1] :
        self.output_silo.write(str(int(cell[2]))+'\n')

# Write the source id
    for cell in self.mesh :
      for vertex in cell[1] :
        self.output_silo.write(str(int(cell[3]))+'\n')

    print 'Number of cells: ',n_cells
    print 'Repartition of the cells:'
    print '-triangles: ',repartition[0]
    print '-quadrilaterals: ',repartition[1]
    print '-pentagons: ',repartition[2]
    print '-hexagons: ',repartition[3]
    print '-heptagons: ',repartition[4]
    print '-octagons: ',repartition[5]
    print '-nonagons: ',repartition[6]
    print '-decagons: ',repartition[7]
    print 'Number of degrees of freedom: ',n_dof

#----------------------------------------------------------------------------#

  def Prepend_n_cells(self) :
    """Prepend the number of cells to the output file."""

    self.output_file = open(self.output_filename+'.inp','r+')
    old = self.output_file.read()
    self.output_file.seek(0)
    self.output_file.write('polygon\n')
    self.output_file.write(str(len(self.mesh))+'\n'+old)
    self.output_file.close()
