# Python code
# Author: Bruno Turcksin
# Date: 2012-03-15 22:11:57.111335

#----------------------------------------------------------------------------#
## Class TRIANGLE                                                           ##
#----------------------------------------------------------------------------#

"""Read the output files created by Triangle and build a polygonal mesh."""

import numpy as np
import utils

class TRIANGLE(object) :
  """Read the output files created by Triangle and build a polygonal mesh. If
  a triangle has an angle superior to pi/2, the generated polygon can be non
  convex. Therefore these triangles are not used to create polygons."""

  def __init__(self,input_filename,output_filename) :

    super(TRIANGLE,self).__init__()

# The nodes file contains:
#   - first line: number of nodes, number of dimension, dimension (always 2), 
#   number of attributes, number of boundary markers (0 or 1)
#   - remaining lines: vertex id, x, y, [attributes], [boundary marker]
    self.node_filename = input_filename+'.node'

# The elements file contains:
#   - first line: number of triangles, nodes per triangle (always 3), number of
#   attributes (always 1)
#   - remaining lines: triangle id, node, node, node, attribute
    self.elem_filename = input_filename+'.ele'

# Output filename
    self.output_filename = output_filename
    self.output_file = open(output_filename+'.inp','w')
    self.output_silo = open('silo_'+output_filename+'.txt','w')

    self.mesh = []
    self.used_triangles = set()

#----------------------------------------------------------------------------#

  def Read_triangle_output_files(self) :
    """Read the triangle output files."""

    self.available_vertices = []
# Read the elements file
    input_file = open(self.elem_filename,'r')
    line = input_file.readline()
    self.n_triangles,read = utils.Read_float(line,1)
    if read==False :
      utils.Abort("Problem reading the number of triangles.")
    self.n_triangles = int(self.n_triangles[0])
    self.triangles = [[[] for j in xrange(2)]\
        for i in xrange(self.n_triangles)] 
    for i in xrange(self.n_triangles) :
      line = input_file.readline()
      values,read = utils.Read_float(line,5)
      if read==False :
        utils.Abort("Problem reading the triangle node ids.")
      vertices = [values[1],values[2],values[3]]
      attribute = values[4]
      self.triangles[int(values[0])][0] = vertices
      self.triangles[int(values[0])][1] = attribute  
    input_file.close()

# Read the nodes file
    input_file = open(self.node_filename,'r')
    line = input_file.readline()
    self.n_vertices,read = utils.Read_float(line,1)
    if read==False :
      utils.Abort("Problem reading the number of vertices.")
    self.n_vertices = int(self.n_vertices[0])
    self.vertices = [[[] for j in xrange(3) ]\
        for i in xrange(self.n_vertices)]
    for i in xrange(self.n_vertices) :
      line = input_file.readline()
      values,read = utils.Read_float(line,4)
      if read==False :
        utils.Abort("Problem reading the vertex coordinates.")
      coord = np.array([values[1],values[2]])
      associated_triangles = []
      for j in xrange(self.n_triangles) :
        if values[0] in self.triangles[j][0] :
          associated_triangles.append(j)
      self.vertices[i][0] = associated_triangles
      self.vertices[i][1] = coord
      self.vertices[i][2] = values[3]
      self.available_vertices.append(i)
    input_file.close()

#----------------------------------------------------------------------------#

  def Generate_polygonal_mesh(self) :
    """Generate a polygonal mesh."""

    while len(self.available_vertices)>0 :
# Take the next available vertex
      triangle_ids = self.vertices[self.available_vertices[0]][0]

# All the triangles need to have the same attribute, i.e. same source id and
# same material id
      rejected = False
      if len(triangle_ids)==1 :
        rejected = True
      for i in xrange(len(triangle_ids)-1) :
        id_1 = triangle_ids[i]
        id_2 = triangle_ids[i+1]
        if self.triangles[id_1][1]!=self.triangles[id_2][1] :
          rejected = True
      
      if rejected==True :
        self.available_vertices.pop(0)
      else :
# Check that the triangle don't have an angle larger than pi/2 and that they
# don't form a corner of the mesh
        rejected = self.Check_triangles(triangle_ids)
        if rejected==True :
          self.available_vertices.pop(0)
        else :
          cell = self.Generate_cell(triangle_ids)
          self.mesh.append(cell)
# Output the generated cell
          self.Write_cell(cell)

# Output the remaining triangular cells
    self.Output_remaining_cells()

# Close output file
    self.output_file.close()

# Prepend the number of cells now that the mesh is built
    self.Prepend_n_cells()

# Create a file readable by apollo
    self.Create_apollo_file()

#----------------------------------------------------------------------------#

  def Check_triangles(self,triangle_ids) :
    """Check that the triangles don't have an angle larger than pi/2 or that
    they form a corner of the mesh (need to suppress this condition later). If 
    they do the vertex is rejected."""
    
    rejected = False
    corner_vertices = []
    corner_ids = set()
    common_coord = self.vertices[self.available_vertices[0]][1]
    for i in xrange(len(triangle_ids)) :
      vertices_id = self.triangles[triangle_ids[i]][0]
      vertex_1 = self.vertices[int(vertices_id[0])][1]
      vertex_2 = self.vertices[int(vertices_id[1])][1]
      vertex_3 = self.vertices[int(vertices_id[2])][1]
      a = np.sqrt((vertex_2[0]-vertex_3[0])**2+(vertex_2[1]-vertex_3[1])**2)
      b = np.sqrt((vertex_1[0]-vertex_3[0])**2+(vertex_1[1]-vertex_3[1])**2)
      c = np.sqrt((vertex_2[0]-vertex_1[0])**2+(vertex_2[1]-vertex_1[1])**2)
      alpha = np.arccos((b**2+c**2-a**2)/(2.*b*c))
      beta = np.arccos((a**2+c**2-b**2)/(2.*a*c))
      gamma = np.arccos((a**2+b**2-c**2)/(2.*a*b))
      if (alpha-np.pi/2)>1e-15 :
        if not np.equal(vertex_1,common_coord).all() :
          rejected = True
      if (beta-np.pi/2)>1e-15 :
        if not np.equal(vertex_2,common_coord).all() :
          rejected = True
      if (gamma-np.pi/2)>1e-15 :
        if not np.equal(vertex_3,common_coord).all() :
          rejected = True
# Check that the vertices don't form a corner. Each vertex can only be present
# once in the list.
      if self.vertices[int(vertices_id[0])][2]==1 and\
          vertices_id[0] not in corner_ids :
        corner_vertices.append(vertex_1)
        corner_ids.add(vertices_id[0])
      if self.vertices[int(vertices_id[1])][2]==1 and\
          vertices_id[1] not in corner_ids :
        corner_vertices.append(vertex_2)
        corner_ids.add(vertices_id[1])
      if self.vertices[int(vertices_id[2])][2]==1 and\
          vertices_id[2] not in corner_ids :
        corner_vertices.append(vertex_3)
        corner_ids.add(vertices_id[2])

    if len(corner_vertices)>=3 :
      if corner_vertices[0][0]==corner_vertices[1][0] :
        j = 0
      elif corner_vertices[0][1]==corner_vertices[1][1] : 
        j = 1
      else :
        rejected = True
      if rejected==False :      
        for i in xrange(2,len(corner_vertices)) :
          if corner_vertices[i][j]!=corner_vertices[i-1][j] :
            rejected = True
    
    return rejected

#----------------------------------------------------------------------------#

  def Generate_cell(self,triangle_ids) :
    """Generate the cell associated to the ids in triangles_ids, remove the
    triangle ids from the remaining cells and remove the vertex ids from the
    available vertices."""

    cell = []
# Reorder the vertices of the triangles and delete the commom vertex
    reordered_vertices = self.Reorder_vertices(triangle_ids)
# Compute the material id and the source id
    mat_id = int(self.triangles[triangle_ids[0]][1]%100)
    src_id = int((self.triangles[triangle_ids[0]][1]-mat_id)/100)
# Build the cell
    cell = [len(reordered_vertices),reordered_vertices,mat_id,src_id]
# Remove the vertices from the available vertices
    for i in xrange(len(triangle_ids)) :
      for j in xrange(3) :
        value = self.triangles[int(triangle_ids[i])][0][j]
        if value in self.available_vertices :
          index = self.available_vertices.index(value)
          self.available_vertices.pop(index)    
# Store the triangle ids in the used_triangles
    for triangle_id in triangle_ids : 
      self.used_triangles.add(triangle_id)
                                    
    return cell

#----------------------------------------------------------------------------#

  def Reorder_vertices(self,triangle_ids) :
    """Reorder the vertices of the triangles and delete the common vertex."""
    
    vertices = [[] for i in xrange(len(triangle_ids))]
    common_coord = self.vertices[self.available_vertices[0]][1]
    for i in xrange(len(triangle_ids)) :
      coord_1 = self.vertices[int(self.triangles[int(triangle_ids[i])][0][0])][1]
      coord_2 = self.vertices[int(self.triangles[int(triangle_ids[i])][0][1])][1]
      coord_3 = self.vertices[int(self.triangles[int(triangle_ids[i])][0][2])][1]
      if np.equal(coord_1,common_coord).all() :
        vertices[i] = [coord_2,coord_3]
      elif np.equal(coord_2,common_coord).all() :
        vertices[i] = [coord_3,coord_1]
      else :
        vertices[i] = [coord_1,coord_2]

# The while loop is used when the cells are on the boundary
    done_reordering = False
    while done_reordering==False :
      reordered_vertices = []
      reordered_vertices.append(vertices[0][0])
      reordered_vertices.append(vertices[0][1])
      reordering = True
      counter = 0
      tmp_remaining_ids = range(1,len(vertices))
      remaining_ids = tmp_remaining_ids[:]
      while len(remaining_ids)>0 :
        for j in remaining_ids :
          if np.equal(vertices[j][0],reordered_vertices[-1]).all() :
            reordered_vertices.append(vertices[j][1])
            tmp_remaining_ids.pop(tmp_remaining_ids.index(j))
        counter += 1
        if counter>=2*len(vertices) :
          tmp = vertices[0]
          for j in xrange(len(vertices)-1) :
            vertices[j] = vertices[j+1]
          vertices[-1] = tmp
          break
        remaining_ids = tmp_remaining_ids[:]
      if len(remaining_ids)==0 :
        done_reordering = True
# If the new cell is inside of the domain, the last vertex is the same than
# the first one and it has to be reomved. If the cell is on the boundary, 
# there is nothing to do. 
    if np.equal(reordered_vertices[0],reordered_vertices[-1]).all() :
      reordered_vertices = reordered_vertices[0:-1]
    
    return reordered_vertices

#----------------------------------------------------------------------------#

  def Write_cell(self,cell) :
    """Write the generated polygonal cell in the output file."""
                     
    output = str(cell[0])+' '
    for i in xrange(len(cell[1])) :
      output += str(cell[1][i][0])+' '+str(cell[1][i][1])+' '
    output += str(cell[-2])+' '+str(cell[-1])+'\n'

    self.output_file.write(output)

#----------------------------------------------------------------------------#

  def Output_remaining_cells(self) :
    """Output the triangular cells which could not be associated with any other
    triangular cell."""
  
    for i in xrange(self.n_triangles) :
      if i not in self.used_triangles :
        output = '3 '
        vertices = []
        for j in xrange(3) :
          vertex = self.vertices[int(self.triangles[i][0][j])][1]
          vertices.append(vertex)
          output += str(vertex[0])+' '+str(vertex[1])+' '
        mat_id = self.triangles[i][1]%100
        src_id = (self.triangles[i][1]-mat_id)/100
        output += str(int(mat_id))+' '+str(int(src_id))+'\n'
        self.mesh.append([3,vertices,mat_id,src_id])

        self.output_file.write(output)

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

    print 'Initial number of triangles: ',self.n_triangles
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
    print 'Initial number of degrees of freedom: ',3*self.n_triangles
    print 'Number of degrees of freedom: ',n_dof

#----------------------------------------------------------------------------#

  def Generate_quadrilateral_mesh(self) :
    """Generate a polygonal mesh."""

    self.output_file.write('polygon\n')
    self.output_file.write(str(3*len(self.triangles))+'\n')

    for triangle in self.triangles :
      mat_id = int(triangle[1]%100)
      src_id = int((triangle[1]-mat_id)/100)
      vertex_1 = self.vertices[int(triangle[0][0])][1]
      vertex_2 = self.vertices[int(triangle[0][1])][1]
      vertex_3 = self.vertices[int(triangle[0][2])][1]
      barycenter,mid_1,mid_2,mid_3 = self.Compute_barycenter_and_midpoints(triangle)
# Build the first cell in the triangle
      vertices = []
      vertices.append(vertex_1)
      vertices.append(mid_1)
      vertices.append(barycenter)
      vertices.append(mid_3)
      cell = [4,vertices,mat_id,src_id]
      self.mesh.append(cell)
      self.Write_cell(cell)
# Build the second cell in the triangle
      vertices = []
      vertices.append(vertex_2)
      vertices.append(mid_2)
      vertices.append(barycenter)
      vertices.append(mid_1)
      cell = [4,vertices,mat_id,src_id]
      self.mesh.append(cell)
      self.Write_cell(cell)
# Build the third cell in the triangle
      vertices = []
      vertices.append(vertex_3)
      vertices.append(mid_3)
      vertices.append(barycenter)
      vertices.append(mid_2)
      cell = [4,vertices,mat_id,src_id]
      self.mesh.append(cell)
      self.Write_cell(cell)

# Close output file
    self.output_file.close()

# Create a file readable by apollo
    self.Create_apollo_file()

#----------------------------------------------------------------------------#

  def Compute_barycenter_and_midpoints(self,triangle) :
    """Compute the coordinates of the barycenter and of the midpoints of each
    edge."""

    vertex_1 = self.vertices[int(triangle[0][0])][1]
    vertex_2 = self.vertices[int(triangle[0][1])][1]
    vertex_3 = self.vertices[int(triangle[0][2])][1]

    barycenter = [(vertex_1[0]+vertex_2[0]+vertex_3[0])/3.,\
        (vertex_1[1]+vertex_2[1]+vertex_3[1])/3.]
    midpoint_1 = [(vertex_1[0]+vertex_2[0])/2.,(vertex_1[1]+vertex_2[1])/2.] 
    midpoint_2 = [(vertex_2[0]+vertex_3[0])/2.,(vertex_2[1]+vertex_3[1])/2.] 
    midpoint_3 = [(vertex_1[0]+vertex_3[0])/2.,(vertex_1[1]+vertex_3[1])/2.] 

    return barycenter,midpoint_1,midpoint_2,midpoint_3

#----------------------------------------------------------------------------#

  def Prepend_n_cells(self) :
    """Prepend the number of cells to the output file."""

    self.output_file = open(self.output_filename+'.inp','r+')
    old = self.output_file.read()
    self.output_file.seek(0)
    self.output_file.write('polygon\n')
    self.output_file.write(str(len(self.mesh))+'\n'+old)
    self.output_file.close()
