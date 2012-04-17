# Python code
# Author: Bruno Turcksin
# Date: 2012-04-02 15:36:56.169688

#----------------------------------------------------------------------------#
## Class CONVERT_INPUT                                                      ##
#----------------------------------------------------------------------------#

"""Convert the python input file to c++ input file."""

import numpy as np

class CONVERT_INPUT(object) :
  """Convert the input file to a input file readable by janus."""

  def __init__(self,input_filename,output_filename) :

    super(CONVERT_INPUT,self).__init__()
    self.input_filename = input_filename
    self.output_filename = output_filename

#----------------------------------------------------------------------------#

  def Search(self,keyword) :
    """Search for the given keyword. The beginning and the end of the
    interesting section are saved in self.begin and self.end"""

    begin = "#BEGIN "
    end = "#END "

# Search the beginning of the section
    self.begin = self.data.find(begin+keyword)
    if self.begin==-1 :
      self.Abort("Cannot find "+begin+keyword)
    else :
      self.begin += len(begin+keyword)+1

# Search the end of the section
    self.end = self.data.find(end+keyword)
    if self.end==-1 :
      self.Abort("Cannot find "+end+keyword)

# Check for an error
    error = self.data[self.begin:self.end-1].find(end)
    if error!=-1 :
      self.abort("Nested sections.")

#----------------------------------------------------------------------------#

  def Abort(self,message) :
    print (message)
    print ("The program aborted.")
    exit()
#----------------------------------------------------------------------------#

  def String_to_np(self) :
    """Convert a string to a numpy array."""
      
    data_list = self.data[self.begin:self.end-1].split()
    data_array = np.zeros(len(data_list))
    for i in xrange(0,len(data_list)) :
      data_array[i] = data_list[i]

    return data_array  

#----------------------------------------------------------------------------#

  def Read_next(self,n) :
    """Read the next n values. If n=1 return the value otherwise return a
    numpy array."""

    self.pos = self.begin
    if n!=1 :
      numbers = np.zeros((n))
    for i in xrange(0,n) :
      value = ' '
# Skip the blanks
      while self.data[self.pos]==' ' or self.data[self.pos]==',' :
        self.pos += 1
# Read the value(s)
      while self.data[self.pos]!=' ' :
        if self.data[self.pos]!='\n' and self.data[self.pos]!=',' :
          value += self.data[self.pos]
          self.pos += 1
        else :
          break
      if value!=' ' :
        if n==1 :
          numbers = float(value)
        else :
          numbers[i] = float(value)
      else :
        utils.Abort('Cannot read element.')

    return numbers    

#----------------------------------------------------------------------------#

  def Read_input(self) :
    """Read the input file and store the data."""

    input_file = open(self.input_filename,'r')
    self.data = input_file.read()

# Read the solver type
    self.Search("SOLVER")
    self.solver = self.data[self.begin:self.end-1].lower()

# Read the tolerance
    self.Search("TOLERANCE")
    self.tolerance = self.Read_next(1)

# Read the maximum number of iteration
    self.Search("MAX ITER")
    self.max_iter = int(self.Read_next(1))

# Read the verbosity level of the code
    self.Search("VERBOSE")
    self.verbose = int(self.Read_next(1))

# Read if Fokker-Planck cross sections are used
    self.Search("FOKKER-PLANCK XS")
    self.fokker_planck = self.data[self.begin:self.end-1]

# Read if transport correction is used
    self.Search("TRANSPORT CORRECTION")
    self.transport_correction = self.data[self.begin:self.end-1].lower()

# If the transport correction is used, check if the optimal transport
# correction is used
    if self.transport_correction=="true" :
      self.Search("OPTIMAL TRANSPORT CORRECTION")
      self.optimal_correction = self.data[self.begin:self.end-1]

# Read if angular multigrid is used
    self.Search("ANGULAR MULTIGRID")
    self.angular_multigrid = self.data[self.begin:self.end-1]

# Read if MIP is used
    self.Search("MIP")
    self.mip = self.data[self.begin:self.end-1].lower()

# If MIP is used, read the type of solver for MIP
    if self.mip=="true" :
      self.Search("MIP SOLVER")
      self.mip_solver = self.data[self.begin:self.end-1].lower()
      if self.mip_solver=="cg_ml" :
        self.Search("AGGREGATION")
        self.aggregation = self.data[self.begin:self.end-1].lower()

# Read the quadrature type
    self.Search("QUADRATURE")
    self.quadrature = self.data[self.begin:self.end-1].lower()

# Read if the quadrature is Galerkin
    self.Search("GALERKIN")
    self.galerkin = self.data[self.begin:self.end-1]

# Read L_max
    self.Search("L MAX")
    self.L_max = int(self.Read_next(1))

# Read Sn order
    self.Search("SN ORDER")
    self.sn_order = int(self.Read_next(1))

# Read the finite elements type
    self.Search("FE")
    self.fe = self.data[self.begin:self.end-1].lower()

# Read the intensity of the source
    self.Search("SOURCE")
    self.source = self.String_to_np()

# Read the intensity of the bottom incoming flux
    self.Search("BOTTOM INCOMING FLUX")
    self.inc_bottom = self.Read_next(1)

# Read the intensity of the right incoming flux
    self.Search("RIGHT INCOMING FLUX")
    self.inc_right = self.Read_next(1)

# Read the intensity of the top incoming flux
    self.Search("TOP INCOMING FLUX")
    self.inc_top = self.Read_next(1)

# Read the intensity of the left incoming flux
    self.Search("LEFT INCOMING FLUX")
    self.inc_left = self.Read_next(1)

# Read the total cross sections
    self.Search("TOTAL XS")
    self.total_xs = self.String_to_np()

# If Fokker-Planck cross section are used, read alpha
    if self.fokker_planck=="true" :
      self.Search("ALPHA")
      self.scattering_xs = self.String_to_np()
    else :
      self.Search("SCATTERING XS")
      self.scattering_xs = self.String_to_np()

    input_file.close()

#----------------------------------------------------------------------------#

  def Convert_input_file(self) :
    """Write the input for janus."""

    output_file = open(self.output_filename,'w')
    
# Write the solver type
    if self.solver=="si" :
      output_file.write("SI ")
    elif self.solver=="bicgstab" :
      output_file.write("BiCGSTAB ")
    elif self.solver=="gmres_condnum" :
      output_file.write("GMRES_CONDNUM ")
    elif self.solver=="gmres" :
      output_file.write("GMRES ")
    else :
      self.Abort("Unknown solver.")

# Write the tolerance
    output_file.write(str(self.tolerance)+" ")

# Write the maximum number of iterations
    output_file.write(str(self.max_iter)+" ")

# Write the verbosity level of the code
    output_file.write(str(self.verbose)+"\n")

# Write Fokker-Planck cross section flag
    output_file.write(self.fokker_planck+"\n")

# Write transport correction flag
    output_file.write(self.transport_correction+"\n")

# If transport correction is used, write the optimal transport correction flag
    if self.transport_correction=="true" :
      output_file.write(self.optimal_correction+"\n")

# Write the angular multigrid flag
    output_file.write(self.angular_multigrid+"\n")

# Write the MIP flag
    output_file.write(self.mip+"\n")

# If MIP is used, write the solver used for MIP
    if self.mip=="true" :
      if self.mip_solver=="agmg" :
        output_file.write("AGMG\n")
      elif self.mip_solver=="cg_ml" :
        output_file.write("CG_ML\n")
        if self.aggregation=="uncoupled" :
          output_file.write("Uncoupled\n")
        elif self.aggregation=="mis" :
          output_file.write("MIS\n")
        elif self.aggregation=="uncoupled-mis" :
          output_file.write("Uncoupled_MIS\n")
        else :
          self.Abort("Unknown ML aggregation.")
      elif self.mip_solver=="cg_sgs" :
        output_file.write("CG_SGS\n")
      elif self.mip_solver=="cg_none" :
        output_file.write("CG_None\n")
      else :
        self.Abort("Unknown MIP solver.")

# Write the quadrature type
    if self.quadrature=="ls" :
      output_file.write("LS\n")
    elif self.quadrature=="glc" :
      output_file.write("GLC\n")
    else :
      self.Abort("Unknown quadrature.")

# Write the galerkin flag
    output_file.write(self.galerkin+"\n")

# Write L_max
    output_file.write(str(self.L_max)+"\n")

# Write Sn order
    output_file.write(str(self.sn_order)+"\n")

# Write the FE type
    if self.fe=="bld" :
      output_file.write("BLD\n")
    elif self.fe=="pwld" :
      output_file.write("PWLD\n")
    else :
      self.Abort("Unknown finite element\n")

# Write the source intensity
    for src in self.source :
      output_file.write(str(src)+" ")
    output_file.write("\n")

# Write the bottom incoming flux
    output_file.write(str(self.inc_bottom)+"\n")

# Write the right incoming flux
    output_file.write(str(self.inc_right)+"\n")

# Write the top incoming flux
    output_file.write(str(self.inc_top)+"\n")

# Write the left incoming flux
    output_file.write(str(self.inc_left)+"\n")

# Write the cross sections
    offset = 0
    for i in xrange(self.total_xs.shape[0]) :
      output_file.write(str(self.total_xs[i]))
      for j in xrange(self.scattering_xs.shape[0]/self.total_xs.shape[0]) :
        output_file.write(" "+str(self.scattering_xs[offset+j]))
      output_file.write("\n")
      offset += self.scattering_xs.shape[0]/self.total_xs.shape[0]

    output_file.close()
