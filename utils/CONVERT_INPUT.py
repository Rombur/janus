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

# Read the inner tolerance
    self.Search("INNER TOLERANCE")
    self.inner_tolerance = self.Read_next(1)

# Read the group tolerance
    self.Search("GROUP TOLERANCE")
    self.group_tolerance = self.Read_next(1)

# Read the maximum number of inner iterations
    self.Search("MAX INNER ITER")
    self.max_inner_iter = int(self.Read_next(1))

# Read the maximum number of supergroup iterations
    self.Search("MAX SUPERGROUP ITER")
    self.max_supergroup_iter = int(self.Read_next(1))

# Read the maximum number of group iterations
    self.Search("MAX GROUP ITER")
    self.max_group_iter = int(self.Read_next(1))

# Read the sum of the weights of the quadrature
    self.Search("WEIGHT SUM")
    self.weight_sum = self.data[self.begin:self.end-1].lower()

# Read the type of cross section file
    self.Search("XS TYPE")
    self.xs_type = self.data[self.begin:self.end-1].lower()

# If the cross sections are not Fokker-Planck cross sections read the type of
# permutation used
    if self.xs_type!="fp" :
      self.Search("PERMUTATION TYPE")
      self.permutation_type = self.data[self.begin:self.end-1].lower() 

# Read the verbosity level of the code
    self.Search("VERBOSE")
    self.verbose = int(self.Read_next(1))

# Read if transport correction is used
    self.Search("TRANSPORT CORRECTION")
    self.transport_correction = self.data[self.begin:self.end-1].lower()

# If the transport correction is used, check if the optimal transport
# correction is used
    if self.transport_correction=="true" :
      self.Search("OPTIMAL TRANSPORT CORRECTION")
      self.optimal_correction = self.data[self.begin:self.end-1].lower()

# Read if angular multigrid is used
    self.Search("ANGULAR MULTIGRID")
    self.angular_multigrid = self.data[self.begin:self.end-1].lower()

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
      elif self.mip_solver=="cg_ssor" :
        self.Search("DAMPING FACTOR")
        self.damping_factor = self.Read_next(1)

# Read the quadrature type
    self.Search("QUADRATURE")
    self.quadrature = self.data[self.begin:self.end-1].lower()

# Read if the quadrature is Galerkin
    self.Search("GALERKIN")
    self.galerkin = self.data[self.begin:self.end-1].lower()

# Read Sn order
    self.Search("SN ORDER")
    self.sn_order = int(self.Read_next(1))

# Read the finite elements type
    self.Search("FE")
    self.fe = self.data[self.begin:self.end-1].lower()
# Read the number of groups
    self.Search("NUMBER OF GROUPS")
    self.n_groups = int(self.Read_next(1))

# Read the intensity of the source
    self.Search("SOURCE")
    self.source = self.String_to_np()

# Read the type of boundary of the bottom side
    self.Search("BOTTOM BOUNDARY CONDITION")
    self.bottom_bc_type = self.data[self.begin:self.end-1].lower()

# Read the intensity of the bottom incoming flux if the boundary condition is
# most normal or isotropic incoming flux
    if self.bottom_bc_type=="most normal" or self.bottom_bc_type=="isotropic":
      self.Search("BOTTOM INCOMING FLUX")
      self.inc_bottom = self.Read_next(1)

# Read the type of boundary of the right side
    self.Search("RIGHT BOUNDARY CONDITION")
    self.right_bc_type = self.data[self.begin:self.end-1].lower()

# Read the intensity of the right incoming flux if the boundary condition is
# most normal or isotropic incoming flux
    if self.right_bc_type=="most normal" or self.right_bc_type=="isotropic":
      self.Search("RIGHT INCOMING FLUX")
      self.inc_right = self.Read_next(1)

# Read the type of boundary of the top side
    self.Search("TOP BOUNDARY CONDITION")
    self.top_bc_type = self.data[self.begin:self.end-1].lower()

# Read the intensity of the top incoming flux if the boundary condition is
# most normal or isotropic incoming flux
    if self.top_bc_type=="most normal" or self.top_bc_type=="isotropic":
      self.Search("TOP INCOMING FLUX")
      self.inc_top = self.Read_next(1)

# Read the type of boundary of the left side
    self.Search("LEFT BOUNDARY CONDITION")
    self.left_bc_type = self.data[self.begin:self.end-1].lower()

# Read the intensity of the left incoming flux if the boundary condition is
# most normal or isotropic incoming flux
    if self.left_bc_type=="most normal" or self.left_bc_type=="isotropic":
      self.Search("LEFT INCOMING FLUX")
      self.inc_left = self.Read_next(1)

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

# Write the tolerances
    output_file.write(str(self.inner_tolerance)+" ")
    output_file.write(str(self.group_tolerance)+" ")

# Write the maximum number of iterations
    output_file.write(str(self.max_inner_iter)+" ")
    output_file.write(str(self.max_supergroup_iter)+" ")
    output_file.write(str(self.max_group_iter)+" ")

# Write the sum of the weight of the quadrature
    output_file.write(self.weight_sum+" ")

# Write the type of cross section file
    output_file.write(self.xs_type+" ")
    if self.xs_type!="fp" :
      output_file.write(self.permutation_type+" ")

# Write the verbosity level of the code
    output_file.write(str(self.verbose)+"\n")

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
      elif self.mip_solver=="cg_ssor" :
        output_file.write("CG_SSOR ")
        output_file.write(str(self.damping_factor)+"\n")
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

# Write Sn order
    output_file.write(str(self.sn_order)+"\n")

# Write the FE type
    if self.fe=="bld" :
      output_file.write("BLD\n")
    elif self.fe=="pwld" :
      output_file.write("PWLD\n")
    else :
      self.Abort("Unknown finite element\n")

# Write the number of energy groups
    output_file.write(str(self.n_groups)+"\n")

# Write the source intensity
    for src in self.source :
      output_file.write(str(src)+" ")
    output_file.write("\n")

# Write the bottom boundary condition type
    if self.bottom_bc_type=="vacuum" :
      output_file.write("vacuum\n")
    elif self.bottom_bc_type=="reflective" :
      output_file.write("reflective\n")
    elif self.bottom_bc_type=="most normal" :
      output_file.write("most_normal ")
    elif self.bottom_bc_type=="isotropic" :
      output_file.write("isotropic ")
    else :
      self.Abort("Unknown boundary condition type for the bottom boundary\n")
# Write the bottom incoming flux if necessary
    if self.bottom_bc_type=="most normal" or self.bottom_bc_type=="isotropic" :
      output_file.write(str(self.inc_bottom)+"\n")

# Write the right boundary condition type
    if self.right_bc_type=="vacuum" :
      output_file.write("vacuum\n")
    elif self.right_bc_type=="reflective" :
      output_file.write("reflective\n")
    elif self.right_bc_type=="most normal" :
      output_file.write("most_normal ")
    elif self.right_bc_type=="isotropic" :
      output_file.write("isotropic ")
    else :
      self.Abort("Unknown boundary condition type for the right boundary\n")
# Write the right incoming flux if necessary
    if self.right_bc_type=="most normal" or self.right_bc_type=="isotropic" :
      output_file.write(str(self.inc_right)+"\n")

# Write the top boundary condition type
    if self.top_bc_type=="vacuum" :
      output_file.write("vacuum\n")
    elif self.top_bc_type=="reflective" :
      output_file.write("reflective\n")
    elif self.top_bc_type=="most normal" :
      output_file.write("most_normal ")
    elif self.top_bc_type=="isotropic" :
      output_file.write("isotropic ")
    else :
      self.Abort("Unknown boundary condition type for the top boundary\n")
# Write the top incoming flux if necessary
    if self.top_bc_type=="most normal" or self.top_bc_type=="isotropic" :
      output_file.write(str(self.inc_top)+"\n")

# Write the left boundary condition type
    if self.left_bc_type=="vacuum" :
      output_file.write("vacuum\n")
    elif self.left_bc_type=="reflective" :
      output_file.write("reflective\n")
    elif self.left_bc_type=="most normal" :
      output_file.write("most_normal ")
    elif self.left_bc_type=="isotropic" :
      output_file.write("isotropic ")
    else :
      self.Abort("Unknown boundary condition type for the left boundary\n")
# Write the top incoming flux if necessary
    if self.left_bc_type=="most normal" or self.left_bc_type=="isotropic" :
      output_file.write(str(self.inc_left)+"\n")

    output_file.close()
