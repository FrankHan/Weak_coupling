#!/usr/bin/env python2.7
# -*- coding: UTF-8 no BOM -*-

#--------------------------------------------------------------------------------------------------
# This code will be redundant once DADF5 with addCalculations exists
#--------------------------------------------------------------------------------------------------
import h5py, re
import os
import numpy as np
import argparse
import math


# ------------------------------------------------------------------
class DADF5():
  """Read and write to DADF5 files"""
  
# ------------------------------------------------------------------
  def __init__(self,
               filename,
               geom_file,
               mode     = 'r',
              ):
    
    if mode not in ['a','r']:
      print('Invalid file access mode')
      with h5py.File(filename,mode):
        pass
      
    with h5py.File(filename,mode) as f:
      r=re.compile('inc[0-9]+')
      self.increments = [{'group':  u,
                          'inc':    int(u.split('inc')[1]),
                          'time':   f[u].attrs['time/s'],
                          'active': True
                          }   for u in f.keys() if r.match(u)]
    
    # ToDo: Fortran we need to ensure that all sections (not only active ones) are written out
    # the current way of skipping inactive phases is too complicated
    
      self.output_types = {}
      
      self.output_types['constituent'] = \
          {'active': True,
           'sections':[
                       {'group':  u,
                        'label':  u.split('_',1)[1],
                        'active': True,
                        }  for u in f['inc00000/constituent']
                      ]
          }
           

      self.output_types['constitutive'] = \
          {'active': True,
           'sections':[
                       {'group':  u,
                        'label':  u.split('_',1)[1],
                        'active': True,
                        }  for u in f['inc00000/constitutive']
                      ]
          }
          
     # self.nodes = f['/geometry/node_0'][()]
     # self.connectivity = f['/geometry/connectivity_cell'][()]-1
     # 
     # self.n_nodes = np.shape(self.nodes)[0]
     # self.n_cells = np.shape(self.connectivity)[0]
     # self.n_nodes_per_cell = np.shape(self.connectivity)[1]

    self.filename = filename
    self.geom_filename = geom_file
    self.mode     = mode
    self.hdf_file = h5py.File(filename)
    
  def select_section(self,section):
    return True

  def select_output_type(self,output_type):
    return True
    
  def __update(self):
    return True
    
  def HDF5_dataset_obj(self,current,existing = None):    #returns an object list which has information for the datasets
    if existing is None:
      existing = []
    for g in current.values():
      if isinstance(g,h5py.Group):
        self.HDF5_dataset_obj(g,existing)
      else:
        existing.append(g)
    return existing
    
  def get_active_groups(self):
    groups = []
    for inc in self.increments:
      if not inc['active']: continue
      for out in self.output_types.keys():
        if not self.output_types[out]['active']: continue
        for i,sec in enumerate(self.output_types[out]['sections']):
          if not self.output_types[out]['sections'][i]['active']: continue
          groups.append(inc['group']+'/'+out+'/'+sec['group'])
    return groups

    
  def get_candidates(self,hdf_file,label,groups = None):  #gives path list of the dataset given as argument
    if groups is None:
      groups = []
    if type(label) is not list:
      print('mist')
    for d in self.HDF5_dataset_obj(hdf_file):
      if os.path.basename(d.name) != label[0]: continue
      try:
        for l in label:
          d.parent[l]
        groups.append(d.parent.name)
      except:
        print('Quantities not in the same group',d.parent)
    return groups
    
  def get_data(self,l,inc):
    i = '{:05d}'.format(inc)
    
    groups = {}
    if type(l) is not list:
      print('mist')
    with h5py.File(self.filename,'r') as f:
      for g in self.get_active_groups():
        if set(l).issubset(f[g].keys()) and g[0:8]== 'inc'+i:
          for j in l:
            groups[j] = f[g+'/'+j][()]
    return groups

  def addCalculation(self,label_arg,formula_arg):  #label,formula similar to how it is done currently in string form eg: 'np.sum(#rho_dip#;#rho_edge#)'
    formula_arg = formula_arg.replace(';',',')
    calc_quants = re.findall(r'#(.+?)#',formula_arg)
    with h5py.File(self.filename,'a') as f:
      for groups in self.get_candidates(f,calc_quants):
        for count,args in enumerate(calc_quants):
          if count == 0:
            work_array = np.array(f[groups + '/' + calc_quants[count]])
          else:
            work_array = np.append(work_array,f[groups + '/' + calc_quants[count]],axis = 1)
        formula_arg = re.sub(r'\((.+?)\)','(work_array,axis=1)',formula_arg)  #this would work only for sum/product
        output = eval(formula_arg)
        f[groups + '/' + label_arg] = output





#class DADF5_displacement():
#
#  def __init__(self,
#               filename,
#               geom_file,
#               mode     = 'r',
#              ):
#    if mode not in ['a','r']:
#      print('Invalid file access mode')
#      with h5py.File(filename,mode):
#        pass
#      
#    with h5py.File(filename,'r') as f:
#      r=re.compile('inc[0-9]+')
#      self.increments = [{'group':  u,
#                          'inc':    int(u.split('inc')[1]),
#                          'time':   f[u].attrs['time/s'],
#                          'active': True
#                          }   for u in f.keys() if r.match(u)]
#    self.filename = filename
#    self.geom_filename = geom_file

  def geom_info(self):
    with open(self.geom_filename,'r') as f:
      lines = f.readlines()
    
    g = re.compile(r'^(grid)')
    s = re.compile(r'^(size)')
    grid = []
    size = []
    for i in lines:
      if re.search(g,i):
        grid.extend((i.split()[2],i.split()[4],i.split()[6]))
      if re.search(s,i):
        size.extend((i.split()[2],i.split()[4],i.split()[6]))
    
    grid = np.array([int(i) for i in grid])
    size = np.array([float(i) for i in size])

    delta = size/grid*0.5

    x, y, z = np.meshgrid(np.linspace(delta[2],size[2]-delta[2],grid[2]),
                          np.linspace(delta[1],size[1]-delta[1],grid[1]),
                          np.linspace(delta[0],size[0]-delta[0],grid[0]),
                          indexing = 'ij')
                          
    coords = np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3)
    coords = coords.reshape([np.product(grid),3])
    self.coords = coords
    self.grid   = grid
    self.size   = size
    return self.coords,grid,size
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
  def cell2node(cellData,grid):

    nodeData = 0.0
    datalen = np.array(cellData.shape[3:]).prod()
    
    for i in range(datalen):
      node = scipy.ndimage.convolve(cellData.reshape(tuple(grid[::-1])+(datalen,))[...,i],
                                    np.ones((2,2,2))/8.,                                              # 2x2x2 neighborhood of cells
                                    mode = 'wrap',
                                    origin = -1,                                                      # offset to have cell origin as center
                                   )                                                                  # now averaged at cell origins
      node = np.append(node,node[np.newaxis,0,:,:,...],axis=0)                                        # wrap along z
      node = np.append(node,node[:,0,np.newaxis,:,...],axis=1)                                        # wrap along y
      node = np.append(node,node[:,:,0,np.newaxis,...],axis=2)                                        # wrap along x

      nodeData = node[...,np.newaxis] if i==0 else np.concatenate((nodeData,node[...,np.newaxis]),axis=-1)

    return nodeData

  #--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
  def displacementAvgFFT(self,F,grid,size,nodal=False,transformed=False):
    """Calculate average cell center (or nodal) displacement for deformation gradient field specified in each grid cell"""
    if nodal:
      x, y, z = np.meshgrid(np.linspace(0,size[2],1+grid[2]),
                            np.linspace(0,size[1],1+grid[1]),
                            np.linspace(0,size[0],1+grid[0]),
                            indexing = 'ij')
    else:
      x, y, z = np.meshgrid(np.linspace(0,size[2],grid[2],endpoint=False),
                            np.linspace(0,size[1],grid[1],endpoint=False),
                            np.linspace(0,size[0],grid[0],endpoint=False),
                            indexing = 'ij')

    origCoords = np.concatenate((z[:,:,:,None],y[:,:,:,None],x[:,:,:,None]),axis = 3) 

    F_fourier = F if transformed else np.fft.rfftn(F,axes=(0,1,2))                                    # transform or use provided data
    Favg = np.real(F_fourier[0,0,0,:,:])/grid.prod()                                                  # take zero freq for average
    avgDisplacement = np.einsum('ml,ijkl->ijkm',Favg-np.eye(3),origCoords)                            # dX = Favg.X

    return avgDisplacement

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
  def displacementFluctFFT(self,F,grid,size,nodal=False,transformed=False):
    """Calculate cell center (or nodal) displacement for deformation gradient field specified in each grid cell"""
    integrator = 0.5j * size / math.pi

    kk, kj, ki = np.meshgrid(np.where(np.arange(self.grid[2])>self.grid[2]//2,np.arange(self.grid[2])-self.grid[2],np.arange(self.grid[2])),
                             np.where(np.arange(self.grid[1])>self.grid[1]//2,np.arange(self.grid[1])-self.grid[1],np.arange(self.grid[1])),
                                      np.arange(self.grid[0]//2+1),
                             indexing = 'ij')
    k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3) 
    k_sSquared = np.einsum('...l,...l',k_s,k_s)
    k_sSquared[0,0,0] = 1.0                                                                           # ignore global average frequency

  #--------------------------------------------------------------------------------------------------
  # integration in Fourier space

    displacement_fourier = -np.einsum('ijkml,ijkl,l->ijkm',
                                      F if transformed else np.fft.rfftn(F,axes=(0,1,2)),
                                      k_s,
                                      integrator,
                                     ) / k_sSquared[...,np.newaxis]

  #--------------------------------------------------------------------------------------------------
  # backtransformation to real space

    displacement = np.fft.irfftn(displacement_fourier,grid[::-1],axes=(0,1,2))

    return cell2node(displacement,grid) if nodal else displacement


  def displacement_data_proc(self):
      self.geom_info()
      for inc in self.get_candidates(self.hdf_file,['F']):
        F_fourier = np.fft.rfftn(np.array(self.hdf_file[inc + '/F']).reshape(self.grid[2],self.grid[1],self.grid[0],3,3),axes=(0,1,2))      # perform transform only once...

        fluctDisplacement = self.displacementFluctFFT(F_fourier,self.grid,self.size,transformed=True)
        avgDisplacement   = self.displacementAvgFFT  (F_fourier,self.grid,self.size,transformed=True)  #options.nodal has been made False
  # ------------------------------------------ output data -------------------------------------------
##adding the displacement values to the constituent group at each increment
##The data is per material point. Therefore, for geometry of 64x64x64, you get displacement array of 64x64x64

        self.hdf_file[inc + '/avg_x'] = avgDisplacement.reshape([np.product(self.grid),3])
        self.hdf_file[inc + '/fluct_x'] = fluctDisplacement.reshape([np.product(self.grid),3])
        self.hdf_file[inc + '/new_pos'] = self.coords + np.array(self.hdf_file[inc + '/avg_x']) + np.array(self.hdf_file[inc + '/fluct_x'])
        self.hdf_file[inc + '/orig_pos'] = self.coords
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

  def gradFFT(self,geomdim,field):
    """Calculate gradient of a vector or scalar field by transforming into Fourier space."""
    self.geom_info()
    shapeFFT = np.array(np.shape(field))[0:3]
    grid     = np.array(np.shape(field)[2::-1])
    N = self.grid.prod()                                                                                    # field size
    n = np.array(np.shape(field)[3:]).prod()                                                           # data size
  
    field_fourier = np.fft.rfftn(field,axes=(0,1,2),s=shapeFFT)
    grad_fourier  = np.empty(field_fourier.shape+(3,),'c16')
  
    # differentiation in Fourier space
    TWOPIIMG = 2.0j*math.pi
    einsums = { 
                1:'ijkl,ijkm->ijkm',                                                                   # scalar, 1 -> 3
                3:'ijkl,ijkm->ijklm',                                                                  # vector, 3 -> 3x3
              }
  
    k_sk = np.where(np.arange(self.grid[2])>self.grid[2]//2,np.arange(self.grid[2])-self.grid[2],np.arange(self.grid[2]))/geomdim[0]
    if self.grid[2]%2 == 0: k_sk[self.grid[2]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)
  
    k_sj = np.where(np.arange(self.grid[1])>self.grid[1]//2,np.arange(self.grid[1])-self.grid[1],np.arange(self.grid[1]))/geomdim[1]
    if self.grid[1]%2 == 0: k_sj[self.grid[1]//2] = 0                                                            # Nyquist freq=0 for even grid (Johnson, MIT, 2011)
  
    k_si = np.arange(self.grid[0]//2+1)/geomdim[2]
  
    kk, kj, ki = np.meshgrid(k_sk,k_sj,k_si,indexing = 'ij')
    k_s = np.concatenate((ki[:,:,:,None],kj[:,:,:,None],kk[:,:,:,None]),axis = 3).astype('c16')                           
    grad_fourier = np.einsum(einsums[n],field_fourier,k_s)*TWOPIIMG
  
    return np.fft.irfftn(grad_fourier,axes=(0,1,2),s=shapeFFT).reshape([N,3*n])

#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
  def addGradient(self,label_arg): #label_arg should be in list form #only works for scalar datasets
      self.geom_info()
      for inc in self.get_candidates(self.hdf_file,label_arg): 
          grad_FFT = self.gradFFT(self.size[::-1],
                                  np.array(self.hdf_file[inc + '/' + label_arg[0]]).reshape(self.grid[::-1].tolist() + [1]))
          self.hdf_file[inc + '/grad_FFT'] = grad_FFT

  def addNorm(self,label_arg): 
      for inc in self.get_candidates(self.hdf_file,label_arg):
          all_vectors = np.array(self.hdf_file[inc + '/' + label_arg[0]])
          row = 0
          norm_array = np.zeros(len(all_vectors))
          for vector in all_vectors:
              norm_array[row] = np.linalg.norm(vector)
              row = row + 1
          self.hdf_file[inc + '/normFrobenius_{}'.format(label_arg[0])] = norm_array    
              

