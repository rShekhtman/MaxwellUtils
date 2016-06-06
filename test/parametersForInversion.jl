# ------- SETUP PARAMETERS FOR THE MODEL AND DATA

# data and topo files
datafile = "data_inv_147.txt"
topofile = "topo.txt"

# of cells in base mesh
n     = vec([ 256  256  512 ])  
# corner of the mesh
x0    = vec([ -536. -540. -1056.0  ])
# total mesh lengths  
meshL = vec([ 1280. 1280. 1536. ]) 

# parameters for meshing
nsmallcells    = vec([3 2 1])  #  # of small cells around each point.
mincellsize    = 1  #  minimum cell size in the data area
depth_core     = vec([60. 50. 100.])  # how far to go down in the core region for the fwd meshes
depth_core_inv = vec([140. 40. 40.])  # how far to go down in the core region for the inv meshes

mincellfactor = 1    # minimum cellsize below topo


# parameters for the forward problem
fname = ""    # leave empty for now
doFV  = true  # use finite volume (other option FEM for finite elements)
doSE  = true  # SE = sensitivity explicit - store sensitivities and not factorizations
 
# reference conductivity 
halfSpaceCond = 1.0e-4

# lower bounds
BL = 1e-4
# Higher bounds
BH = 1.0

# Regularization function
#regfun = wdiffusionReg
# parameters for the regularization function
#regparams = [3.0, 3.0, 1.0, 1e-5]  # alphax  alphay  alphaz  alphas

# For TVp regularization
regfun = wPTV
regparams = [1e-3, 2.0, 2.0, 1.0, 0.1] # epsilon, alphax  alphay  alphaz,  p

#regfun    =  wTVReg
#regparams = [1e-4, 2.0, 2.0, 1.0]

# regularization parameter
beta = 1e-30
# misfit function
misfun = SSDFun
#  inner CG iter
cgit = 10
# maximum iter for the inversion
maxit = 30


# approximate mesh interpolation matrix (inv -> fwd) using [2^ninterp]^3 quadrature points
# (set ninterp = [] to use full interpolation matrix)
ninterp = 3
# ninterp = []

# store global mesh to local mesh interpolation matrix in compact form
compact = true

# model parameter is log conductivity
modfun = expMod

