
using maxwellUtils
@everywhere begin
	using Mesh
	using Utils
	using MaxwellFrequency
	using LinearSolvers
	using InverseSolve
		
end


# ------- SETUP PARAMETERS FOR THE MODEL AND DATA

# data and topo files
datafile = "data_prob_bounds.txt"
topofile = "topo.txt"

# of cells in base mesh
n     = vec([ 128  128  256 ])
# corner of the mesh
x0    = vec([ 475111. 7326905. -2585.  ])
# total mesh lengths  
meshL = vec([ 2560. 2560. 3840. ])




# parameters for meshing
nsmallcells    = vec([3 2 1])  #  # of small cells around each point.
mincellsize    = 1  #  minimum cell size in the data area
depth_core     = vec([60. 100. 300.])  # how far to go down in the core region for the fwd meshes
depth_core_inv = vec([300. 100. 100.])  # how far to go down in the core region for the inv meshes

mincellfactor = 1    # minimum cellsize below topo


# parameters for the forward problem
fname = ""    # leave empty for now
doFV  = true  # use finite volume (other option FEM for finite elements)
doSE  = true  # SE = sensitivity explicit - store sensitivities and not factorizations
 
# reference conductivity 
halfSpaceCond = 1.0e-4

# lower bounds
BL = 1e-5
# Higher bounds
BH = 1.0

# Regularization function
regfun = wdiffusionReg
regparams = [3.0, 3.0, 1.0, 1e-5]  # alphax  alphay  alphaz  alphas


## For TVp regularization
#regfun = wPTV
#regparams = [1e-3, 2.0, 2.0, 1.0, 0.1] # epsilon, alphax  alphay  alphaz,  p



# regularization parameter
beta = 1e-32
# misfit function
misfun = SSDFun
#  inner CG iter
cgit = 15
# maximum iter for the inversion
maxit = 30

 
