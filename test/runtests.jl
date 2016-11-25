using MaxwellUtils
using jInv.Mesh
using jInv.Utils
using MaxwellFrequency
using jInv.LinearSolvers

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

hasJOcTree = false
try 
  using JOcTree
  hasJOcTree = true
catch
  warn("JOcTree unavailable, no testing will be performed")
end

if hasJOcTree
include("testReadData.jl")

@testset "Miscellaneous tests" begin
#sub2ind = sub2indFast
#ind2sub = ind2subFast

datafile = ["testFileForE3DFormat.txt"]
topofile = "topo.txt"

n     = vec([ 128  64  64 ])  # # of cells
x0    = vec([ -64. -64. -80  ])  # corner
meshL = vec([ 200. 128. 128. ])   # mesh lengths

nsmallcells = vec([3 2 1])  #  # of small cells around each point.
mincellsize = 1  #  minimum cell size in the data area
depth_core = vec([4. 8. 20.])  # how far to go down in the core region for the fwd meshes
mincellfactor = 1    # minimum cellsize below topo


trx, h, itopo = setupMeshParam( datafile, topofile, n,x0,meshL )

meshingParam = cell(5)
meshingParam[1] = nsmallcells 
meshingParam[2] = mincellsize
meshingParam[3] = itopo
meshingParam[4] = depth_core
meshingParam[5] = mincellfactor

M = createSmallMeshFromTX(trx[1], vec(h),vec(n),vec(x0), 
                     nsmallcells, mincellsize, itopo, depth_core, mincellfactor )
println("Mesh size: ", nnz(M.S))

#exportOcTreeMeshRoman("mesh.txt",M)

Srcs, Recs = getSxRxFromData(trx[1])


nEx,nEy,nEz = getEdgeNumbering(M.S)

# transmitters
s   = zeros(Complex128,sum(M.ne),length(Srcs))
for k=1:length(Srcs)
		s[:,k] = sparse(getEdgeIntegralOfPolygonalChain(M,Srcs[k],nEx,nEy,nEz))
end

   # receivers
Obs = spzeros(sum(M.ne),length(Recs))
for k=1:length(Recs)
		Obs[:,k] = sparse(getEdgeIntegralOfPolygonalChain(M,Recs[k],nEx,nEy,nEz))
end

freq    = trx[1].omega
Sources = s

end #End test set
end #End hasJOcTree if block