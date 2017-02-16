module MaxwellUtils

using jInv.Utils
using jInv.Mesh
using jInv.ForwardShare
hasJOcTree = false
try
  using JOcTree
  hasJOcTree = true
catch
  hasJOcTree = false
end

# ----------------------------------------------------------
type Transmitter
   trxpts      # (3,npts) points making up the transmitter
   omega       # 2*pi*frequency
   rcvpts      # (3,nrcv) receiver locations
   # The following is used only when only_loc=false
   data        # (12,nrcv)
   sd          # (12,nrcv) standard deviation
   data_exist  # (12,nrcv) true if there is data
   ndata       # number of data
end # type Transmitter
export Transmitter


# The 12 data columns are ordered:
#  1   2    3   4    5   6    7   8    9   10   11  12
# Exr Exi  Eyr Eyi  Ezr Ezi  Hxr Hxi  Hyr Hyi  Hzr Hzi

# ---------------------------------------------------------

include("readDataFiles.jl")
include("getTrxOmega.jl")
include("getSxRxFromData.jl")
include("getInitialModel.jl")
include("getDobsWdFromTrx.jl")
include("arrayUtils.jl")
include("diagM.jl")
include("QuickHull.jl")

if hasJOcTree
  include("readUBCData.jl")
  include("writeUBCData.jl")
  include("setupMeshParam.jl")
  include("readTopo.jl")
  include("createSmallMeshFromTX.jl")
  include("exportMesh.jl")
  include("importMesh.jl")
  include("setupBigOctreeMeshPolygon.jl")
  include("prepareMesh2MeshOT.jl") 
  include("createOcTreeFromTRX.jl")
end

end
