using jInv.Utils
using MaxwellUtils

@testset "Data reading tests" begin
datafile = "data_locations.txt"
only_loc = true   # set to true to read only locations for fwd modeling (no data)

println("Reading transmitter and receiver locations")
tic()
trx = read_datafile(datafile, only_loc )
toc()

# mesh file (output)
meshFile = "octree.txt"

# cell size
cellSize = (25.0, 25.0, 25.0)

# minimum distance survey to boundary (padding)
padding = 1000.0 # 1 km
#
# extend fine cells below survey to depth
depthFine = 500.0 # 100 m

# create OcTree mesh
println("Creating OcTree mesh from transmitter and receiver locations")
tic()
M = getOcTreeFromTRX(trx,cellSize,padding,depthFine,numFineLayers=3,numCoarseLayers=2,outFile=meshFile)
toc()

end #End of data reading test set