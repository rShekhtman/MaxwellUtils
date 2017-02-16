export exportOcTreeMeshRoman, exportOcTreeModelRoman, outputOctreeMesh,
       exportOcTreeMeshToVTK,CreateVTKtimeSeriesCollectionOcTree

using jInv.Utils


function exportOcTreeMeshRoman(name::AbstractString, mesh::OcTreeMesh)
# exportOcTree(name::AbstractString,mesh::OcTreeMesh)
# Export OcTree for use with Roman's Fortran codes.

m1,m2,m3     = mesh.n
i1,i2,i3,bsz = find3(mesh.S)
h1,h2,h3     = mesh.h
x1,x2,x3     = mesh.x0

# Roman's code starts the OcTree at the top corner. Change from bottom
# corner.
i3 = m3 + 2 .- i3 - bsz
x3 = x3 + m3 * h3

#S   = sortrows([i3 i2 i1 bsz])
S = sub2ind( (m1,m2,m3), i1,i2,i3 )
p = sortpermFast(S)[1]

n = length(bsz)


# Write OcTree mesh
f = open(name, "w")
println(f, m1, " ", m2, " ", m3, " ! # of cells in underlying mesh")
println(f, x1, " ", x2, " ", x3, " ! top corner")
println(f, h1, " ", h2, " ", h3, " ! cell size")
println(f, n, " ! size of octree mesh")
for i = 1:n
   idx = p[i]
   @printf(f,"%i %i %i %i\n", i1[idx], i2[idx], i3[idx], bsz[idx])
end

close(f)

end

function exportOcTreeModelRoman(name::AbstractString, mesh::OcTreeMesh, u::Array{Float64,1})
# exportModel(name::AbstractString,mesh::OcTreeMesh,u::Array{AbstractFloat,1})
# Export OcTree cell property for use with Roman's Fortran codes.

m1,m2,m3     = mesh.n
i1,i2,i3,bsz = find3(mesh.S)

# Roman's code starts the OcTree at the top corner. Change from bottom
# corner.
i3 = m3 + 2 .- i3 - bsz

n = nnz(mesh.S)

S = sub2ind( (m1,m2,m3), i1,i2,i3 )

p = sortpermFast(S)[1]

v = u[p]

# Write model vector
f = open(name, "w")
for i = 1:n
	println(f, v[i])
end
close(f)

end

function exportOcTreeModelRoman(name::AbstractString, mesh::OcTreeMesh, u::Array{Int64,1})
# exportModel(name::AbstractString,mesh::OcTreeMesh,u::Array{Int,1})
# Export OcTree cell property for use with Roman's Fortran codes.

m1,m2,m3     = mesh.n
i1,i2,i3,bsz = find3(mesh.S)

# Roman's code starts the OcTree at the top corner. Change from bottom
# corner.
i3 = m3 + 2 .- i3 - bsz

n = nnz(mesh.S)

S = sub2ind( (m1,m2,m3), i1,i2,i3 )

p = sortpermFast(S)[1]

v = u[p]

# Write model vector
f = open(name, "w")
for i = 1:n
	println(f, v[i])
end
close(f)

end

#-----------------------------------------------

function outputOctreeMesh(name::AbstractString, mesh::OcTreeMesh)
i1,i2,i3,bsz = find3(mesh.S)
n = nnz(mesh.S)

f = open(name, "w")
for i=1:n
	println(f, i1[i], " ", i2[i], " ", i3[i], " ", bsz[i])
end
close(f)

end

#-----------------------------------------------

function exportOcTreeMeshToVTK(fnameNoExt::AbstractString,mesh::OcTreeMesh,data...)
  #writeOctreeAsVTKunstructured(fnameNoExt,mesh,data...)
  
  #Converts an octree mesh to an unstructured mesh and writes the result to an xml file
  #with compressed binary data storage.
  #Input:
  #  fnameNoExt: file name of output xml file without extension
  #  mesh: jInv octree mesh, can be FEM or FV mesh.
  #  data: Data to be attached to mesh. An arbitrary number of data fields can be attached
  #to the mesh, each as a separate argument. Each field must be written as a tuple in 
  #the format (dataArray,name), where dataArray is a dense array holding the data, it 
  #can be a vector or a multidimensional array. name is a string that gives the name
  #the data field will have in paraview.
  
  #Convert OcTree mesh to vtk unstructured grid format (get node coordinates and 
  #cell connectivity
  pts,cells = convertOctreeToVTKunstructured(mesh)
  
  #Generate xml file with mesh data
  vtk = vtk_grid(fnameNoExt,pts,cells)
  
  #Extract data from varargs and add data fields to xml file
  nnodes = size(pts,2)
  for idat in data
    println(idat[2])
    if (length(idat) == 2)
      addVTKdataOcTree(vtk,idat[1],idat[2],nnodes,mesh)
    else
      addVTKdataOcTree(vtk,idat[1],idat[2],nnodes,mesh,idat[3])
    end
  end
  
  #Write xml file to disk and close it
  vtk_save(vtk)
  println("Wrote mesh to vtk file with data fields:")
  for idat in data
     println(idat[2])
  end
  return pts,cells
end

#-----------------------------------------------

function convertOctreeToVTKunstructured(mesh::OcTreeMesh)
  S = mesh.S
  ncells = mesh.nc
  n1,n2,n3 = mesh.n
  ic,jc,kc,bsz = find3(S)
  inn,jnn,knn,Snnum = getNodalNumberingVTKlocal(S)
  nnodes = nnz(Snnum.SV)
  cells = Array(MeshCell,ncells)
  pts   = zeros(Float32,3,nnodes)
  const celltype = VTKCellType.VTK_HEXAHEDRON
  
  #Get connectivity
  #=  Nodes within a cell numbered as shown (Converted to vtk hexahdron numbering for
      input to WriteVTK functions
     /5----/7
    / |   /|
    1----3 |
    | /6 | /8
    |/   |/
    2----4
  =#
  for cIdx = 1:ncells
    inds    = Array(Int32, 8)
    regNd   = (kc[cIdx]-1)*(n1+1)*(n2+1) + (jc[cIdx]-1)*(n1+1) + ic[cIdx]  #regNd is the
                                                    #underlying regular mesh node number
                                                    #of node 1 in current cell.
    hc      = S.SV.nzval[cIdx]
    inds[1] = Snnum.SV[regNd]
    inds[2] = Snnum.SV[regNd+hc]
    inds[3] = Snnum.SV[regNd+hc*(n1+1)]
    inds[4] = Snnum.SV[regNd+hc*(n1+1)+hc]
    inds[5] = Snnum.SV[regNd+hc*(n1+1)*(n2+1)]
    inds[6] = Snnum.SV[regNd+hc*(n1+1)*(n2+1)+hc]
    inds[7] = Snnum.SV[regNd+hc*(n1+1)*(n2+1) + hc*(n1+1)]
    inds[8] = Snnum.SV[regNd+hc*(n1+1)*(n2+1) + hc*(n1+1)+hc]
    newInd  = [inds[1];inds[2];inds[4];inds[3];inds[5];inds[6];inds[8];inds[7]]
    cells[cIdx] = MeshCell(celltype,newInd)
  end
  
  #Get node coordinates
  pts = getNodalGrid(mesh)'
  
  return pts, cells
end

#-----------------------------------------------

function getNodalNumberingVTKlocal(S::SparseArray3D)
  #Same as getNodalNumbering function in JOcTree/getNodalNumbering.jl except
  #that this function returns 3D indices of nodes.
  
  m1,m2,m3 = S.sz
  i,j,k,bsz = find3(S); bsz = round(Int64,bsz)

  nind = [i        j       k;
          i        j       k+bsz;
          i        j+bsz   k;
          i        j+bsz   k+bsz;
          i+bsz    j       k;
          i+bsz    j       k+bsz;
          i+bsz    j+bsz   k;
          i+bsz    j+bsz   k+bsz ]

  Ni = sparse3(nind[:,1],nind[:,2],nind[:,3],nind[:,3],[m1+1,m2+1,m3+1])
  i,j,k = find3(Ni)
  N = sparse3(i,j,k,1:length(i), [m1+1,m2+1,m3+1]);
  return i,j,k,N
end
