export prepareMesh2MeshOT


function prepareMesh2MeshOT(pF::ForwardProbType, Minv::OcTreeMesh, N::Integer, compact::Bool=true)

	P    = getInterpolationMatrix(Minv,pF.Mesh,N)'
	if compact
		Ps = SparseMatrixCSC(P.m,P.n,round(UInt32,P.colptr),round(UInt32,P.rowval),round(Int8,log2(P.nzval)/3))
	else
		Ps = SparseMatrixCSC(P.m,P.n,round(UInt32,P.colptr),round(UInt32,P.rowval),P.nzval)
	end
	return Ps
	
end


function prepareMesh2MeshOT(pFor::RemoteRef{Channel{Any}}, MinvRef::RemoteRef{Channel{Any}}, N::Integer, compact::Bool=true)
	pF   = fetch(pFor)
	Minv = fetch(MinvRef)
	return prepareMesh2MeshOT(pF, Minv, N, compact)
end


function prepareMesh2MeshOT(pFor::Array{RemoteRef{Channel{Any}}},Minv::OcTreeMesh,N::Integer,compact::Bool=true)

	Mesh2Mesh = Array(RemoteRef{Channel{Any}},length(pFor))

	# find out which workers are involved
	workerList = []
	for k=1:length(pFor)
		push!(workerList,pFor[k].where)
	end
	workerList = unique(workerList)
	# send sigma to all workers
	MinvRef = Array(RemoteRef{Channel{Any}},maximum(workers()))
	
	tic()
	@sync begin
		for p=workerList
			@async begin
				MinvRef[p] = remotecall_wait(p,identity,Minv)   # send model to workers
			end
		end
	end
	sendTime=toq()
	println("Time for sending out Minv $sendTime")
		
	@sync begin
		for p=workerList
			@async begin
				for idx=1:length(pFor)
					if p==pFor[idx].where
						Mesh2Mesh[idx] = remotecall_wait(p,prepareMesh2MeshOT,pFor[idx],MinvRef[p],N,compact)
					end
				end
			end
		end
	end
	return Mesh2Mesh
end