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


function prepareMesh2MeshOT(pFor::RemoteChannel, MinvRef::Future, N::Integer, compact::Bool=true)
	pF   = fetch(pFor)
	Minv = fetch(MinvRef)
	return prepareMesh2MeshOT(pF, Minv, N, compact)
end


function prepareMesh2MeshOT(pFor::Array{RemoteChannel},Minv::OcTreeMesh,N::Integer,compact::Bool=true)

	Mesh2Mesh = Array(RemoteChannel,length(pFor))

	# find out which workers are involved
	workerList = []
	for k=1:length(pFor)
		push!(workerList,pFor[k].where)
	end
	workerList = unique(workerList)
	# send sigma to all workers
	MinvRef = Array(Future,maximum(workers()))
	
	tic()
	@sync begin
		for p=workerList
			@async begin
				MinvRef[p] = remotecall_wait(identity,p,Minv)   # send model to workers
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
						Mesh2Mesh[idx] = remotecall_wait(prepareMesh2MeshOT,p,pFor[idx],MinvRef[p],N,compact)
					end
				end
			end
		end
	end
	return Mesh2Mesh
end