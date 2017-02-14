export write_datafile, write_data

function write_datafile( trx::Array{Transmitter,1}, datafile::String, only_loc::Bool )
	
	const ignore = "*"
    const ndataline = 12     # # of data per line: real and imag (Ex Ey Ez Hx Hy Hz)
	
	ntrx = length(trx)
	
	f = open(datafile,"w")
	
	# flag indicating which data should be ignored
	println(f, "IGNORE ", ignore)
	println(f)
	
	# number of transmitters or frequencies
	println(f, "N_TRX ", ntrx)
	println(f)
	
	for itrx = 1:ntrx
		
		# transmitter type
		println(f, "TRX_ORIG")
		
		# number of transmitter points
		npts = size(trx[itrx].trxpts,2)
		println(f, npts)
		
		# transmitter coordinates
		for ipts = 1:npts
			println(f, trx[itrx].trxpts[1,ipts], " ", trx[itrx].trxpts[2,ipts], " ", trx[itrx].trxpts[3,ipts])
		end
		println(f)
		
		# frequency
		println(f, "FREQUENCY ", trx[itrx].omega * 0.5 / pi)
		println(f)
		
		# number of receiver points
		nrcv = size(trx[itrx].rcvpts,2)
		println(f, "N_RECV ", nrcv)
		
		if only_loc
			
			# receiver coordinates
			for ircv = 1:nrcv
				println(f, 
					trx[itrx].rcvpts[1,ircv], " ", trx[itrx].rcvpts[2,ircv], " ", trx[itrx].rcvpts[3,ircv]
					)
			end
			
		else
			
			# receiver coordinates, data and standard deviation
			for ircv = 1:nrcv
				print(f, trx[itrx].rcvpts[1,ircv], " ", trx[itrx].rcvpts[2,ircv], " ", trx[itrx].rcvpts[3,ircv])
				for i = 1:ndataline
					if trx[itrx].data_exist[i,ircv]
						print(f, " ", trx[itrx].data[i,ircv], " ", trx[itrx].sd[i,ircv])
					else
						print(f, " ", ignore, " ", ignore)
					end
				end
				println(f)
			end
		end
		println(f)
	end
	
	close(f)
	
end # write_datafile

function write_data( trx::Array{Transmitter,1}, datafile::String )
	
    const ndataline = 12     # # of data per line: real and imag (Ex Ey Ez Hx Hy Hz)
	
	# number of transmitters and/or frequencies
	ntrx = length(trx)
	
	f = open(datafile,"w")
	
	# header
	@printf(f, "%%%+15s %+15s %+15s", "Xr", "Yr", "Zr")
	@printf(f, " %+13s %+13s %+13s %+13s %+13s %+13s", "", "Ex", "", "Ey", "", "Ez")
	@printf(f, " %+13s %+13s %+13s %+13s %+13s %+13s", "", "Hx", "", "Hy", "", "Hz")
	@printf(f, "\n")
	
	# for each transmitter and/or frequency
	for itrx = 1:ntrx
		
		# frequency
		@printf(f, "\n")
		@printf(f, "%% frequency = %13.6E\n", trx[itrx].omega * 0.5 / pi)
		
		# number of receiver points
		nrcv = size(trx[itrx].rcvpts,2)
		
		# for each receiver
		for ircv = 1:nrcv
			
			# receiver location
			@printf(f, " %15.8E %15.8E %15.8E", trx[itrx].rcvpts[1,ircv], trx[itrx].rcvpts[2,ircv], trx[itrx].rcvpts[3,ircv])
			
			# data
			for idat = 1:ndataline
				@printf(f, " %13.6E", trx[itrx].data[idat,ircv])
			end
			@printf(f, "\n")

		end
		
	end
	
	close(f)
	
end