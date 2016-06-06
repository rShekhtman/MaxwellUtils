
export getSxRxFromData

function getSxRxFromData(trx::Transmitter, recLen::Array{Float64,1} = ones(3))

#ntrx = length(trx)  # total # of transmitters or frequencies

# Setup sources
Srcs = Array(Array{Float64}, 1)
#for i = 1:ntrx
Srcs[1] = trx.trxpts'
#end # i


# Count total # of receivers
#nrcv_total = 0
#for i = 1:ntrx
nrcv_total = size(trx.rcvpts, 2) 
#end # i
nrcv_total *= 6   # for Ex,Ey,Ez, Hx,Hy,Hz


# setup dipole receivers for E and/or square loops for H
dx = recLen[1] / 2
dy = recLen[2] / 2
dz = recLen[3] / 2
Recs = Array(Array{Float64}, nrcv_total)
k = 1  #  receiver counter

#for i = 1:ntrx
	rcv = trx.rcvpts
	data_exist = trx.data_exist
   doAllRcv = data_exist==0  # for only fwd, do all components

   nrcv = size(rcv, 2)

   for j = 1:nrcv
   	if doAllRcv || data_exist[1,j] || data_exist[2,j]
	      Recs[k] = [ rcv[1,j]-dx rcv[2,j] rcv[3,j] ;
	                  rcv[1,j]+dx rcv[2,j] rcv[3,j] ]  # Ex
	      k += 1              
	   end

   	if doAllRcv || data_exist[3,j] || data_exist[4,j]
	      Recs[k] = [ rcv[1,j] rcv[2,j]-dy rcv[3,j] ;
	                  rcv[1,j] rcv[2,j]+dy rcv[3,j] ]  # Ey
	      k += 1         
      end   

   	if doAllRcv || data_exist[5,j] || data_exist[6,j]
	      Recs[k] = [ rcv[1,j] rcv[2,j] rcv[3,j]-dz ;
	                  rcv[1,j] rcv[2,j] rcv[3,j]+dz ]  # Ez
	      k += 1              
      end   
                    

   	if doAllRcv || data_exist[7,j] || data_exist[8,j]
	      Recs[k] = [ rcv[1,j] rcv[2,j]-dy rcv[3,j]-dz ;
	                  rcv[1,j] rcv[2,j]-dy rcv[3,j]+dz ;
	                  rcv[1,j] rcv[2,j]+dy rcv[3,j]+dz ;
	                  rcv[1,j] rcv[2,j]+dy rcv[3,j]-dz ;
					  rcv[1,j] rcv[2,j]-dy rcv[3,j]-dz  ;]  # Hx
	      k += 1              
      end   
            
   	if doAllRcv || data_exist[9,j] || data_exist[10,j]
	      Recs[k] = [ rcv[1,j]-dx rcv[2,j] rcv[3,j]-dz ;
	                  rcv[1,j]-dx rcv[2,j] rcv[3,j]+dz ;
	                  rcv[1,j]+dx rcv[2,j] rcv[3,j]+dz ;
	                  rcv[1,j]+dx rcv[2,j] rcv[3,j]-dz ;
					  rcv[1,j]-dx rcv[2,j] rcv[3,j]-dz ;]  # Hy
	      k += 1              
      end   

   	if doAllRcv || data_exist[11,j] || data_exist[12,j]
	      Recs[k] = [ rcv[1,j]-dx rcv[2,j]-dy rcv[3,j] ;
	                  rcv[1,j]-dx rcv[2,j]+dy rcv[3,j] ;
	                  rcv[1,j]+dx rcv[2,j]+dy rcv[3,j] ;
	                  rcv[1,j]+dx rcv[2,j]-dy rcv[3,j] ;
					  rcv[1,j]-dx rcv[2,j]-dy rcv[3,j] ;]  # Hz
	      k += 1              
      end   
            
                    
   end # j
#end # i

if k < nrcv_total
   deleteat!(Recs, k:nrcv_total)
end


return Srcs, Recs

end  # getSurvey
