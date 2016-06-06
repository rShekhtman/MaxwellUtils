export rm_cmnts, read_receiver_data, read_receiver, read_transmitter, read_datafile

#----------------------------------------------------------

function skipcmnts(f::IOStream)
# Skip lines that start with '!' or '#' or are empty.
   while true
      if eof(f)
         return ""  # end of file reached
      end

      line = readline(f)
      line = split(line)
      emptyline = length(line) == 0
      commented = !emptyline && (startswith(line[1],"!") || startswith(line[1],"#"))
      if !emptyline && !commented
         return line
      end
   end # while true
end  # function skipcmnts

#----------------------------------------------------------

function rm_cmnts( s::ASCIIString )
# Remove comments (everything after '!')
   f = searchindex(s,"!")
   if  f > 0 
      s = s[1:f-1]
   end
   return s
end  # rm_cmnts

#----------------------------------------------------------

function read_receiver_data( f::IOStream, nrcv::Int64,
                             ignore::SubString{ASCIIString} )
# Read receiver locations, data and standard deviations.

   const ndataline = 12     # # of data per line: real and imag (Ex Ey Ez Hx Hy Hz)
   
   rcvpts = zeros(3,nrcv)
   
   data = zeros(ndataline,nrcv)
   sd   = zeros(ndataline,nrcv)   # standard deviation
   data_exist = falses(ndataline,nrcv)
   ndata = 0   # counter for the number of data
   
   for i = 1:nrcv
      line = skipcmnts(f)
      for j = 1:3
         rcvpts[j,i] = parse(Float64,line[j])  # receiver coordinates
      end
      k = 4   # column counter
      
      for j = 1:ndataline
         if (line[k] != ignore)
            data[j,i] = parse(Float64,line[k])
              sd[j,i] = parse(Float64,line[k+1])
            data_exist[j,i] = true
            ndata = ndata + 1
         end
         k = k + 2
      end  # j
      
   end  # i
   
   return rcvpts, data, sd, data_exist, ndata
end  # read_receiver_data

#----------------------------------------------------------

function read_receiver( f::IOStream, nrcv::Int64,
                        ignore::SubString{ASCIIString} )
# Read only receiver locations (no data).

   rcvpts = Array(Float64, 3,nrcv)
   
   for i = 1:nrcv
      line = skipcmnts(f)
      for j = 1:3
         rcvpts[j,i] = parse(Float64,line[j])  # receiver coordinates
      end
   end  # i
   
   return rcvpts
end  # read_receiver

#----------------------------------------------------------

function read_transmitter( f::IOStream )
# Read transmitter and frequency information.

   line = skipcmnts(f)  # transmitter type
   trxtp = line[1]

   if trxtp == "TRX_LOOP" || trxtp == "TRX_LOOP_X" || trxtp == "TRX_LOOP_Y" 

      line = skipcmnts(f)  # loop centre (x,y,z), radius
      xc = parse(Float64,line[1])
      yc = parse(Float64,line[2])
      zc = parse(Float64,line[3])
      r  = parse(Float64,line[4])
      const npt = 17  # number of loop points

      if trxtp == "TRX_LOOP"
         x,y = makeloop( xc, yc, r, npt )
         z = fill(zc, npt)
      elseif trxtp == "TRX_LOOP_X"
         y,z = makeloop( yc, zc, r, npt )
         x = fill(xc, npt)
      elseif trxtp == "TRX_LOOP_Y"
         x,z = makeloop( xc, zc, r, npt )
         y = fill(yc, npt)
      end

      trxpts = [ x y z ]'

   else  # read line segments
      line = skipcmnts(f)
      npt = parse(Int,line[1])  # # of points making up the transmitter
      
      trxpts = Array(Float64, 3,npt)
      
      # Read the transmitter coordinates
      for i = 1:npt
         line = skipcmnts(f)
         for j = 1:3
            trxpts[j,i] = parse(Float64,line[j])
         end
      end  # i
   
   end


   line = skipcmnts(f)  # FREQUENCY f
   frequency = parse(Float64,line[2])
   omega = 2*pi*frequency

   return trxpts, omega
end  # read_transmitter

#----------------------------------------------------------

function makeloop( xc, yc, r, n )
# Make a loop around (xc,yc), radius r,
# with n points, first point==last point

const d = 2*pi/(n-1)

x = Array(Float64, n)
y = Array(Float64, n)

for i = 1:n-1
	angle = d*(i-1)
   x[i] = xc  +  r * cos(angle)
   y[i] = yc  +  r * sin(angle)
end # i

x[n] = x[1]
y[n] = y[1]

return x, y
end # function makeloop

#----------------------------------------------------------

function read_datafile( datafile::ASCIIString,  # data filename
                        only_loc::Bool )  # true to only read locations, false for data and sd

   f = open(datafile,"r")
   
   line = skipcmnts(f)  # IGNORE ii
   ignore = line[2]  # flag indicating which data should be ignored
   
   line = skipcmnts(f)  # N_TRX n
   ntrx = parse(Int,line[2])  # # of transmitters or frequencies
   
   trx = Array{Transmitter}(ntrx)
   
   for itrx = 1:ntrx
      
      (trxpts, omega) = read_transmitter(f)
      
      line = skipcmnts(f)  # N_RECV n
      nrcv = parse(Int,line[2])  # # of receivers
      
      if only_loc
         rcvpts = read_receiver( f, nrcv, ignore )
         trx[itrx] = Transmitter( trxpts, omega, rcvpts,
                                  0,0,0,0 )
      else
         ( rcvpts, data, sd, data_exist, ndata) = read_receiver_data( f, nrcv, ignore )
         trx[itrx] = Transmitter( trxpts, omega, rcvpts,
                                  data, sd, data_exist, ndata )
      end
   
   end  # itrx
   
   close(f)
   
   return trx
   
end  # read_datafile
