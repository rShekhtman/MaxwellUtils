
export getDobsWdFromTrx

function getDobsWdFromTrx( trx::Transmitter )

nrcv_total = size(trx.data, 2) 
nrcv_total *= 6   # for Ex,Ey,Ez, Hx,Hy,Hz

Dobs = Array(Complex128, nrcv_total)
Wd   = Array(Complex128, nrcv_total)

nrcv = size(trx.data, 2)
k = 1  #  receiver counter

for j = 1:nrcv

   for i = 1:2:12  # components
   	if trx.data_exist[i,j] || trx.data_exist[i+1,j]

         if trx.data_exist[i,j]  # real component
            d_r  = trx.data[i,j]
            Wd_r = 1.0 / trx.sd[i,j]
         else
            d_r  = 0.0  # dummy value
            Wd_r = 0.0
         end
         
         if trx.data_exist[i+1,j]  # imaginary component
            d_i  = trx.data[i+1,j]
            Wd_i = 1.0 / trx.sd[i+1,j]
         else
            d_i  = 0.0  # dummy value
            Wd_i = 0.0
         end

         Dobs[k] = complex(d_r, d_i)
         Wd[k] = complex(Wd_r, Wd_i)
	      k += 1              
         
      end
   end  # i

end  # j


deleteat!(Dobs, k:nrcv_total)
deleteat!(Wd,   k:nrcv_total)

return Dobs, Wd
end # function getDobsWdFromTrx	