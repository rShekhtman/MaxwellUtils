export spdiags, sdiag, fwdTriSolve, bwdTriSolve



function sdiag(a::Vector)
# S = sdiag(s) builds sparse diagonal matrix
	n = length(a)
	i = collect(1:n+1) # colptr
	j = collect(1:n)   # rowval
	return SparseMatrixCSC(n,n,i,j,a)
end


function spdiags(B,d,m,n)
#   spdiags(B,d,m,n)
# creates a sparse matrix from its diagonals
d = d[:]
p = length(d)

len = zeros(p+1,1)
for k = 1:p
    len[k+1] = round(Int64,len[k]+length(max(1,1-d[k]):min(m,n-d[k])))
end
a = zeros(round(Int64,len[p+1]),3)
for k = 1:p
    # Append new d[k]-th diagonal to compact form
    i = max(1,1-d[k]):min(m,n-d[k])
    a[(round(Int64,len[k])+1):round(Int64,len[k+1]),:] = [i i+d[k] B[i+(m>=n)*d[k],k]]
end

A = sparse(round(Int64,a[:,1]),round(Int64,a[:,2]),a[:,3],m,n);

return A

end


function fwdTriSolve(A::SparseMatrixCSC,b::Vector{Float64})
# forward substitution for CSC matrices
n = length(b)
ncol, nrow = size(A)
if ncol!=nrow; error("A must be square"); end;
if n!=ncol; error("A must be nXn"); end;

x = zeros(n)
if typeof(b[1]) == Complex{Float64}; x = complex( x); end

aa = A.nzval
ja = A.rowval
ia = A.colptr

for j=1:n
   i1 = ia[j]
   i2 = ia[j+1]-1
   for i = i1:i2
      if ja[i] == j
         x[j] = (b[j]-x[j])/aa[i]
      elseif ja[i] > j
         x[ja[i]] += aa[i]*x[j]
      end
   end
end
return x
end


function bwdTriSolve(A::SparseMatrixCSC,b::Vector{Float64})
# backward substitution for CSC matrices
n = length(b)
ncol, nrow = size(A)
if ncol!=nrow; error("A must be square"); end;
if n!=ncol; error("A must be nXn"); end;

x = zeros(n)
if typeof(b[1]) == Complex{Float64}; x = complex(x); end


aa = A.nzval
ja = A.rowval
ia = A.colptr

for j=n:-1:1
   i1 = ia[j]
   i2 = ia[j+1]-1
   for i = i2:-1:i1
       if ja[i]==j
          x[j] = (b[j]-x[j])/aa[i]
       elseif ja[i]<j
          x[ja[i]] += aa[i]*x[j]
       end
   end
end
return x
end
