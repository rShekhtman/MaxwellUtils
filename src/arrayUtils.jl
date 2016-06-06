export ComplexHadamard

function ComplexHadamard(x::Array,y::Array)	
	return real(x).*real(y) + im*(imag(x).*imag(y))
end