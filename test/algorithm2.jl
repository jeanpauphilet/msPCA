n = 5::Int
r = 3::Int
Sigma = rand(Float64, n, n)
Sigma = Sigma*Sigma' #Make it PSD 

ks = [1, 2, 3]::Vector{Int}
numIters = 200::Int
verbose = 1::Int
violation_tolerance = 1e-4::Float64

ofv_best = Ref{Cdouble}(0)
violation_best = Ref{Cdouble}(0)
runtime = Ref{Cdouble}(0)
x_best = zeros(Cdouble, n * r)
ptr_best = pointer(x_best)

@ccall "./algorithm2.so".findmultPCs_deflation(
   ofv_best :: Ptr{Cdouble},
   violation_best :: Ptr{Cdouble},
   runtime :: Ptr{Cdouble},
   pointer(x_best) :: Ptr{Cdouble},
   n :: Cint,
   r :: Cint,
   pointer(reinterpret(Cdouble, reshape(Sigma, n * n))) :: Ptr{Cdouble},
   pointer(reinterpret(Cint, ks)) :: Ptr{Cint},
   numIters :: Cint,
   verbose :: Cint,
   violation_tolerance :: Cdouble
)::Cvoid;

ofv_best = ofv_best[]::Float64
violation_best = violation_best[]::Float64
runtime = runtime[]::Float64
println(x_best)
x_best = reshape(reinterpret(Float64, x_best), n, r)
