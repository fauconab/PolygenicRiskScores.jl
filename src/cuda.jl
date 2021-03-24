import .CUDA
import .CUDA: CuArray

CUDA.allowscalar(false)

cuarray() = CuArray

backend(::Type{CuArray}) = KernelAbstractions.CUDADevice()
mzeros(::Type{CuArray{T,N}}, args...) where {T,N} = CUDA.zeros(T, args...)
mones(::Type{CuArray{T,N}}, args...) where {T,N} = CUDA.ones(T, args...)
mrand(::Type{CuArray{T,N}}, args...) where {T,N} = CUDA.rand(T, args...)
mrandn(::Type{CuArray{T,N}}, args...) where {T,N} = CUDA.randn(T, args...)
dot(A::Transpose{T,CuArray{T,N}}, B::CuArray{T,N}) where {T,N} = A*B
