# Run with
#
#     julia --project=benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#     julia --project=benchmark scratch/benchmark_rectangular_product_weights.jl
#
# Pass `--codegen` to print the complete LLVM for the three `N = 3` scalar implementations.

using BenchmarkTools
using InteractiveUtils
using Printf
using Statistics

using NSEBase
using ReSolverRectangularGrids: RectangularProductWeights

println("Julia $VERSION on $(Sys.CPU_NAME) ($(Sys.ARCH), $(Threads.nthreads()) thread(s))")

# Original one-line implementation retained here as the benchmark baseline.
struct OriginalProductWeights{N, T, V<:AbstractVector{T}} <: AbstractArray{T, N}
    ws::NTuple{N, V}
end

Base.IndexStyle(::Type{<:OriginalProductWeights}) = IndexCartesian()
Base.size(A::OriginalProductWeights) = map(length, A.ws)
Base.axes(A::OriginalProductWeights) = map(w -> axes(w, 1), A.ws)

Base.getindex(A::OriginalProductWeights{N}, I::Vararg{Int, N}) where {N} =
    prod(ntuple(i -> A.ws[i][I[i]], N))

# Candidate implementation with a separately compiled method for every supported dimensionality.
struct UnrolledProductWeights{N, T, V<:AbstractVector{T}} <: AbstractArray{T, N}
    ws::NTuple{N, V}
end

Base.IndexStyle(::Type{<:UnrolledProductWeights}) = IndexCartesian()
Base.size(A::UnrolledProductWeights) = map(length, A.ws)
Base.axes(A::UnrolledProductWeights) = map(w -> axes(w, 1), A.ws)

Base.@propagate_inbounds function Base.getindex(A::UnrolledProductWeights{1}, i::Int)
    @boundscheck checkbounds(A, i)
    return @inbounds A.ws[1][i]
end

Base.@propagate_inbounds function Base.getindex(A::UnrolledProductWeights{2}, i::Int, j::Int)
    @boundscheck checkbounds(A, i, j)
    return @inbounds A.ws[1][i] * A.ws[2][j]
end

Base.@propagate_inbounds function Base.getindex(
        A::UnrolledProductWeights{3}, i::Int, j::Int, k::Int)
    @boundscheck checkbounds(A, i, j, k)
    return @inbounds A.ws[1][i] * A.ws[2][j] * A.ws[3][k]
end

function cartesian_sum(A)
    value = zero(eltype(A))
    @inbounds for I in CartesianIndices(A)
        value += A[I]
    end
    return value
end

function weighted_energy(A, u)
    value = zero(promote_type(eltype(A), eltype(u)))
    @inbounds for I in CartesianIndices(A)
        value += A[I] * abs2(u[I])
    end
    return value
end

function report(label, original, current, unrolled)
    estimates = median.((original, current, unrolled))
    current_ratio = estimates[2].time / estimates[1].time
    unrolled_ratio = estimates[3].time / estimates[1].time
    @printf("%-26s %10.2f ns  %10.2f ns  %8.3f×  %10.2f ns  %8.3f×  %d/%d/%d allocations\n",
            label, estimates[1].time, estimates[2].time, current_ratio, estimates[3].time,
            unrolled_ratio, estimates[1].allocs, estimates[2].allocs, estimates[3].allocs)
end

function llvm_summary(label, A, signature)
    llvm = sprint(io -> code_llvm(io, getindex, signature; raw=false, debuginfo=:none))
    @printf("%-26s %d floating multiplies; ntuple=%s; prod=%s\n", label,
            length(collect(eachmatch(r"fmul", llvm))), occursin("ntuple", llvm), occursin("prod", llvm))
    return llvm
end

w1 = collect(range(0.5, 1.5; length=257))
w2 = (collect(range(0.5, 1.5; length=65)), collect(range(0.75, 1.75; length=67)))
w3 = (collect(range(0.5, 1.5; length=25)), collect(range(0.75, 1.75; length=27)),
      collect(range(1.0, 2.0; length=29)))

current1 = RectangularProductWeights((w1,))
current2 = RectangularProductWeights(w2)
current3 = RectangularProductWeights(w3)
original1 = OriginalProductWeights((w1,))
original2 = OriginalProductWeights(w2)
original3 = OriginalProductWeights(w3)
unrolled1 = UnrolledProductWeights((w1,))
unrolled2 = UnrolledProductWeights(w2)
unrolled3 = UnrolledProductWeights(w3)

for (original, current, unrolled) in ((original1, current1, unrolled1),
                                      (original2, current2, unrolled2),
                                      (original3, current3, unrolled3))
    @assert original == current == unrolled
    @assert cartesian_sum(original) == cartesian_sum(current) == cartesian_sum(unrolled)
end

u1 = reshape(sin.(range(0, 2π; length=length(current1))), size(current1))
u2 = reshape(sin.(range(0, 2π; length=length(current2))), size(current2))
u3 = reshape(sin.(range(0, 2π; length=length(current3))), size(current3))

dot_size1 = (size(current1)..., 9)
dot_size2 = (size(current2)..., 9)
dot_size3 = (size(current3)..., 9)
dot_u1 = reshape(sin.(range(0, 2π; length=prod(dot_size1))), dot_size1)
dot_u2 = reshape(sin.(range(0, 2π; length=prod(dot_size2))), dot_size2)
dot_u3 = reshape(sin.(range(0, 2π; length=prod(dot_size3))), dot_size3)
dot_v1 = reshape(cos.(range(0, 2π; length=prod(dot_size1))), dot_size1)
dot_v2 = reshape(cos.(range(0, 2π; length=prod(dot_size2))), dot_size2)
dot_v3 = reshape(cos.(range(0, 2π; length=prod(dot_size3))), dot_size3)

i1 = 129
i2, j2 = 33, 34
i3, j3, k3 = 13, 14, 15

scalar_current1 = @benchmark getindex($current1, $i1) seconds=1
scalar_original1 = @benchmark getindex($original1, $i1) seconds=1
scalar_unrolled1 = @benchmark getindex($unrolled1, $i1) seconds=1
scalar_current2 = @benchmark getindex($current2, $i2, $j2) seconds=1
scalar_original2 = @benchmark getindex($original2, $i2, $j2) seconds=1
scalar_unrolled2 = @benchmark getindex($unrolled2, $i2, $j2) seconds=1
scalar_current3 = @benchmark getindex($current3, $i3, $j3, $k3) seconds=1
scalar_original3 = @benchmark getindex($original3, $i3, $j3, $k3) seconds=1
scalar_unrolled3 = @benchmark getindex($unrolled3, $i3, $j3, $k3) seconds=1

sum_current1 = @benchmark cartesian_sum($current1) seconds=1
sum_original1 = @benchmark cartesian_sum($original1) seconds=1
sum_unrolled1 = @benchmark cartesian_sum($unrolled1) seconds=1
sum_current2 = @benchmark cartesian_sum($current2) seconds=1
sum_original2 = @benchmark cartesian_sum($original2) seconds=1
sum_unrolled2 = @benchmark cartesian_sum($unrolled2) seconds=1
sum_current3 = @benchmark cartesian_sum($current3) seconds=1
sum_original3 = @benchmark cartesian_sum($original3) seconds=1
sum_unrolled3 = @benchmark cartesian_sum($unrolled3) seconds=1

energy_current1 = @benchmark weighted_energy($current1, $u1) seconds=1
energy_original1 = @benchmark weighted_energy($original1, $u1) seconds=1
energy_unrolled1 = @benchmark weighted_energy($unrolled1, $u1) seconds=1
energy_current2 = @benchmark weighted_energy($current2, $u2) seconds=1
energy_original2 = @benchmark weighted_energy($original2, $u2) seconds=1
energy_unrolled2 = @benchmark weighted_energy($unrolled2, $u2) seconds=1
energy_current3 = @benchmark weighted_energy($current3, $u3) seconds=1
energy_original3 = @benchmark weighted_energy($original3, $u3) seconds=1
energy_unrolled3 = @benchmark weighted_energy($unrolled3, $u3) seconds=1

dot_current1 = @benchmark NSEBase._dot($dot_u1, $dot_v1, $current1, Val((2,))) seconds=1
dot_original1 = @benchmark NSEBase._dot($dot_u1, $dot_v1, $original1, Val((2,))) seconds=1
dot_unrolled1 = @benchmark NSEBase._dot($dot_u1, $dot_v1, $unrolled1, Val((2,))) seconds=1
dot_current2 = @benchmark NSEBase._dot($dot_u2, $dot_v2, $current2, Val((3,))) seconds=1
dot_original2 = @benchmark NSEBase._dot($dot_u2, $dot_v2, $original2, Val((3,))) seconds=1
dot_unrolled2 = @benchmark NSEBase._dot($dot_u2, $dot_v2, $unrolled2, Val((3,))) seconds=1
dot_current3 = @benchmark NSEBase._dot($dot_u3, $dot_v3, $current3, Val((4,))) seconds=1
dot_original3 = @benchmark NSEBase._dot($dot_u3, $dot_v3, $original3, Val((4,))) seconds=1
dot_unrolled3 = @benchmark NSEBase._dot($dot_u3, $dot_v3, $unrolled3, Val((4,))) seconds=1

println("\nMedian time: original, production, production/original, explicit, explicit/original")
println("-"^118)
report("scalar getindex N=1", scalar_original1, scalar_current1, scalar_unrolled1)
report("scalar getindex N=2", scalar_original2, scalar_current2, scalar_unrolled2)
report("scalar getindex N=3", scalar_original3, scalar_current3, scalar_unrolled3)
report("Cartesian sum N=1", sum_original1, sum_current1, sum_unrolled1)
report("Cartesian sum N=2", sum_original2, sum_current2, sum_unrolled2)
report("Cartesian sum N=3", sum_original3, sum_current3, sum_unrolled3)
report("weighted energy N=1", energy_original1, energy_current1, energy_unrolled1)
report("weighted energy N=2", energy_original2, energy_current2, energy_unrolled2)
report("weighted energy N=3", energy_original3, energy_current3, energy_unrolled3)
report("NSEBase._dot N=1", dot_original1, dot_current1, dot_unrolled1)
report("NSEBase._dot N=2", dot_original2, dot_current2, dot_unrolled2)
report("NSEBase._dot N=3", dot_original3, dot_current3, dot_unrolled3)

println("\nLLVM summaries for scalar getindex")
println("-"^118)
llvm_original1 = llvm_summary("original N=1", original1, Tuple{typeof(original1), Int})
llvm_current1 = llvm_summary("current N=1", current1, Tuple{typeof(current1), Int})
llvm_unrolled1 = llvm_summary("explicit N=1", unrolled1, Tuple{typeof(unrolled1), Int})
llvm_original2 = llvm_summary("original N=2", original2, Tuple{typeof(original2), Int, Int})
llvm_current2 = llvm_summary("current N=2", current2, Tuple{typeof(current2), Int, Int})
llvm_unrolled2 = llvm_summary("explicit N=2", unrolled2, Tuple{typeof(unrolled2), Int, Int})
llvm_original3 = llvm_summary("original N=3", original3,
                              Tuple{typeof(original3), Int, Int, Int})
llvm_current3 = llvm_summary("current N=3", current3, Tuple{typeof(current3), Int, Int, Int})
llvm_unrolled3 = llvm_summary("explicit N=3", unrolled3, Tuple{typeof(unrolled3), Int, Int, Int})

if "--codegen" in ARGS
    println("\nOriginal N=3 LLVM")
    println("-"^118)
    print(llvm_original3)
    println("\nCurrent N=3 LLVM")
    println("-"^118)
    print(llvm_current3)
    println("\nExplicit N=3 LLVM")
    println("-"^118)
    print(llvm_unrolled3)
end
