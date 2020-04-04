using VirulenceEvolution
using VirulenceEvolution:getindex, setindex!
using Test

# test ind2sub
t = reshape(1:8, 2, 2, 2)
@test t[VirulenceEvolution.ind2sub(axes(t), 6)...] == 6

# test Dynamics
d = Dynamics(0, 1)
@test d.history == [(0, 1)]
record!(d, 1, 2)
@test d.history == [(0, 1), (1, 2)]

# test zeroatdiag
m = [1 2
     3 4]
VirulenceEvolution.zeroatdiag!(m)
@test m[1, 1] == m[2, 2] == 0

# test @apply
x = 1
y = 2
VirulenceEvolution.@apply x->x+1 x y
@test x == 2
@test y == 3

# test slicecol
@test m[:, 2] == VirulenceEvolution.slicecol(m, 2)

# test getindex by tuple
@test m[(1, 1)] == m[1, 1]

# test setindex by tuple
m[(1, 2)] = 5.0
@test m[1, 2] == 5

# test gillespie work
f(x::T) where {T} = rand(T)

gillespie_single(f, 50; S=1000, I=[1], R=[0], v=[0.3],
    β=0.01, b=0.2, d=0.1, r=0.1, μ=0.1, m=0.01, maxepoch=30_000)

gillespie_meta(f, 50; S=[1000, 100], I=[1 0], R=[0 0], v=[0.3],
    β=0.01, b=0.2, d=0.1, r=0.1, μ=0.1, mt=0.01, mg=0.01, maxepoch=50_000)
