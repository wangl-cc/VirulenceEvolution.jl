using VirulenceEvolution

f(x::T) where {T} = rand(T)

gillespie_single(f, 50; S=1000, I=[1], R=[0], v=[0.3],
    β=0.01, b=0.2, d=0.1, r=0.1, μ=0.1, m=0.01, maxepoch=30_000)

gillespie_meta(f, 50; S=[1000, 100], I=[1 0], R=[0 0], v=[0.3],
    β=0.01, b=0.2, d=0.1, r=0.1, μ=0.1, mt=0.01, mg=0.01, maxepoch=30_000)
