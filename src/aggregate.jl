type bracket{n}
    i::Vector{Int}
    w::Float64
end

function bracket(v::Float64, x::Vector)
    dx = (x[end] - x[1]) / (length(x) - 1)
    i,d = divrem(v - x[1], dx)
    if d == 0
        return [bracket{1}([Int(i + 1)], (1 - d / dx))]
    elseif d == 1
        return [bracket{1}([Int(i + 2)], d / dx)]
    else
        return [bracket{1}([Int(i + 1)], (1 - d / dx)),bracket{1}([Int(i + 2)], d / dx)]
    end
end

promote_type{n1, n2}(::Type{bracket{n1}}, ::Type{bracket{n2}}) = Type(bracket{n1 + n2})
*{n1, n2}(a::bracket{n1}, b::bracket{n2}) = bracket{n1 + n2}(vcat(a.i, b.i), a.w * b.w)
*{n2}(a::bracket{0}, b::bracket{n2}) = b
kron{n}(a::Array{bracket{n}, 1}) = a

convert(::Type{Tuple}, x::Vector{Int}) = ntuple(i -> x[i], length(x))
spzeros(n::Int) = SparseVector(n, Int[], Float64[])



"
   updateT(M)

Compute transition matrix for model M using tensor grid.
"
function updateT(M::Model)
    N = length(M.distribution.G[1])
    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    Pf = Vector{Float64}[clamp!(M(s.name, hcat([g[:] for g in M.distribution.G]...)), s.bounds...) for s in Endogenous(M.variables)]
    Gdims = size(M.distribution.G[1])
    Sdim   = tuple([length(s) for s in Stochastic(M.variables)]...)

    ns = length(Stoc)
    nT = prod(Sdim) * length(Endo) * 2
    vals = zeros(N, nT)
    rloc = zeros(UInt32, N, nT)
    cloc = zeros(UInt32, N, nT)
    ccnt = zeros(UInt32, N)

    ne = length(Endo)

    @threads for i = 1:N
        brackets = [bracket(Pf[ie][i], Endo[ie].x) for ie = 1:length(Endo)]
        S1 = collect(ind2sub(Gdims, i))[ne + 1:end]
        abrackets = kron(brackets...)

        cnt = 0
        for j = 1:length(abrackets)
            for jj ∈ CartesianRange(Sdim)
                t = 1.0
                for S2 = 1:ns
                    if isag(Stoc[S2])
                        t *= Stoc[S2].T[jj.I[S2], S1[S2]]
                    else
                        t *= Stoc[S2].T[S1[S2], jj.I[S2]]
                    end
                end
                cnt += 1
                rloc[i,cnt] = sub2ind(Gdims, abrackets[j].i..., jj.I...)
                vals[i,cnt] = abrackets[j].w * t
                cloc[i,cnt] = i
                ccnt[i] = cnt
            end
        end
    end
    rloc, cloc, vals = vec(rloc), vec(cloc), vec(vals)
    id = (vals .== 0)
    M.distribution.T = sparse(rloc[!id], cloc[!id], vals[!id], N, N)
    return
end


function updated(M::Model)
    d0 = M.distribution.d[:] + .25
    for j = 1:5000
        d1 = M.distribution.T * d0
        mean(abs.((d1[:] - d0[:]))) < 1e-10 && break
        d0[:] = d1[:] / sum(d1)
        if j == 5000
            warn("Model distribution did not converge")
        end
    end

    M.distribution.d[:] = d0
    return
end

function updatedistribution(M::Model)
    updateT(M)
    updated(M)
end

function ∫(M::Model, v::Symbol, t = 0)
    isa(M.variables[v], Dependant) && M.variables[v].update(M)
    ∫(M, M[v, t])
end

function ∫(M::Model, v::Vector{Float64})
    ∫(M, reshape(M.G(v, hcat([vec(x) for x in M.distribution.G]...)), size(M.distribution.G[1])))
end

function ∫(M::Model, V::Array{Float64})
    d = V .* M.distribution.d
    if all(!Bool[isag(s) for s in State(M.variables)])
        return sum(d)
    else
        ns = length(State(M.variables))
        id = find([isag(s) for s in State(M.variables)])
        D = zeros(size(M.distribution.d, id...))
        for i = 1:prod(size(D))
            id1 = ind2sub(size(D), i)
            inds = Any[Colon() for i = 1:ns]
            cnt = 1
            for ii = 1:ns
                if [isag(s) for s in State(M.variables)][ii]
                    inds[ii] = id1[cnt]
                    cnt += 1
                end
            end

            D[id1...] = sum(d[inds...]) / max(sum(M.distribution.d[inds...]), 1e-16)
        end
        return D
    end
end

function updateAggregates(M::Model, ϕ = 1.0)
    if any(Bool[isag(s) for s in State(M.variables)])
        agg = Aggregate(M.variables)
        for s ∈ Aggregate(M.variables)
            ag  = ∫(M, M[s.target.name, timeof(s.target)])
            eag = filter(x -> isag(x), State(M.variables))
            t = [isa(e, Endogenous) ? -1 : 0 for e in eag]

            X   = ndgrid([e.x for e in eag]...)
            Xid = findfirst(M.variables, s)
            # tX = zeros(size(M.X,1))
            for ii = 1:length(X[1])
                id = BitArray(.*([M[eag[ie].name,t[ie]] .== X[ie][ii] for ie = 1:length(eag)]...))
                M.X[id,Xid] *= (1 - ϕ)
                M.X[id,Xid] += ϕ * ag[ii] * ones(sum(id))
            end
        end
    else
        for s ∈ Aggregate(M.variables)
            # M.X[:,findfirst(M.variables,s)] *= (1-ϕ)
            # M.X[:,findfirst(M.variables,s)] += ϕ*∫(M,s.target.name,timeof(s.target))
            M.X[:,findfirst(M.variables, s)] = (1 - ϕ) * M.X[:, findfirst(M.variables, s)] + ϕ * ∫(M, s.target.name, timeof(s.target))
        end
    end
end

function setaggregate!(M::Model, s::Symbol, x)
    if all(!Bool[isag(s) for s in State(M.variables)])
        M[s,0] = x
    else
        eag = filter(x -> isag(x), State(M.variables))
        t = [isa(e, Endogenous) ? -1 : 0 for e in eag]
        X   = ndgrid([e.x for e in eag]...)
        @assert size(X[1]) == size(x) "size(x) must equal size of aggregate state space"
        # Xid = findfirst(M.variables,M[s])
        Xid = findfirst(x -> x.name == s && !isa(x, State), M.variables)
        for ii = 1:length(X[1])
            id = BitArray(.*([M[eag[ie].name,t[ie]] .== X[ie][ii] for ie = 1:length(eag)]...))
            M.X[id,Xid] = x[ii]
        end
    end
end

function updateA(M, ϕ = 1.0)
    updateallD(M)
    updatedistribution(M)
    updateAggregates(M, ϕ)
end
