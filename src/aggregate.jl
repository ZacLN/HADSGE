type bracket{n}
    i::Vector{Int}
    w::Float64
end

function bracket(v::Float64,x::Vector)
    dx = (x[end]-x[1])/(length(x)-1)
    i,d = divrem(v-x[1],dx)
    if d==0
        return [bracket{1}([Int(i+1)],(1-d/dx))]
    elseif d==1
        return [bracket{1}([Int(i+2)],d/dx)]
    else
        return [bracket{1}([Int(i+1)],(1-d/dx)),bracket{1}([Int(i+2)],d/dx)]
    end
end

promote_type{n1,n2}(::Type{bracket{n1}},::Type{bracket{n2}}) = Type(bracket{n1+n2})
*{n1,n2}(a::bracket{n1},b::bracket{n2}) = bracket{n1+n2}(vcat(a.i,b.i),a.w*b.w)
*{n2}(a::bracket{0},b::bracket{n2}) = b
kron{n}(a::Array{bracket{n},1})= a

convert(::Type{Tuple},x::Vector{Int})=ntuple(i->x[i],length(x))
spzeros(n::Int) = SparseVector(n,Int[],Float64[])

function updateT(M::Model)
    N = length(M.distribution.G[1])
    Pf = Vector{Float64}[clamp!(M(s.name,hcat([g[:] for g in M.distribution.G]...)),s.bounds...) for s in Endogenous(M.variables)]

    dims1  = ntuple(i->length(Endogenous(M.variables)[i]),length(Endogenous(M.variables)))
    pdims1 = prod(dims1)
    Gdims = size(M.distribution.G[1])

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne = length(Endo)
    for i = 1:N
        brackets = [bracket(Pf[ie][i],Endo[ie].x) for ie = 1:length(Endo)]
        S1= collect(ind2sub(Gdims,i))[ne+1:end]
        abrackets = kron(brackets...)

        w = spzeros(pdims1)
        for b in abrackets
            w[sub2ind(dims1,b.i...)]+=b.w
        end
        for S2 = 1:length(Stoc)
            if isag(Stoc[S2])
                w = vcat([w[:]*p for p ∈ Stoc[S2].T[:,S1[S2]]]...)
            else
                w = vcat([w[:]*p for p ∈ Stoc[S2].T[S1[S2],:]]...)
            end
        end

        for j = 1:length(w.nzind)
            M.distribution.T[w.nzind[j],i] += w.nzval[j]
        end
    end
    return
end

function updateT1(M::Model)
    N = length(M.distribution.G[1])
    Pf = Vector{Float64}[clamp!(M(s.name,hcat([g[:] for g in M.distribution.G]...)),s.bounds...) for s in Endogenous(M.variables)]
    Gdims = size(M.distribution.G[1])
    Sdim   = tuple([length(s) for s in Stochastic(M.variables)]...)

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne = length(Endo)
    for i = 1:N
        brackets = [bracket(Pf[ie][i],Endo[ie].x) for ie = 1:length(Endo)]
        S1= collect(ind2sub(Gdims,i))[ne+1:end]
        abrackets = kron(brackets...)

        for j = 1:length(abrackets)
            for jj ∈ CartesianRange(Sdim)
                t = 1.0
                for S2 = 1:length(Stoc)
                    if isag(Stoc[S2])
                        t*= Stoc[S2].T[jj.I[S2],S1[S2]]
                    else
                        t*= Stoc[S2].T[S1[S2],jj.I[S2]]
                    end
                end
                M.distribution.T[sub2ind(Gdims,abrackets[j].i...,jj.I...),i]+=abrackets[j].w*t
            end
        end
    end
    return
end

function updateT2(M::Model)
    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne,ns = length(Endo),length(Stoc)
    n = ne+ns
    Gstoc = hcat([x[:] for x ∈ ndgrid([v.x for v in Stoc]...)]...)
    dims = size(M.distribution.G[1])
    Sdims = ([length(v.x) for v ∈ Stochastic(M.variables)]...)

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    for i ∈ CartesianRange(dims)
        Ei = CartesianIndex{ne}(i.I[1:ne]...)
        Si = CartesianIndex{ns}(i.I[ne+1:end]...)
        X = hcat(Float64[Endo[j].x[i.I[j]] for k=1:size(Gstoc,1),j = 1:ne],Gstoc)
        XP = hcat([clamp!(M(v.name,X),v.bounds...) for v ∈ Endo]...,Gstoc)
        for j ∈ 1:size(Gstoc,1)
            B = bracket(XP[j,1], M.variables[1].x)
            for k = 2:n
                v = M.variables[k]
                B = kron(B,bracket(XP[j,k],v.x))
                if isa(v,Stochastic)
                    [(b.w*=M.variables[k].T[i.I[k],ind2sub(Sdims,j)[k-ne]]) for b in B]
                end
            end
            for b in B
                M.distribution.T[sub2ind(dims,b.i...),sub2ind(dims,i.I...)]+=b.w
            end
        end
    end
    return
end

function updateT3(M::Model)
    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne,ns = length(Endo),length(Stoc)
    n = ne+ns
    Gendo = hcat([x[:] for x ∈ ndgrid([v.x for v in Endo]...)]...)
    Gstoc = hcat([x[:] for x ∈ ndgrid([v.x for v in Stoc]...)]...)

    dims = size(M.distribution.G[1])
    Edims = ([length(v.x) for v ∈ Endo]...)
    Sdims = ([length(v.x) for v ∈ Stoc]...)
    XP =  [reshape(clamp!(M(s.name,hcat([g[:] for g in M.distribution.G]...)),s.bounds...),prod(size(Gendo)),prod(size(Gstoc))) for s in Endo]

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    for i ∈ CartesianRange(dims)
        Ei = CartesianIndex{ne}(i.I[1:ne]...)
        Si = CartesianIndex{ns}(i.I[ne+1:end]...)
        xP = hcat(hcat([XP[j][sub2ind(Edims,Ei.I...),:] for j = 1:ne]...),vec(Gstoc))
        for j ∈ 1:size(Gstoc,1)
            B = [bracket{0}([0],0)]
            for k = 1:n
                v = M.variables[k]
                B = kron(B,bracket(xP[j,k],v.x))
                if isa(v,Stochastic)
                    for b in B
                        b.w*=v.T[i.I[k],ind2sub(Sdims,j)[k-ne]]
                    end
                end
            end
            for b in B
                M.distribution.T[sub2ind(dims,b.i...),sub2ind(dims,i.I...)]+=b.w
            end
        end
    end
    return
end

function updated(M::Model)
    d0 = M.distribution.d[:]+.25
    for j = 1:5000
        d1=M.distribution.T*d0
        mean(abs((d1[:]-d0[:])))<1e-10 && break
        d0[:]=d1[:]/sum(d1)
        if j==5000
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

function ∫(M::Model,v::Symbol,t=0)
    isa(M.variables[v],Dependant) && M.variables[v].update(M)
    ∫(M,M[v,t])
end

function ∫(M::Model,v::Vector{Float64})
    ∫(M,reshape(M.G(v,hcat([vec(x) for x in M.distribution.G]...)),size(M.distribution.G[1])))
end

function ∫(M::Model,V::Array{Float64})
    d = V.*M.distribution.d
    if all(!Bool[isag(s) for s in State(M.variables)])
        return sum(d)
    else
        ns=length(State(M.variables))
        id = find([isag(s) for s in State(M.variables)])
        D= zeros(size(M.distribution.d,id...))
        for i = 1:prod(size(D))
            id1 = ind2sub(size(D),i)
            inds = Any[Colon() for i = 1:ns]
            cnt = 1
            for ii = 1:ns
                if [isag(s) for s in State(M.variables)][ii]
                    inds[ii] = id1[cnt]
                    cnt+=1
                end
            end

            D[id1...] = sum(d[inds...])/max(sum(M.distribution.d[inds...]),1e-16)
        end
        return D
    end
end

function updateAggregates(M::Model,ϕ=1.0)
    if any(Bool[isag(s) for s in State(M.variables)])
        agg=Aggregate(M.variables)
        for s ∈ Aggregate(M.variables)
            ag  = ∫(M,M[s.target.name,timeof(s.target)])
            eag = filter(x->isag(x),State(M.variables))
            t = [isa(e,Endogenous) ? -1 : 0 for e in eag]

            X   = ndgrid([e.x for e in eag]...)
            Xid = findfirst(M.variables,s)
            for ii = 1:length(X[1])
                id=BitArray(.*([M[eag[ie].name,t[ie]].==X[ie][ii] for ie = 1:length(eag)]...))
                M.X[id,Xid] *= (1-ϕ)
                M.X[id,Xid] += ϕ*ag[ii]*ones(sum(id))
            end
        end
    else
        for s ∈ Aggregate(M.variables)
            M.X[:,findfirst(M.variables,s)] *= (1-ϕ)
            M.X[:,findfirst(M.variables,s)] += ϕ*∫(M,s.target.name,timeof(s.target))
        end
    end
end

function setaggregate!(M::Model,s::Symbol,x)
    if all(!Bool[isag(s) for s in State(M.variables)])
        M[s,0]=x
    else
        eag = filter(x->isag(x),State(M.variables))
        t = [isa(e,Endogenous) ? -1 : 0 for e in eag]
        X   = ndgrid([e.x for e in eag]...)
        @assert size(X[1]) == size(x) "size(x) must equal size of aggregate state space"
        Xid = findfirst(M.variables,M[s])
        for ii = 1:length(X[1])
            id=BitArray(.*([M[eag[ie].name,t[ie]].==X[ie][ii] for ie = 1:length(eag)]...))
            M.X[id,Xid] =x[ii]
        end
    end
end

function updateA(M,ϕ=1.0)
    updateallD(M)
    updatedistribution(M)
    updateAggregates(M,ϕ)
end
