type brack
    i::Vector{Int}
    w::Float64
end

convert(::Type{Tuple},x::Vector{Int})=ntuple(i->x[i],length(x))
*(a::brack,b::brack) = (x= deepcopy(a);push!(x.i,b.i[1]);x.w=a.w*b.w;x)
kron(a::Array{brack,1})= a
spzeros(n::Int) = SparseVector(n,Int[],Float64[])

function findbracket(v,x::Vector)
    i=searchsortedfirst(x,v)
    if i==1
        return [brack([1],.5);brack([1],.5)]
    else
        dx = x[i]-x[i-1]
        return [brack([i-1],(x[i]-v)/dx);brack([i],(v-x[i-1])/dx)]
    end
end

function updateT(M::Model)
    N = length(M.distribution.tG[1])
    Pf = Vector{Float64}[clamp!(M(s.name,hcat([g[:] for g in M.distribution.tG]...)),s.bounds...) for s in Endogenous(M.variables)]

    dims1  = ntuple(i->length(Endogenous(M.variables)[i]),length(Endogenous(M.variables)))
    pdims1 = prod(dims1)
    Gdims = size(M.distribution.tG[1])

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne = length(Endo)
    for i = 1:N
        brackets = [findbracket(Pf[ie][i],Endo[ie].x) for ie = 1:length(Endo)]
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
    N = length(M.distribution.tG[1])
    Pf = Vector{Float64}[clamp!(M(s.name,hcat([g[:] for g in M.distribution.tG]...)),s.bounds...) for s in Endogenous(M.variables)]
    Gdims = size(M.distribution.tG[1])
    Sdim   = tuple([length(s) for s in Stochastic(M.variables)]...)

    M.distribution.T.colptr = ones(Int,size(M.distribution.T,1)+1)
    M.distribution.T.nzval=Float64[]
    M.distribution.T.rowval=Float64[]

    Endo = Endogenous(M.variables)
    Stoc = Stochastic(M.variables)
    ne = length(Endo)
    for i = 1:N
        brackets = [findbracket(Pf[ie][i],Endo[ie].x) for ie = 1:length(Endo)]
        # S1 = Int[findfirst(s.x,M.distribution.tG[findfirst(names(State(M.variables)),s.name)][i]) for s ∈ Stoc]
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

function updatetd(M::Model)
    d0 = M.distribution.td[:]+.25
    for j = 1:5000
        d1=M.distribution.T*d0
        mean(abs((d1[:]-d0[:])))<1e-10 && break
        d0[:]=d1[:]/sum(d1)
        if j==5000
            warn("Model distribution did not converge")
        end
    end

    M.distribution.td[:] = d0
    return
end

function updatedistribution(M::Model)
    updateT(M)
    updatetd(M)
end





function ∫(M::Model,v::Symbol,t=0)
    isa(M.variables[v],Dependant) && M.variables[v].update(M)
    ∫(M,M[v,t])
end

function ∫(M::Model,v::Vector{Float64})
    ∫(M,reshape(M.G(v,hcat([vec(x) for x in M.distribution.tG]...)),size(M.distribution.tG[1])))
end

function ∫(M::Model,V::Array{Float64})
    d = V.*M.distribution.td
    if all(!Bool[isag(s) for s in State(M.variables)])
        return sum(d)
    else
        ns=length(State(M.variables))
        id = find([isag(s) for s in State(M.variables)])
        D= zeros(size(M.distribution.td,id...))
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

            D[id1...] = sum(d[inds...])/max(sum(M.distribution.td[inds...]),1e-16)
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
            Xid = findfirst(Static(M.variables),s)
            for ii = 1:length(X[1])
                id=BitArray(.*([M[eag[ie].name,t[ie]].==X[ie][ii] for ie = 1:length(eag)]...))
                M.X[id,Xid] *= (1-ϕ)
                M.X[id,Xid] += ϕ*ag[ii]*ones(sum(id))
            end
        end
    else
        for s ∈ Aggregate(M.variables)
            M.X[:,findfirst(Static(M.variables),s)] *= (1-ϕ)
            M.X[:,findfirst(Static(M.variables),s)] += ϕ*∫(M,s.target.name,timeof(s.target))
        end
    end
end


function updateA(M,ϕ=1.0)
    updateallD(M)
    updatedistribution(M)
    updateAggregates(M,ϕ)
end
