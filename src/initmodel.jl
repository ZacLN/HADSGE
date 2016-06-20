import SparseGrids:NGrid,Quadratic,Linear
function initS(slist::Vector{Variable},B=Quadratic,G1 = nothing)
    if G1==nothing
        G = NGrid(Int[SparseGrids.cc_iM(length(x))-1 for x in slist],hcat(Vector{Float64}[x.bounds for x in slist]...),B=B)
    else
        G = G1
    end
    S = values(G)
    for s in slist
        s.x = sort(unique(S[:,find(slist,s.name)]))
    end

    exogenousfuture=map(x->repmat(x[:]',length(G),1),ndgrid([s.x for s in Stochastic(slist)]...))
    nP=prod(map(length,Stochastic(slist)))
    ProbWeights =  Float64[prod([Stochastic(slist)[N].T[findfirst(Stochastic(slist)[N].x,S[i,findfirst(names(State(slist)),names(Stochastic(slist))[N])]),findfirst(Stochastic(slist)[N].x,exogenousfuture[N][i,j])] for N = 1:length(Stochastic(slist))]) for i = 1:length(G),j=1:nP]
    SP = zeros(length(G)*nP,size(S,2))
    for i = 1:length(Stochastic(slist))
        SP[:,findfirst(names(slist),Stochastic(slist)[i].name)] = exogenousfuture[i][:]
    end
    return G,S,SP,ProbWeights
end




function initX(M::Model)
    for i = 1:length(M.variables)
        v = M.variables[i]
        if isa(v,Policy) || isa(v,Exogenous) || isa(v,Aggregate)
            if isa(v.init,Number)
                M.X[:,i] = v.init
            elseif isa(v.init,Expr)
                init = addindex!(v.init)
                initv = getv(init)
                @assert all(map(x->in(x.args[1],names(State(M.variables))),initv)) "Initialisation of $(v.name) includes a non state variable"
                for j = 1:length(M)
                    initD=[x=>M.X[j,find(State(M.variables),x.args[1])] for x in initv]
                    M.X[j,i]=eval(current_module(),subs(init,initD))
                end
            end
        end
    end
end


function initD(M::Model)
    hloc = hardloc(M.variables)
    for v in Dependant(M.variables)
        F = quote
                function $(gensym(v.name))(M)
                    for i = 1:length(M)
                        $(M.X)[i,$(find(M.variables,v.name))] = $(subs(v.e,hloc))
                    end
                end
            end
        v.update=eval(current_module(),F)
    end
end
