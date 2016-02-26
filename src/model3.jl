module Model3
import HADSGE:Distribution,NGrid,ndgrid,buildJ
import HADSGE:Variable,State,Endogenous,Stochastic,Policy,Static,Exogenous,Dependant,Aggregate,Future,AR,Markov,Idiosyncratic,Collective
import HADSGE:addindex!,tchange,subs,subs!,parseex,getexpectation,jacobian,getv,timeof,simplifyindices!

type Model
    G::NGrid
    SP::Array{Float64,2}
    ProbWeights::Array{Float64,2}
    X::Array{Float64,2}
    XP::Array{Float64,2}
    F::Function
    J::Function
    Fval::Array{Float64,2}
    Jval::Array{Float64,2}
    variables::Vector{Variable}
    distribution::Distribution
end
Base.length(M::Model)=length(M.G)

function parsevars(foc::Expr,states::Expr,vars::Expr,params::Expr)
    parameters   = Dict{Symbol,Float64}(zip([x.args[1] for x in params.args],[x.args[2] for x in params.args]))
    variables = Variable[]
    Dlist = Dict{Expr,Expr}()

    for v ∈ vars.args
        if isa(v.args[2],Expr) && v.args[2].args[1]!=:∫ && v.args[2].head !=:tuple
            p = addindex!(v.args[1])=>subs(addindex!(subs(v.args[2],parameters)),Dlist)
            push!(Dlist,p)
            push!(Dlist,tchange(p[1],1)=>tchange(p[2],1))
            push!(variables,Dependant(p[1].args[1],p[2],x->x,[-Inf,Inf]))
        end
    end

    for v ∈ states.args
        if length(v.args[2].args)==3 && reduce(&,map(x->isa(x,Number),v.args[2].args))
            for v1 in vars.args
                if v1.args[1]==v.args[1]
                    push!(variables,Endogenous(v,v1))
                    push!(variables,variables[end].lom)
                end
            end
        elseif length(v.args[2].args)==3
            push!(variables,Markov{Idiosyncratic}(v))
        elseif length(v.args[2].args)==4 && reduce(&,map(x->isa(x,Number),v.args[2].args))
            push!(variables,AR{Idiosyncratic}(v))
        end
    end

    for v in vars.args
        if !in(v.args[1],names(variables)) && !in(parseex(v),[Dependant,Aggregate])
            push!(variables,parseex(v)(v))
        elseif parseex(v)==Aggregate
            if isa(v.args[2].args[2],Symbol)
                if in(v.args[2].args[2],names(State(variables)))
                    push!(variables,Aggregate(v.args[1],State(variables)[v.args[2].args[2]],v.args[2].args[3],[-Inf,Inf]))
                    push!(variables,Aggregate(v.args[1],State(variables,true)[v.args[2].args[2]],v.args[2].args[3],[-Inf,Inf]))
                elseif in(v.args[2].args[2],names(Dependant(variables)))
                    push!(variables,Aggregate(v.args[1],Dependant(variables)[v.args[2].args[2]],v.args[2].args[3],[-Inf,Inf]))
                elseif in(v.args[2].args[2],names(Policy(variables)))
                    push!(variables,Aggregate(v.args[1],Policy(variables)[v.args[2].args[2]],v.args[2].args[3],[-Inf,Inf]))
                else
                    error("Integral target $(v.args[2].args[2]) for variable $(v.args[1]) not found.")
                end
            elseif isa(v.args[2].args[2],Expr)
                p = gensym(v.args[1])=>subs(addindex!(subs(v.args[2].args[2],parameters)),Dlist)
                push!(variables,Dependant(p[1],p[2],x->x,[-Inf,Inf]))
                push!(variables,Aggregate(v.args[1],Dependant(variables)[end],v.args[2].args[3],[-Inf,Inf]))
            end
        end
    end
    variables = vcat(Endogenous(variables),Stochastic(variables),Policy(variables),Exogenous(variables),Aggregate(variables),Dependant(variables))

    f2,j2=parsefoc(foc,variables,Dlist,parameters)

    return variables,parameters,f2,j2
end

function buildfunc(ex::Expr,targ)
    F = Expr(:for,:(i=1:length(M)),Expr(:block))
    for i = 1:length(ex.args)
        push!(F.args[2].args,:($targ[$i,i] = $(ex.args[i])))
    end
    return F
end

function addpweights!(x,nP)
  if typeof(x) == Expr
    if (x.head == :call) && (x.args[1] == :*) && (x.args[2]==:ProbWeights)
      return sum([subs(x.args[3],:j=>j)*:(M.ProbWeights[$j,i]) for j = 1:nP])
    else
        for i = 1:length(x.args)
            x.args[i]=addpweights!(x.args[i],nP)
        end
    end
  end
  return x
end

addpweights(x,nP) = addpweights!(deepcopy(x),nP)

function parsefoc(foc::Expr,variables::Vector{Variable},Dlist::Dict,parameters::Dict)
    f=subs(addindex!(subs(foc,parameters)),Dlist)
    @assert f.head==:vcat || f.head==:vect
    list = :([])
    for i = 1:length(f.args)
        f.args[i],list= getexpectation(f.args[i],list,length(list.args)+1)
    end
    f2 = deepcopy(f)
    subs!(f2,Dict(zip([Expr(:ref,:Expect,i) for i = 1:length(list.args)],list.args)))
    filter!((k,v)->k.args[2]<1,Dlist)
    j2  = jacobian(f2,[Expr(:ref,v,0) for v in names(Policy(variables))])

    for v in sort(setdiff(Symbol[x.args[1] for x in filter(x->x.args[2]==1,getv(f2))],names(Stochastic(variables))))
        push!(variables,Future(v,State(variables,true)[v]))
    end

    hloc = hardloc(variables)
    nP  = reduce(*,map(length,Stochastic(variables)))
    F = simplifyindices!(addpweights!(subs(f2,hloc),nP))
    J = simplifyindices!(vec(addpweights!(subs(j2,hloc),nP)))
    return F,J
end


function hardloc(variables::Vector{Variable})
    ns = length(State(variables))
    hloc = Dict{Expr,Expr}()

    for v in variables
        if isa(v,Endogenous)
            push!(hloc,Expr(:ref,v.name,-1)=>:(M.X[$(findfirst(variables,v)),i]))
        elseif isa(v,Stochastic)
            push!(hloc,Expr(:ref,v.name,0)=>:(M.X[$(findfirst(variables,v)),i]))
            push!(hloc,Expr(:ref,v.name,1)=>:(M.SP[$(findfirst(State(variables),v)),i+(j-1)*length(M)]))
        elseif isa(v,Policy)
            push!(hloc,Expr(:ref,v.name,0)=>:(M.X[$(ns+findfirst(State(variables,true),v)),i]))
        elseif isa(v,Static)
            push!(hloc,Expr(:ref,v.name,timeof(v))=>:(M.X[$(ns+findfirst(State(variables,true),v)),i]))
        elseif isa(v,Future)
            push!(hloc,Expr(:ref,v.name,1)=>:(M.XP[$(findfirst(Future(variables),v)),i+(j-1)*length(M)]))
        end
    end
    return hloc
end

function Model(foc,states,vars,params)
    variables,parameters,F,J=parsevars(foc,states,vars,params)
    G,S,SP,ProbWeights=initS(State(variables))
    X   = [S;zeros(length(Policy(variables))+length(Static(variables)),length(G))]
    XP  = zeros(length(Future(variables)),length(G)*size(ProbWeights,1))

    tG = ndgrid(Vector{Float64}[s.x for s in State(variables)]...)
    distribution = Distribution(tG,zeros(size(tG[1]))+1/length(tG[1]),spzeros(length(tG[1]),length(tG[1])))

    subs!(F,:(length(M))=>length(G))
    F = :($(gensym(:F))(M::Model) = @fastmath $(buildfunc(F,:(M.Fval))))
    subs!(J,:(length(M))=>length(G))
    J = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(buildJ(J)))

    M = Model(G,SP,ProbWeights,X,XP,eval(F),eval(J),zeros(length(Policy(variables)),length(G)),zeros(length(Policy(variables)),length(Policy(variables))),variables,distribution)

    initX(M)
    initD(M)
    return M
end

import SparseGrids:NGrid,QuadraticBF,CC
function initS(slist::Vector{Variable})
    G = NGrid(Int[CC.iM(length(x))-1 for x in slist],hcat(Vector{Float64}[x.bounds for x in slist]...),B=QuadraticBF)
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
    return G,S',SP',ProbWeights'
end

function initX(M::Model)
    for i = 1:length(M.variables)
        v = M.variables[i]
        if isa(v,Policy) || isa(v,Exogenous) || isa(v,Aggregate)
            if isa(v.init,Number)
                M.X[i,:] = v.init
            elseif isa(v.init,Expr)
                init = addindex!(v.init)
                initv = getv(init)
                @assert all(map(x->in(x.args[1],names(State(M.variables))),initv)) "Initialisation of $(v.name) includes a non state variable"
                for j = 1:length(M)
                    initD=[x=>M.X[find(State(M.variables),x.args[1]),j] for x in initv]
                    M.X[i,j]=eval(current_module(),subs(init,initD))
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
                        $(M.X)[$(find(M.variables,v.name)),i] = $(subs(v.e,hloc))
                    end
                end
            end
        v.update=eval(F)
    end
end

getindex(M::Model,v::Symbol) = M.variables[v]
function getindex(M::Model,x::Symbol,t::Int)
    for i ∈ 1:length(M.variables)
        if M.variables[i].name==x && timeof(M.variables[i])==t
            return M.X[i,:]
        end
    end
end



function updateSP(M::Model)
    ns = length(M.G.L)
    for s in Endogenous(M.variables)
        sid = find(State(M.variables),s.name)
        nG     = length(M.G)
        uid = find(State(M.variables,true),s.name)
        for j = 1:size(M.ProbWeights,1)
            for i = 1:nG
                @inbounds M.SP[sid,i+(j-1)*nG] = M.X[ns+uid,i]
            end
        end
    end
end

function updateXP(M::Model)
    fvar = Future(M.variables)
    A = hcat([M[v.target.name,0] for v in fvar]...)
    M.XP[:,:] = M.G(A,M.SP')'
    for i in 1:length(fvar)
        clamp!(M.XP[i,:],fvar[i].target.bounds[1],fvar[i].target.bounds[2])
    end
end

function updateU(M::Model,ϕ=0.9)
    M.F(M)
    n = size(M.Fval,1)
    ns = length(M.G.L)
    for i =1:length(M.G)
        M.J(M,i)
        x = M.Jval\M.Fval[:,i]
        @fastmath @simd for j = 1:n
            @inbounds M.X[ns+j,i] -= ϕ*x[j]
        end
    end
end

function forceboundsU(M::Model)
    ns = length(State(M.variables))
    for v ∈ Policy(M.variables)
        lb,ub = v.bounds[1],v.bounds[2]
        vid=findnext(M.variables,v,ns)
        @inbounds @simd for i = 1:length(M)
            if M.X[vid,i] < lb
                M.X[vid,i] = lb
            elseif M.X[vid,i] > ub
                M.X[vid,i] = ub
            end
        end
    end
end

function solve(M::Model,n::Int=1000,ϕ::Float64=0.8;disp::Int=10000,crit=1e-6)
    MaxError = vec(maximum(abs(M.Fval),2))
    SError = vec(sum(abs(M.Fval).>crit*5,2))/length(M)
    updateallD(M)

    for iter = 1:n
        (mod(iter,20)==0 || iter < 5) && (MaxError = vec(maximum(abs(M.Fval),2));SError = vec(sum(abs(M.Fval).>crit*5,2))/length(M))
        updateSP(M)
        updateXP(M)
        updateU(M,ϕ)
        forceboundsU(M)
        mod(iter,disp)==0 && println(round(log10(MaxError),2),"  ",round(SError,2))
        iter>min(n,5)  && (maximum(MaxError)<crit || all(SError.<0.05)) && (print(iter);break)
    end
    updateallD(M)
end


function updateallD(M::Model)
    for v in Dependant(M.variables)
        v.update(M)
    end
end


end
