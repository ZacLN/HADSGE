
type Distribution
    tG::Tuple{Vararg{Array{Float64}}}
    td::Array{Float64}
    T::SparseMatrixCSC
end
type Model
    G::NGrid
    S::Array{Float64,2}
    SP::Array{Float64,2}
    ProbWeights::Array{Float64,2}
    U::Array{Float64,2}
    X::Array{Float64,2}
    XP::Array{Float64,2}
    F::Function
    J::Function
    Fval::Array{Float64,2}
    Jval::Array{Float64,2}
    variables::Vector{Variable}
    distribution::Distribution
end

display(M::Model) = display(M.variables)
length(M::Model) =  length(M.G)


function Model(foc::Expr,states::Expr,vars::Expr,params::Expr)
    variables,parameters,F,J=parsevars(foc,states,vars,params)

    G,S,SP,ProbWeights=initS(State(variables))
    X   = zeros(length(G),length(Static(variables)))
    U   = zeros(length(G),length(Policy(variables)))
    XP  = zeros(length(G)*size(ProbWeights,2),length(Future(variables)))

    tG = ndgrid(Vector{Float64}[s.x for s in State(variables)]...)
    distribution = Distribution(tG,zeros(size(tG[1]))+1/length(tG[1]),spzeros(length(tG[1]),length(tG[1])))

    F = :( $(gensym(:F))(M::Model) = @fastmath $(buildfunc(F,:(M.Fval))))
    J = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(buildJ(J)))
    subs!(F,:(length(M))=>length(G))
    subs!(J,:(length(M))=>length(G))

    M = Model(G,S,SP,ProbWeights,U,X,XP,eval(F),eval(J),zeros(length(G),length(Policy(variables))),zeros(length(Policy(variables)),length(Policy(variables))),variables,distribution)

    initU(M)
    initX(M)
    initD(M)
    return M
end

(M::Model)(s::Symbol,x::Array{Float64,2}) = M.G(M[s,0],x)

function getindex(M::Model,x::Symbol,t::Int)
    if t==-1
        lst = Endogenous(M.variables)
        for i ∈ 1:length(lst)
            if lst[i].name==x
                return M.S[:,i]
            end
        end
    end
    lst = Policy(M.variables)
    for i ∈ 1:length(lst)
        if lst[i].name==x
            return M.U[:,i]
        end
    end
    lst = Static(M.variables)
    for i ∈ 1:length(lst)
        if lst[i].name==x && timeof(lst[i])==t
            return M.X[:,i]
        end
    end

    lst = State(M.variables)
    for i ∈ 1:length(lst)
        if lst[i].name==x
            return M.S[:,i]
        end
    end
end

getindex(M::Model,v::Symbol) = M.variables[v]


function buildfunc(ex::Expr,targ,t=0)
    F = Expr(:for,:(i=1:length(M.G)),Expr(:block))
    for i = 1:length(ex.args)
        push!(F.args[2].args,:($targ[i,$(i+t)] = $(ex.args[i])))
    end
    return F
end

function buildJ(vJ)
    ex = Expr(:block)
    for i = 1:length(vJ.args)
        push!(ex.args,:(M.Jval[$i] = $(vJ.args[i])))
    end
    push!(ex.args,:(return))
    return ex
end
