using Blink,Base.Threads


type Distribution
    G::Tuple{Vararg{Array{Float64}}}
    d::Array{Float64}
    T::SparseMatrixCSC
end
type ModelSummary
    equations::Expr
    parameters::Dict
    displays::Vector{Window}
end

type Model
    G::NGrid
    SP::Array{Float64,2}
    ProbWeights::Array{Float64,2}
    X::Array{Float64,2}
    XP::Array{Float64,2}
    F::Function
    Fall::Function
    Fval::Array{Float64,2}
    variables::Vector{Variable}
    distribution::Distribution
    temp
    summary::ModelSummary
end

display(M::Model) = display(M.variables)
length(M::Model) =  length(M.G)


function Model(foc,states,vars,params,B=Quadratic)
    foc1 = deepcopy(foc)
    variables,parameters,F,J = parsevars(foc,states,vars,params)
    G,S,SP,ProbWeights       = initS(State(variables),B)
    F1 = deepcopy(F)
    Fall                     = buildSolver(F,J,G,variables)
    X       = [S zeros(length(G),length(Policy(variables))+length(Static(variables)))]
    XP      = zeros(length(G)*size(ProbWeights,2),length(Future(variables)))

    tG      = ndgrid(Vector{Float64}[s.x for s in State(variables)]...)
    distribution = Distribution(tG,zeros(size(tG[1]))+1/length(tG[1]),spzeros(length(tG[1]),length(tG[1])))

    F1       = :($(gensym(:F))(M::Model) = @fastmath $(buildfunc(F1,:(M.Fval))))

    M = Model(G,SP,ProbWeights,X,XP,eval(F1),eval(Fall),zeros(length(G),length(Policy(variables))),variables,distribution,(Fall,),ModelSummary(foc1,parameters,[]))
    precompile(M.Fall,(Model,))

    initX(M)
    initD(M)

    solve((M),2)

    return M
end

(M::Model)(s::Symbol,x::Array{Float64,2}) = M.G(M[s,0],x)

function getindex(M::Model,x::Symbol,t::Int)
    for i ∈ 1:length(M.variables)
        if M.variables[i].name==x && timeof(M.variables[i])==t
            return M.X[:,i]
        end
    end
end

function setindex!(M::Model,x,v::Symbol,t::Int)
    for i ∈ 1:length(M.variables)
        if M.variables[i].name==v && timeof(M.variables[i])==t
            M.X[:,i] = x
        end
    end
end

getindex(M::Model,v::Symbol) = M.variables[v]

function buildfunc(ex::Expr,targ,t=0)
    F = :(for i = 1:length(M) end)
    for i = 1:length(ex.args)
        push!(F.args[2].args,:($targ[i,$(i+t)] = $(ex.args[i])))
    end
    return F
end


@inline pow(a,b) = exp(b*log(a))
function subspow!(x)
    if isa(x,Expr)
        if x.head==:call && x.args[1]==:^
            x.args[1] = :pow
            for i = 2:length(x.args)
                x.args[i] = subspow!(x.args[i])
            end
        else
            for i = 1:length(x.args)
                x.args[i] = subspow!(x.args[i])
            end
        end
    end
    return x
end
# function to update policy variables
# 1 - defineds equation errors
# 2 - finds jacobian
# 3 - 1 Newton step (x = x-ϕ*J(x)\F(x))
# 4 - strip ~long repeated expressions

function buildSolver(F,J,G,variables)
    Fexpr=Expr(:function,:($(gensym("updateX"))(M::Model,ϕ=0.8)),Expr(:block,:(Fval=zeros($(length(Policy(variables))))),:(Jval=zeros($(length(Policy(variables))),$(length(Policy(variables)))))))
    preloopblock = Expr(:block)
    mainloopblock = Expr(:block)
    for i = 1:length(F.args)
        push!(mainloopblock.args,:(Fval[$i]=$(F.args[i])))
    end

    for i = 1:length(J.args)
        if isa(J.args[i],Number)
            J.args[i] !=0 && push!(Fexpr.args[2].args,:(Jval[$i]=$(J.args[i])))
        else
            push!(mainloopblock.args,:(Jval[$i]=$(J.args[i])))
        end
    end

    push!(mainloopblock.args,:(A_ldiv_B!(lufact(Jval),Fval)))
    push!(mainloopblock.args,Expr(:for,:(j=1:$(length(Policy(variables)))),Expr(:macrocall,Symbol("@inbounds"),:(M.X[i,$(length(State(variables)))+j]-=ϕ*Fval[j]))))
    mainloop = Expr(:macrocall,Symbol("@threadsfixed"),:([Fval,Jval]),Expr(:for,:(i=1:(length(M))),mainloopblock))
    # mainloop = Expr(:for,:(i=1:(length(M))),mainloopblock)
    push!(Fexpr.args[2].args,mainloop)
    # Fexpr = subspow!(Fexpr)
    return Fexpr
end


# subs!(mainloopblock,:(length(M))=>length(G))

# list = findrepeated(mainloopblock);
# v = [collect(keys(list)) collect(values(list)) map(length,collect(keys(list)))];
# id = (v[:,2].>5) & (v[:,3].>20)
# v   = v[id,:]
# v = v[reverse(sortperm(v[:,3])),:]
#
# for v in v[:,1]
#     gsn = Symbol("tempvar"*randstring(3))
#     subs!(mainloopblock,v=>gsn)
#     unshift!(mainloopblock.args,Expr(:(=),gsn,v))
# end
