
type Distribution
    G::Tuple{Vararg{Array{Float64}}}
    d::Array{Float64}
    T::SparseMatrixCSC
end
type Model
    G::NGrid
    SP::Array{Float64,2}
    ProbWeights::Array{Float64,2}
    X::Array{Float64,2}
    XP::Array{Float64,2}
    F::Function
    J::Function
    Fall::Function
    Fval::Array{Float64,2}
    Jval::Array{Float64,2}
    variables::Vector{Variable}
    distribution::Distribution
    temp
end

display(M::Model) = display(M.variables)
length(M::Model) =  length(M.G)


function Model(foc,states,vars,params)
    variables,parameters,F,J=parsevars(foc,states,vars,params)
    G,S,SP,ProbWeights=initS(State(variables))
    X   = [S zeros(length(G),length(Policy(variables))) zeros(length(G),length(Static(variables)))]
    XP  = zeros(length(G)*size(ProbWeights,2),length(Future(variables)))

    tG = ndgrid(Vector{Float64}[s.x for s in State(variables)]...)
    distribution = Distribution(tG,zeros(size(tG[1]))+1/length(tG[1]),spzeros(length(tG[1]),length(tG[1])))

    F = :($(gensym(:F))(M::Model) = @fastmath $(buildfunc(F,:(M.Fval))))
    subs!(F,:(length(M))=>length(G))
    Jp = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(buildJ2(J)))
    J = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(buildJ(J)))
    subs!(J,:(length(M))=>length(G))

    M = Model(G,SP,ProbWeights,X,XP,eval(F),eval(J),x->x,zeros(length(G),length(Policy(variables))),zeros(length(Policy(variables)),length(Policy(variables))),variables,distribution,(F,J,eval(Jp)))


    M.Fall = buildFall(M)
    initX(M)
    initD(M)
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
    # F = Expr(:for,:(i=1:length(M)),Expr(:block))
    F = :(for i = 1:length(M) end)
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

function buildJ2(vJ)
    n=div(length(vJ.args),2)
    ex = Expr(:block)
    push!(ex.args,:(Jval = zeros($n,$n)))
    for i = 1:length(vJ.args)
        push!(ex.args,:(Jval[$i] = $(vJ.args[i])))
    end
    push!(ex.args,:(return Jval))
    return ex
end





function buildFall(M)
    N,npf = length(M),size(M.ProbWeights,2)
    F = deepcopy(M.temp[1])
    HADSGE.vecind!(F,N,npf)
    HADSGE.simplifyindices!(F)
    F.args[1] = Expr(:call,gensym(:F),:Fval,:Jval,:X,:XP,:SP,:ProbWeights,:i)
    F.args[2].args = F.args[2].args[2].args[2].args[2].args
    for i = 1:length(F.args[2].args)
        F.args[2].args[i] = Expr(:macrocall,symbol("@inbounds"),F.args[2].args[i])
    end

    J = deepcopy(M.temp[2])
    for i = 1:size(M.Fval,2)^2
        HADSGE.subs!(J,:(M.Jval[$i])=>:(Jval[$i+(i-1)*$(size(M.Fval,2)^2)]))
    end
    HADSGE.vecind!(J,N,npf)
    J = J.args[2].args[2].args[2].args[1:end-1]
    for j in J
        push!(F.args[2].args,j)
    end

    push!(F.args[2].args,:(return))
    F3 = eval(F)
    cf3 = cfunction(F3,Void,(Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Int))

    F(M) = ccall((:_Z9c_updateUPvS_S_S_S_S_lPFvS_S_S_S_S_S_lE,"/home/zac/.julia/v0.5/SparseGrids/deps/libsparse.so"),Void,(Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64},Int32,Ptr{Void}),vec(M.Fval),vec(M.Jval),vec(M.X),vec(M.XP),vec(M.SP),vec(M.ProbWeights),length(M),cf3)
end

function vecind!(x::Expr,N,npf)
    if x.head == :ref
        if length(x.args)==3
            if x.args[1]==:(M.X)
                x.args[1]=:(X)
                x.args[2]=x.args[2]+(pop!(x.args)-1)*N
            elseif x.args[1]==:(M.Fval)
                x.args[1]=:(Fval)
                x.args[2]=x.args[2]+(pop!(x.args)-1)*N
            elseif x.args[1]==:(M.XP)
                x.args[1]=:(XP)
                x.args[2]=x.args[2]+(pop!(x.args)-1)*N*npf
            elseif x.args[1]==:(M.SP)
                x.args[1]=:(SP)
                x.args[2]=x.args[2]+(pop!(x.args)-1)*N*npf
            elseif x.args[1]==:(M.ProbWeights)
                x.args[1]=:(ProbWeights)
                x.args[2]=x.args[2]+(pop!(x.args)-1)*N
            end
        end
    else
        for i = 1:length(x.args)
            if isa(x.args[i],Expr)
                x.args[i] = vecind!(x.args[i],N,npf)
            end
        end
    end
    x
end
