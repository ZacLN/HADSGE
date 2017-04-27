import SparseGrids.cc_g
abstract Variable
abstract Idiosyncratic
abstract Collective
abstract State <: Variable

type Endogenous <: State
    name::Symbol
    x::Vector{Float64}
    lom::Variable
    bounds::Vector{Float64}
    Endogenous(v::Expr, v1::Expr) = new(v.args[1], cc_g(v.args[2].args[3] + 1) * (v.args[2].args[2] - v.args[2].args[1]) + v.args[2].args[1], parse(v1), v.args[2].args[1:2])
end

abstract Stochastic <: State
type AR{S} <: Stochastic
    name::Symbol
    μ::Float64
    ρ::Float64
    σ::Float64
    x::Vector{Float64}
    T::Array{Float64, 2}
    bounds::Vector{Float64}

    function AR(v::Expr)
        x, T = rouwenhorst(v.args[2].args[1:3]..., SparseGrids.M(1 + v.args[2].args[4]))
        new{v.head == :(:=) ? Collective : Idiosyncratic}(v.args[1], v.args[2].args[1], v.args[2].args[2], v.args[2].args[3], x, T, collect(extrema(x)))
    end
end

type Markov{S} <: Stochastic
    name::Symbol
    x::Vector{Float64}
    T::Array{Float64, 2}
    bounds::Vector{Float64}
    function Markov(v::Expr)
        new{v.head == :(:=) ? Collective : Idiosyncratic}(v.args[1], cc_g(v.args[2].args[3] + 1) * (v.args[2].args[1].args[2] - v.args[2].args[1].args[1]) + v.args[2].args[1].args[1], eval(current_module(), v.args[2].args[2]), eval(current_module(), v.args[2].args[1]))
    end
end

type Policy <: Variable
    name::Symbol
    init::Union{Float64, Expr}
    bounds::Vector{Float64}
    Policy(v::Expr) = new(v.args[1], v.args[2].args[3], (x -> x == :Inf ? 1e10 : x).(v.args[2].args[1:2]))
end

abstract Static <: Variable
type Exogenous <: Static
    name::Symbol
    init::Float64
    bounds::Vector{Float64}
    Exogenous(v::Expr) = new(v.args[1], v.args[2], [-Inf,Inf])
end

type Dependant <:Static
    name::Symbol
    e::Expr
    update::Function
    bounds::Vector{Float64}
end

type Aggregate <: Static
    name::Symbol
    target::Variable
    init::Float64
    bounds::Vector{Float64}
end

type Future <: Variable
    name::Symbol
    target::Variable
    Future(name::Symbol, target::Variable) = new(name, target)
end





(::Type{T}){T<:Variable}(V::Array{Variable}, rev::Bool = false) = rev ?  filter(x -> !isa(x, T), V) : filter(x -> isa(x, T), V)
names(V::Vector{Variable}) = map(x -> x.name, V)
function getindex(list::Vector{Variable}, v::Symbol)
    for i = 1:length(list)
        list[i].name == v && (return list[i])
    end
end
function find(list::Vector{Variable}, v::Symbol)
    for i = 1:length(list)
        list[i].name == v && (return i)
    end
end
display(v::Variable) = print((typeof(v) == Aggregate ? "∫" : ""), v.name, ((typeof(v) == Endogenous || (typeof(v) == Aggregate && typeof(v.target) == Endogenous)) ? "[-1]" : "[0]"))
display(v::Future) = print(v.name, "[1]")
display(V::Vector{Variable}) = [(length(T(V)) > 0 && (print(split(string(T), ".")[end][1:4], ":");[(print("  ");display(v)) for v  in T(V)]);println()) for T in [Endogenous, Stochastic, Policy, Exogenous, Aggregate, Dependant, Future]]


length{T<:State}(v::T) = length(v.x)

isag{T}(x::Markov{T}) = T == Collective
isag{T}(x::AR{T}) = T == Collective
isag(x) = false

function timeof(v::Variable)
    if isa(v, Endogenous)
        return -1
    elseif isa(v, Stochastic) || isa(v, Policy) ||  isa(v, Exogenous) || isa(v, Dependant)
        return 0
    elseif isa(v, Aggregate)
        return timeof(v.target)
    elseif isa(v, Future)
        return 1
    end
end

function parseex(ex::Expr)
    @assert isa(ex.args[1], Symbol)
    if isa(ex.args[2], Number)
        return Exogenous
    elseif ex.args[2].head == :tuple
        return Policy
    elseif isa(ex.args[2], Expr)
        if ex.args[2].head == :call && ex.args[2].args[1] == :∫
            return Aggregate
        else
            return Dependant
        end
    else
        error("can't parse $ex")
    end
end
parse(ex::Expr) = parseex(ex)(ex)


"""
    parsevars(foc,states,vars,params)
Parses model equations and variables.
"""
function parsevars(foc::Expr, states::Expr, vars::Expr, params::Expr)
    parameters   = Dict{Symbol, Float64}(zip([x.args[1] for x in params.args], [x.args[2] for x in params.args]))
    variables = Variable[]
    Dlist = Dict{Expr, Expr}()

    for v ∈ vars.args
        if isa(v.args[2], Expr) && v.args[2].args[1] != :∫ && v.args[2].head != :tuple
            p = addindex!(v.args[1]) => subs(addindex!(subs(v.args[2], parameters)), Dlist)
            push!(Dlist, p)
            push!(Dlist, tchange(p[1], 1) => tchange(p[2], 1))
            push!(variables, Dependant(p[1].args[1], p[2], x -> x, [-Inf,Inf]))
        end
    end

    for v ∈ states.args
        if length(v.args[2].args) == 3 && reduce(&, map(x -> isa(x, Number), v.args[2].args))
            for v1 in vars.args
                if v1.args[1] == v.args[1]
                    push!(variables, Endogenous(v, v1))
                    push!(variables, variables[end].lom)
                end
            end
        elseif length(v.args[2].args) == 3
            push!(variables, Markov{Idiosyncratic}(v))
        elseif length(v.args[2].args) == 4 && reduce(&, map(x -> isa(x, Number), v.args[2].args))
            push!(variables, AR{Idiosyncratic}(v))
        end
    end

    for v in vars.args
        if !in(v.args[1], names(variables)) && !in(parseex(v), [Dependant,Aggregate])
            push!(variables, parseex(v)(v))
        elseif parseex(v) == Aggregate
            if isa(v.args[2].args[2], Symbol)
                if in(v.args[2].args[2], names(State(variables)))
                    push!(variables, Aggregate(v.args[1], State(variables)[v.args[2].args[2]], v.args[2].args[3], [-Inf,Inf]))
                    push!(variables, Aggregate(v.args[1], State(variables, true)[v.args[2].args[2]], v.args[2].args[3], [-Inf,Inf]))
                elseif in(v.args[2].args[2], names(Dependant(variables)))
                    push!(variables, Aggregate(v.args[1], Dependant(variables)[v.args[2].args[2]], v.args[2].args[3], [-Inf,Inf]))
                elseif in(v.args[2].args[2], names(Policy(variables)))
                    push!(variables, Aggregate(v.args[1], Policy(variables)[v.args[2].args[2]], v.args[2].args[3], [-Inf,Inf]))
                else
                    error("Integral target $(v.args[2].args[2]) for variable $(v.args[1]) not found.")
                end
            elseif isa(v.args[2].args[2], Expr)
                p = gensym(v.args[1]) => subs(addindex!(subs(v.args[2].args[2], parameters)), Dlist)
                push!(variables, Dependant(p[1], p[2], x -> x, [-Inf,Inf]))
                push!(variables, Aggregate(v.args[1], Dependant(variables)[end], v.args[2].args[3], [-Inf,Inf]))
            end
        end
    end
    variables = vcat(Endogenous(variables), Stochastic(variables), Policy(variables), Exogenous(variables), Aggregate(variables), Dependant(variables))

    f2, j2 = parsefoc(foc, variables, Dlist, parameters)

    return variables, parameters, f2, j2
end



function hardloc(variables::Vector{Variable})
    ns = length(State(variables))
    hloc = Dict{Expr, Expr}()

    for v in variables
        if isa(v, Endogenous)
            push!(hloc, Expr(:ref, v.name, -1) => :(M.X[i, $(findfirst(variables, v))]))
        elseif isa(v, Stochastic)
            push!(hloc, Expr(:ref, v.name, 0) => :(M.X[i, $(findfirst(variables, v))]))
            push!(hloc, Expr(:ref, v.name, 1) => :(M.SP[i + (j - 1) * length(M), $(findfirst(State(variables), v))]))
        elseif isa(v, Policy)
            push!(hloc, Expr(:ref, v.name, 0) => :(M.X[i, $(ns + findfirst(State(variables, true), v))]))
        elseif isa(v, Static)
            push!(hloc, Expr(:ref, v.name, timeof(v)) => :(M.X[i, $(ns + findfirst(State(variables, true), v))]))
        elseif isa(v, Future)
            push!(hloc, Expr(:ref, v.name, 1) => :(M.XP[i + (j - 1) * length(M), $(findfirst(Future(variables), v))]))
        end
    end
    return hloc
end

function parsefoc(foc::Expr, variables::Vector{Variable}, Dlist::Dict, parameters::Dict)
    f = subs(addindex!(subs(foc, parameters)), Dlist)
    @assert f.head == :vcat || f.head == :vect
    list = :([])
    for i = 1:length(f.args)
        f.args[i], list = getexpectation(f.args[i], list, length(list.args) + 1)
    end
    f2 = deepcopy(f)
    subs!(f2, Dict(zip([Expr(:ref, :Expect, i) for i = 1:length(list.args)], list.args)))
    filter!((k, v) -> k.args[2] < 1, Dlist)
    j2  = jacobian(f2, [Expr(:ref, v, 0) for v in names(Policy(variables))])

    for v in sort(setdiff(Symbol[x.args[1] for x in filter(x -> x.args[2] == 1, getv(f2))], names(Stochastic(variables))))
        push!(variables, Future(v, State(variables, true)[v]))
    end

    hloc = hardloc(variables)
    nP  = reduce(*, map(length, Stochastic(variables)))
    F = simplifyindices!(addpweights!(subs(f2, hloc), nP))
    J = simplifyindices!(vec(addpweights!(subs(j2, hloc), nP)))
    return F, J
end

function getexpectation(x, list, ieq)
    if isa(x, Expr)
        if x.head == :call && x.args[1] == :Expect
            push!(list.args, :ProbWeights * x.args[2])
            x.head = :ref
            x.args[2] = ieq
        else
            for i = 1:length(x.args)
                x.args[i], list = getexpectation(x.args[i], list, ieq)
            end
        end
    end
    return x, list
end
