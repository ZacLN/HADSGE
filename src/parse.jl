# Provides functions to parse model equations. Array access syntax is used to indicate the time period of functions.


# List of mathematical operators supported in model equations
operators = [:-
            :+
            :*
            :/
            :^
            :exp
            :log
            :Expect
            :max
            :min
            :∫
            :≤]


"""
    subs!(x,d)

Replaces instances of the keys of dictionary d with the values in expression x.
"""
function subs!(x::Expr,d::Dict)
    for i = 1:length(x.args)
        if in(x.args[i],keys(d))
            x.args[i] = d[x.args[i]]
        elseif isa(x.args[i],Expr)
            subs!(x.args[i],d)
        end
    end
end

function subs(x1::Expr,d::Dict)
    x=deepcopy(x1)
    subs!(x,d)
    return x
end

subs!(x::Expr,p::Pair) = subs!(x,Dict(p))
subs(x::Expr,p::Pair) = subs(x,Dict(p))

# function subs!(x::Expr,s::Pair)
#     for i = 1:length(x.args)
#         if x.args[i]==s.first
#             x.args[i] = s.second
#         elseif isa(x.args[i],Expr)
#             subs!(x.args[i],s)
#         end
#     end
# end
#
# function subs(x1::Expr,s::Pair)
#     x = deepcopy(x1)
#     subs!(x,s)
#     return x
# end

"""
    addindex!(x)

Adds a time index [0] to unindexed variables in x.
e.g. c^-2.5 -> c[0]^-2.5
"""
function addindex!(x,ignore=operators)
    if typeof(x) == Expr
        if x.head == :ref
            return x
        else
            for i = 1:length(x.args)
                x.args[i]=addindex!(x.args[i],ignore)
            end
        end
    elseif typeof(x) == Symbol
        if !in(x,ignore)
            x= :($x[0])
            return x
        end
    end
    return x
end

addindex(x,ignore=operators) = addindex!(deepcopy(x),ignore)

"""
    removeindex!(x)
Removes all time indices from model expression x.
"""
function removeindex!(x)
    if isa(x,Expr)
        if x.head == :ref
            return x.args[1]
        else
            for i = 1:length(x.args)
                x.args[i]=removeindex!(x.args[i])
            end
        end
    end
    return x
end

removeindex(x,ignore=operators) = removeindex!(deepcopy(x),ignore)

"""
    tchange!(x,t)
Either shifts the time indices in model expression x forward or backwards depending on
whether t = 1 or t = -1 respectively.
"""
function tchange!(x::Expr,t::Int,ignore=operators)
    t != -1 && t!=1 && return
    if x.head == :ref
        if (x.args[2] == -1 && t==1) || (x.args[2] == 1 && t==-1)
            x = :($(x.args[1])[0])
        elseif x.args[2]==0
            x.args[2] = t
        end
    else
        for i = 1:length(x.args)
            if isa(x.args[i],Expr)
                x.args[i] = tchange!(x.args[i],t,ignore)
            end
        end
    end
    x
end

tchange(x::Expr,t::Int,ignore=operators) = tchange!(deepcopy(x),t,ignore)


"""
    getv(x)
Returns all symbols in the expression x, excluding mathematical operators.
"""
function getv(x,list=Expr[],ignore=operators)
    if isa(x,Expr)
        if x.head==:ref
            push!(list,x)
        else
            for i = 1:length(x.args)
                list = getv(x.args[i],list,ignore)
            end
        end
    end
    return unique(list)
end


# function vecind!(x::Expr,n::Int)
#     if x.head == :ref
#         if length(x.args)==3
#             x.args[2]=x.args[2]+(pop!(x.args)-1)*n
#         end
#     else
#         for i = 1:length(x.args)
#             if isa(x.args[i],Expr)
#                 x.args[i] = vecind!(x.args[i],n)
#             end
#         end
#     end
#     x
# end


"""
    simplifyindices!(x)
Simplifies expressions in array access portions of an arbitrary expression x.
"""
function simplifyindices!(x)
    if isa(x,Expr)
        if x.head==:ref
            for i = 2:length(x.args)
                x.args[i] = simplify(x.args[i])
            end
        else
            for a in x.args
                a = simplifyindices!(a)
            end
        end
    end
    return x
end


"""
    addpweights!(x,nP)
Replaces instances of ProbWeights*(...) in model equations x with a weighted sum over the probabilities of future realisations of the models stochastic variables. nP defines the number of possible future realisations.
"""
function addpweights!(x,nP)
  if typeof(x) == Expr
    if (x.head == :call) && (x.args[1] == :*) && (x.args[2]==:ProbWeights)
      return sum([subs(x.args[3],:j=>j)*:(M.ProbWeights[i,$j]) for j = 1:nP])
    else
      for i = 1:length(x.args)
        x.args[i]=addpweights!(x.args[i],nP)
      end
    end
  end
  return x
end

addpweights(x,nP) = addpweights!(deepcopy(x),nP)

"""
    findrepeated(x)
Returns a list of repeated sub-expressions in the model equations x.
"""
function findrepeated(x,list=Dict())
    if isa(x,Expr) && x.head!=:ref
        if in(x,keys(list))
            list[x]=list[x]+1
        else
            list[x] = 1
        end

        for i = 1:length(x.args)
            list = findrepeated(x.args[i],list)
        end
    end
    return (list)
end

function Base.length(ex,n=0)
    if isa(ex,Expr)
        n+=sum(map(length,ex.args))
    else
        n+=1
    end
    return n
end




# The following section allows arithmetic between expressions, symbols and numbers simplifying julia function generation from the parsed model equations.

import Base:+,-,*,/,vec

# When applied to an expression that defines a matrix will return the elements as an expression defining the column vector of all elements.
function vec(ex::Expr)
	nex = :([])
	nc = length(ex.args)
	nr = length(ex.args[1].args)
	@assert ex.head==:vcat

	for r in ex.args
		@assert r.head == :row && length(r.args) == nr
	end

    for ri = 1:nr
    	for ci = 1:nc
			push!(nex.args,ex.args[ci].args[ri])
		end
	end
	return nex
end



+{T<:Union{Expr,Symbol,Number}}(x::T) = x
+{T<:Union{Expr,Symbol,Number},S<:Union{Expr,Symbol,Number}}(x1::T,x2::S) =    :(+($x1,$x2))

function +{T<:Union{Expr,Symbol},S<:Float64}(x1::T,x2::Vector{S})
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]+x1))
    end
    return out
end

function +{T<:Union{Expr,Symbol},S<:Float64}(x2::Vector{S},x1::T)
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]+x1))
    end
    return out
end

-{T<:Union{Expr,Symbol,Number}}(x::T) = -1*x

function -{T<:Union{Expr,Symbol,Number},S<:Union{Expr,Symbol,Number}}(x1::T,x2::S)
    :(-($x1,$x2))
end

function -{T<:Union{Expr,Symbol},S<:Float64}(x1::T,x2::Vector{S})
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]-x1))
    end
    return out
end

function -{T<:Union{Expr,Symbol},S<:Float64}(x2::Vector{S},x1::T)
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]-x1))
    end
    return out
end



function *{T<:Union{Expr,Symbol,Number},S<:Union{Expr,Symbol,Number}}(x1::T,x2::S)
    :(*($x1,$x2))
end

function *{T<:Union{Expr,Symbol},S<:Float64}(x1::T,x2::Vector{S})
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]*x1))
    end
    return out
end

function *{T<:Union{Expr,Symbol},S<:Float64}(x2::Vector{S},x1::T)
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]*x1))
    end
    return out
end


function /{T<:Union{Expr,Symbol,Number},S<:Union{Expr,Symbol,Number}}(x1::T,x2::S)
    :(/($x1,$x2))
end

function /{T<:Union{Expr,Symbol},S<:Float64}(x1::T,x2::Vector{S})
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x1/x2[i]))
    end
    return out
end

function /{T<:Union{Expr,Symbol},S<:Float64}(x2::Vector{S},x1::T)
    out = :([])
    for i = 1:length(x2)
        push!(out.args,simplify(x2[i]/x1))
    end
    return out
end
