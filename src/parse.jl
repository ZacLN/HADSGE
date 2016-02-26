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

function subs!(x::Expr,list::Dict)
    for i = 1:length(x.args)
        if in(x.args[i],keys(list))
            x.args[i] = list[x.args[i]]
        elseif isa(x.args[i],Expr)
            subs!(x.args[i],list)
        end
    end
end

function subs(x1::Expr,list::Dict)
    x = deepcopy(x1)
    subs!(x,list)
    return x
end

function subs!(x::Expr,s::Pair)
    for i = 1:length(x.args)
        if x.args[i]==s.first
            x.args[i] = s.second
        elseif isa(x.args[i],Expr)
            subs!(x.args[i],s)
        end
    end
end

function subs(x1::Expr,s::Pair)
    x = deepcopy(x1)
    subs!(x,s)
    return x
end


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

function tchange!(x::Expr,t::Int,ignore=operators)
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

function addpweights!(x,nP)
  if typeof(x) == Expr
    if (x.head == :call) && (x.args[1] == :*) && (x.args[2]==:ProbWeights)
    #   return Expr(:call,:+,[subs(x.args[3],:j=>j)*:(M.ProbWeights[i,$j]) for j = 1:nP]...)
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



import Base:+,-,*,/,vec

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

function +{T<:Union{Expr,Symbol,Number},S<:Union{Expr,Symbol,Number}}(x1::T,x2::S)
    :(+($x1,$x2))
end

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
