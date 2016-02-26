import Calculus:differentiate,jacobian,simplify,SymbolParameter,isminus


function simplify(::SymbolParameter{:-}, args)
	if length(args) == 2 && args[2]==0
		return args[1]
	end
	# Remove any 0's in a subtraction
    args = map(simplify, filter(x -> x != 0, args))
    if length(args) == 0
        return 0
    # Special Case: simplify(:(x - x)) == 0
    elseif length(args) == 2 && args[1] == args[2]
        return 0
    # Special Case: simplify(:(x - (-y))) == x + y
    elseif length(args) == 2 && isminus(args[2])
        return Expr(:call, :+, args[1], args[2].args[2])
    else
	        return Expr(:call, :-, args...)
    end
end






function differentiate(ex::Expr,wrt::Expr)
	if ex.head==:vect || ex.head ==:vcat
		return differentiate(SymbolParameter(:vect), ex.args[1:end], wrt)
	elseif ex.head == :ref
        return ex==wrt ? 1 : 0
    elseif ex.head != :call
		error("Unrecognized expression $ex")
	end
    return simplify(differentiate(SymbolParameter(ex.args[1]), ex.args[2:end], wrt))
end

function differentiate(ex::Expr, targets::Vector{Expr})
	n = length(targets)
	exprs = Array(Any, n)
	for i in 1:n
		exprs[i] = differentiate(ex, targets[i])
	end
	return exprs
end

differentiate(ex::Union{Number,Symbol},wrt::Expr) = 0

function differentiate(::SymbolParameter{:vect}, args, wrt)
	for i = 1:length(args)
		args[i] = differentiate(args[i],wrt)
	end
	return Expr(:vect,args...)
end

function jacobian(ex::Expr, targets::Vector)
    @assert ex.head==:vect || ex.head==:vcat
    exprs = Expr(:vcat)
    for i = 1:length(ex.args)
        push!(exprs.args,Expr(:row))
        for j = 1:length(targets)
            push!(exprs.args[i].args,simplify(differentiate(ex.args[i],targets[j])))
        end
    end
    return exprs
end

differentiate(::SymbolParameter{:max}, args, wrt) = Expr(:if,:($(args[1])>$(args[2])),differentiate(args[1],wrt),differentiate(args[2],wrt))
differentiate(::SymbolParameter{:min}, args, wrt) = Expr(:if,:($(args[1])<$(args[2])),differentiate(args[1],wrt),differentiate(args[2],wrt))



function simplify(ex::Expr)
	if in(ex.head,[:vcat,:row])
		for i = 1:length(ex.args)
			ex.args[i]=simplify(ex.args[i])
		end
		return ex
	end

	if ex.head != :call
        return ex
    end

    if all(Calculus.isnumber, ex.args[2:end]) && length(ex.args) > 1
        return eval(current_module(), ex)
    end
    new_ex = simplify(SymbolParameter(ex.args[1]), ex.args[2:end])
    while new_ex != ex
        new_ex, ex = simplify(new_ex), new_ex
    end
    return new_ex
end
