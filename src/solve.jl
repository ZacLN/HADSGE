function updateSP(M::Model)
    ns = length(M.G.L)
    for s in Endogenous(M.variables)
        sid = find(State(M.variables),s.name)
        nG     = length(M.G)
        uid = find(State(M.variables,true),s.name)
        for j = 1:size(M.ProbWeights,2)
            for i = 1:nG
                @inbounds M.SP[i+(j-1)*nG,sid] = M.X[i,ns+uid]
            end
        end
    end
end

function updateXP(M::Model)
    fvar = Future(M.variables)
    A = hcat([M[v.target.name,0] for v in fvar]...)
    M.XP[:,:] = M.G(A,M.SP)
    for i in 1:length(fvar)
        clamp!(M.XP[:,i],fvar[i].target.bounds[1],fvar[i].target.bounds[2])
    end
end

function updateU(M::Model,ϕ=0.9)
    M.F(M)
    n = size(M.Fval,2)
    ns = length(M.G.L)
    for i =1:length(M.G)
        M.J(M,i)
        x = M.Jval\M.Fval[i,:]
        @fastmath @simd for j = 1:n
            @inbounds M.X[i,ns+j] -= ϕ*x[j]
        end
    end
end

function forceboundsU(M::Model)
    ns = length(State(M.variables))
    for v ∈ Policy(M.variables)
        lb,ub = v.bounds[1],v.bounds[2]
        vid=findnext(M.variables,v,ns)
        @inbounds @simd for i = 1:length(M)
            if M.X[i,vid] < lb
                M.X[i,vid] = lb
            elseif M.X[i,vid] > ub
                M.X[i,vid] = ub
            end
        end
    end
end


function solve(M::Model,n::Int=1000,ϕ::Float64=0.8;disp::Int=10000,crit=1e-6)
    MaxError = maximum(abs(M.Fval),1)
    SError = sum(abs(M.Fval).>crit*5,1)/length(M)
    updateallD(M)

    for iter = 1:n
        (mod(iter,20)==0 || iter < 5) && (MaxError = maximum(abs(M.Fval),1);SError = sum(abs(M.Fval).>crit*5,1)/length(M))
        updateSP(M)
        updateXP(M)
        updateU(M,ϕ)
        forceboundsU(M)
        mod(iter,disp)==0 && println(round(log10(MaxError),2),"  ",round(SError,2))
        iter>min(n,5)  && (maximum(MaxError)<crit || all(SError.<0.025)) && (print(iter);break)
    end
    updateallD(M)
end


function updateallD(M::Model)
    for v in Dependant(M.variables)
        v.update(M)
    end
end
