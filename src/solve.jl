function updateSP(M::Model)
    for s in Endogenous(M.variables)
        sid = find(State(M.variables),s.name)
        nG     = length(M.G)
        if isa(s.lom,Policy)
            uid = find(Policy(M.variables),s.name)
            for j = 1:size(M.ProbWeights,2)
                for i = 1:nG
                    @inbounds M.SP[i+(j-1)*nG,sid] = M.U[i,uid]
                end
            end
        elseif isa(s.lom,Static)
            xid = find(Static(M.variables),s.name)
            isa(s.lom,Dependant) && s.lom.update(M)
            for j = 1:size(M.ProbWeights,2)
                for i = 1:nG
                    @inbounds M.SP[i+(j-1)*nG,sid] = M.X[i,xid]
                end
            end
        end
    end
end

function updateXP(M::Model)
    fvar = HADSGE.Future(M.variables)
    A = hcat([M[v.target.name,0] for v in fvar]...)
    M.XP[:,:] = M.G(A,M.SP)
    for i in 1:length(fvar)
        clamp!(M.XP[:,i],fvar[i].target.bounds[1],fvar[i].target.bounds[2])
    end
end

function updateU(M::Model,ϕ=0.9)
    M.F(M)
    n = size(M.U,2)
    for i =1:length(M.G)
        M.J(M,i)
        x = M.Jval\M.Fval[i,:]
        @fastmath @simd for j = 1:n
            @inbounds M.U[i,j] -= ϕ*x[j]
        end
    end
end

function forceboundsU(M::Model)
    for j = 1:size(M.U,2)
        lb,ub = Policy(M.variables)[j].bounds[1],Policy(M.variables)[j].bounds[2]
        @inbounds @simd for i = 1:length(M)
            if M.U[i,j] < lb
                M.U[i,j] = lb
            elseif M.U[i,j] > ub
                M.U[i,j] = ub
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
        iter>min(n,5)  && (maximum(MaxError)<crit || all(SError.<0.05)) && (print(iter);break)
    end
    updateallD(M)
end


function updateallD(M::Model)
    for v in Dependant(M.variables)
        v.update(M)
    end
end
