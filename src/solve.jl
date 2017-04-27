function updateXP(M::Model)
    fvar = Future(M.variables)
    A = hcat([M[v.target.name, 0] for v in fvar]...)
    M.G(A, M.SP, M.XP)
    for i in 1:length(fvar)
        lb = fvar[i].target.bounds[1]::Float64
        ub = fvar[i].target.bounds[2]::Float64
        @inbounds @simd for j = 1:size(M.XP, 1)
            if M.XP[j,i] < lb
                M.XP[j,i] = lb
            elseif M.XP[j,i] > ub
                M.XP[j,i] = ub
            end
        end
    end
end

function forceboundsU(M::Model)
    ns = length(State(M.variables))
    for v ∈ Policy(M.variables)
        lb = v.bounds[1]::Float64
        ub = v.bounds[2]::Float64
        vid = findnext(M.variables, v, ns + 1)
        @inbounds @simd for i = 1:length(M)
            if isnan(M.X[i,vid])
                error("NaN in $(v.name) at $i")
            elseif M.X[i,vid] < lb
                M.X[i,vid] = lb
            elseif M.X[i,vid] > ub
                M.X[i,vid] = ub
            end
        end
    end
end


function solve(M::Model, n::Int = 1000, ϕ::Float64 = 0.8;disp::Int = 10000, crit = 1e-6)
    MaxError = maximum(abs.(M.Fval), 1)
    SError = sum(abs.(M.Fval) .> crit * 5, 1) / length(M)
    updateallD(M)
    SPids  = [(find(State(M.variables), s.name), find(State(M.variables, true), s.name))::Tuple{Int, Int} for s in Endogenous(M.variables)]
    ns = length(M.G.L)
    nG = length(M)

    ne = length(MaxError)
    Fval = zeros(ne)
    Jval = zeros(ne, ne)

    for iter = 1:n
        (iter <= 5 || mod(iter, 20) == 0 || iter < 5) && (M.F(M);MaxError = maximum(abs.(M.Fval), 1);SError = sum(abs.(M.Fval) .> crit * 5, 1) / length(M))
        for s in SPids
            for j = 1:size(M.ProbWeights, 2)
                for i = 1:nG
                    @inbounds M.SP[i + (j - 1) * nG, s[1]] = M.X[i, ns + s[2]]
                end
            end
        end
        updateXP(M)
        # M.Fall(M,ϕ)
        @threadsfixed [Fval, Jval] for i = 1:nG
            M.Finner(M, i, nG, Fval, Jval, ϕ)
        end
        forceboundsU(M)
        mod(iter, disp) == 0 && println(round(log10.(MaxError), 2), "  ", round(SError, 2))
        iter > min(n, 5) && (maximum(MaxError) < crit || all(SError .< 0.025)) && (print(iter);break)
    end
    # displaysolve(M,"done",n,MaxError,SError,time()-tstart)
    updateallD(M)
end


function updateallD(M::Model)
    for v in Dependant(M.variables)
        v.update(M)
    end
end
