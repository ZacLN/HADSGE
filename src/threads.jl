module pHADSGE
using Base.Threads
import HADSGE


function Model(foc,states,vars,params)
    variables,parameters,F,J=HADSGE.parsevars(foc,states,vars,params)
    G,S,SP,ProbWeights=HADSGE.initS(HADSGE.State(variables))
    X   = [S zeros(length(G),length(HADSGE.Policy(variables))) zeros(length(G),length(HADSGE.Static(variables)))]
    XP  = zeros(length(G)*size(ProbWeights,2),length(HADSGE.Future(variables)))

    tG = HADSGE.ndgrid(Vector{Float64}[s.x for s in HADSGE.State(variables)]...)
    distribution = HADSGE.Distribution(tG,zeros(size(tG[1]))+1/length(tG[1]),spzeros(length(tG[1]),length(tG[1])))

    HADSGE.subs!(F,:(length(M))=>length(G))
    F = :($(gensym(:F))(M::Model) = @fastmath $(HADSGE.buildfunc(F,:(M.Fval))))
    HADSGE.subs!(J,:(length(M))=>length(G))
    Jp = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(HADSGE.buildJ2(J)))
    J = :( $(gensym(:J))(M::Model,i::Int) = @fastmath $(HADSGE.buildJ(J)))

    M = HADSGE.Model(G,SP,ProbWeights,X,XP,eval(F),eval(J),zeros(length(G),length(HADSGE.Policy(variables))),zeros(length(HADSGE.Policy(variables)),length(HADSGE.Policy(variables))),variables,distribution,(F,J,eval(Jp)))

    HADSGE.initX(M)
    HADSGE.initD(M)
    return M
end

function updateU(M::HADSGE.Model,ϕ=0.9)
    M.F(M)
    n = size(M.Fval,2)
    ns = length(M.G.L)
    @threads for i =1:length(M)
        Jval=M.temp[3](M,i)
        x = Jval\M.Fval[i,:]
        @fastmath @simd for j = 1:n
            @inbounds M.X[i,ns+j] -= ϕ*x[j]
        end
    end
end

function solve(M::HADSGE.Model,n::Int=1000,ϕ::Float64=0.8;disp::Int=10000,crit=1e-6)
    MaxError = maximum(abs(M.Fval),1)
    SError = sum(abs(M.Fval).>crit*5,1)/length(M)
    HADSGE.updateallD(M)

    for iter = 1:n
        (mod(iter,20)==0 || iter < 5) && (MaxError = maximum(abs(M.Fval),1);SError = sum(abs(M.Fval).>crit*5,1)/length(M))
        HADSGE.updateSP(M)
        HADSGE.updateXP(M)
        updateU(M,ϕ)
        HADSGE.forceboundsU(M)
        mod(iter,disp)==0 && println(round(log10(MaxError),2),"  ",round(SError,2))
        iter>min(n,5)  && (maximum(MaxError)<crit || all(SError.<0.025)) && (print(iter);break)
    end
    HADSGE.updateallD(M)
end

end  # module
