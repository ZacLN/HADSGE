import SparseGrids:nXtoU,cc_dM,cc_M,cc_dg,cc_bf_q,cc_bf_l
import HADSGE:State,Stochastic
function updateTsg(M)
    D         = length(M.G.L)
    mL        = maximum(M.G.L)
    J         = ones(Int, mL + 1, D)
    B         = zeros(mL + 1, D)
    T = spzeros(length(M), length(M))
    nstoc = size(M.ProbWeights, 2)
    SP = clamp!(reshape(nXtoU(M.SP, M.G.bounds), length(M), size(M.ProbWeights, 2), size(M.SP, 2)), 0, 1)
    isstoc = [isa(x, Stochastic) for x ∈ State(M.variables)]
    for i ∈ 1:length(M)
        X = SP[i,:, :]
        for k ∈ 1:nstoc
            x = X[k,:]
            for d = (1:D)[isstoc]
                for l = 1:mL + 1
                    j     = clamp(round(Int, x[d] * (cc_dM(l)) + 1 / 2), 1, cc_dM(l))
                    B[l,d]     = cc_bf_l(x[d], cc_dg(l, j), Int16(cc_M(15)))
                    J[l,d]  = j
                end
            end
            mLstoc = Int[findfirst(B[:, isstoc][:, ii]) for ii = 1:length(sum(isstoc))]

            for d = (1:D)[!isstoc]myd
                for l = 1:mL + 1
                    j     = clamp(round(Int, x[d] * (cc_dM(l)) + 1 / 2), 1, cc_dM(l))
                    # B[l,d]     = cc_bf_q(x[d],cc_dg(l,j),Int16(cc_M(l)))
                    B[l,d]     = cc_bf_l(x[d], cc_dg(l, j), Int16(cc_M(mL + D - sum(mLstoc))))
                    J[l,d]  = j
                end
            end

            res = []
            for ii = 1:size(M.G.coverings, 1)
                bv, jv = [B[M.G.coverings[ii,d],d] for d = 1:D], [J[M.G.coverings[ii,d],d] for d = 1:D]
                b  = B[M.G.coverings[ii,D],D] * B[M.G.coverings[ii,1],1]
                id1 = J[M.G.coverings[ii,D],D] - 1
                for d = D - 1:-1:2
                    b *= B[M.G.coverings[ii,d],d]
                    id1 = id1 * M.G.coverings_dM[ii,d] + (J[M.G.coverings[ii,d],d] - 1)
                end
                id1 = (J[M.G.coverings[ii,1],1] - 1) + M.G.coverings_dM[ii,1] * id1 + 1 + M.G.coveringsloc[1][ii] - 1

                M.G.coverings[ii,isstoc] == mLstoc && b > 0 && push!(res, (id1, b))
            end
            for r in res
                 T[r[1],i] += r[2] * M.ProbWeights[i,k]
            end
        end
    end
    T
end
