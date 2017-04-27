import SparseGrids.grow!

"""
    Grow state space grid from node M.G[i].
"""
function grow!(M::Model, id, bounds = ones(Int, size(M.SP, 2)) * 10)

    bounds[map(v -> isa(v, Stochastic), State(M.variables))] = M.G.L[map(v -> isa(v, Stochastic), State(M.variables))]

    G = deepcopy(M.G)
    grow!(G, id, bounds)
    X = M.G(M.X, values(G))
    G, S, SP, ProbWeights = initS(State(M.variables), typeof(M.G).parameters[2], G)

    M.G     = G
    M.SP    = SP
    M.ProbWeights = ProbWeights
    M.X     = X
    M.XP    = zeros(size(SP, 1), size(M.XP, 2))
    M.Fval  = zeros(length(G), size(M.Fval, 2))
    initD(M)
end
