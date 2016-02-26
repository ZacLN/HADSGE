ARSim(x::Vector{Float64},e::AR,s=randn(size(x))) = clamp((1-e.ρ)*e.μ + e.ρ*x+s*e.σ,e.x[1],e.x[end])
function MarkovSim!(ID::Vector{Int},s::Stochastic,r::Vector{Float64}=rand(length(ID)))
    csT=  cumsum(s.T,2)
    for i = 1:length(ID)
        j = 1
        while csT[ID[i],j]<r[i]
            j+=1
        end
        ID[i] = j
    end
    return ID
end


(s::AR)(ID::Vector{Int},r::Vector{Float64}=rand(length(ID))) = MarkovSim!(ID,s,r)
(s::AR)(ID::Int,r::Vector{Float64}=rand(length(ID))) = MarkovSim!([ID],s,r)[1]
(s::Markov)(ID::Vector{Int},r::Vector{Float64}=rand(length(ID))) = MarkovSim!(ID,s,r)
