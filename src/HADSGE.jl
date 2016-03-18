module HADSGE
using Calculus,Base.Threads
import Base:names,length,getindex,display,find,parse,spzeros,kron,*,promote_type,setindex!
import SparseGrids:NGrid,CC,ndgrid
include("utils.jl")
include("parse.jl")
include("calculus.jl")
include("parsemodel.jl")
include("model.jl")
include("initmodel.jl")
include("solve.jl")
include("aggregate.jl")
include("simulation.jl")
include("threads.jl")
# include("genF.jl")

export Model,solve,∫,updateA

end
