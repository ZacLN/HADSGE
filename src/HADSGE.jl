module HADSGE
using Calculus,Base.Threads
import Base:names,length,getindex,display,find,parse,spzeros,kron,*,promote_type,setindex!,Future
import SparseGrids:NGrid,ndgrid,@threadsfixed
using SparseGrids
include("utils.jl")
include("parse.jl")
# include("calculus.jl")
include("parsemodel.jl")
include("model.jl")
include("initmodel.jl")
include("solve.jl")
include("aggregate.jl")
include("simulation.jl")
include("adapt.jl")
# include("display.jl")


export Model,solve,∫,updateA

end
