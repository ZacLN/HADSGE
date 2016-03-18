type Result
    lb::Int
    M::Model
    State::Array{Float64,2}
end

results = Result[]
for Lb ∈ 4:8
M=Model(:[ 1-R*β*Expect(λ[+1])/λ
           Uh-λ*η
],:[b       = (-2,10.,$Lb)
    η       = (1,0.9,0.1,1)
],:[b       = (-2,20.,b*0.95)
    h       = (0,1,0.7)
    c       = h*η+R*b[-1]-b
    λ       = c^-σc
    Uh      = ϕh*(1-h)^-σh
    B       = ∫(b,0.0)
    H       = ∫(h*η,0.3)
],:[β       = 0.98
    σc      = 2.5
    ϕh      = 2.0
    σh      = 2.0
    R       = 1.0166])
solve(M)
updateA(M)


@time HADSGE.updateT(M);HADSGE.updatetd(M);∫(M,:b)
@time HADSGE.updateT1(M);HADSGE.updatetd(M);∫(M,:b)
@time HADSGE.updateT2(M);HADSGE.updatetd(M);∫(M,:b)
@profile HADSGE.updateT2(M)


T,nh = 15000,100
state = zeros(nh)*[0 1]
M[:η](zid)
zid = ones(Int,nh)
State = zeros(T,2)
for t = 1:T
    M[:η](zid)
    state= [M(:b,state) M[:η].x[zid]]
    State[t,:]=mean(state,1)
    if mod(t,2000)==0
        state = repmat(state,4)
        zid=repmat(zid,4)
        println(mean(State[t-1500:t,1]),"  ",std(State[t-1500:t,1]))
    end
end
push!(results,Result(lb,M,State))
end
