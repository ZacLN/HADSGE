if false
M=Model(:[1-R*β*Expect(λ[+1])/λ
           Uh/λ-η
],:[b       = (-2,10.,8)
    η       = (1,0.9,0.1,1)
],:[b       = (-2,10.,b*0.95)
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

HADSGE.initdisplay(M)
solve(M,disp=50,crit=1e-12)

end



function initdisplay(M::Model)
    W = Window()
    Blink.addstyle!(W,"""
    body {background-color:#2d2d2d;font-size: 10;color: white}
    .tooltip {
        position: relative;display: inline-block;border-bottom: 1px dotted black;}
    .tooltip .tooltiptext {visibility: hidden;width: 120px;background-color: black;color: #fff;text-align:center;padding: 5px 0;border-radius:6px;position:absolute;z-index:1;}
    .tooltip:hover .tooltiptext {visibility: visible;}
    .field-tip {position:relative;cursor:default}
    .field-tip .tip-content {position:absolute;top:-10px;right:9999px;width:200px;margin-right:-220px;padding:10px;color:#fff;background:#333;opacity:0;
        -webkit-transition:opacity 250ms ease-out;
           -moz-transition:opacity 250ms ease-out;
            -ms-transition:opacity 250ms ease-out;
             -o-transition:opacity 250ms ease-out;
                transition:opacity 250ms ease-out;}
    .field-tip .tip-content:before {
        content:' ';position:absolute;top:50%;
        left:-16px;width:0;height:0;margin-top:-8px;border:8px solid transparent;border-right-color:#333;}
    .field-tip:hover .tip-content {right:-20px;opacity:1;}
    """)
    @js W resizeTo(495,400)
    @js W moveTo(866,5)




    # Summary
    io=IOBuffer()
    for v in State(M.variables)
        if isa(v,Stochastic)
            if isa(v,AR)
                desc="[$(v.name) = $(round(1-v.ρ*v.μ,4)) + $(v.ρ)*$(v.name)[-1] + $(v.σ)*ϵ]  $(v.bounds)"
            elseif isa(v,Markov)
                desc="[$(v.name) = T($(v.name)[-1])] $(v.bounds)"
            end
        else
            desc="Endogenous: $(v.bounds)"
        end
        desc*= " $(length(v.x))"

        print(io,"""<span class="field-tip">
            $(v.name)
            <span class="tip-content">$desc</span>
        </span>  """)
    end
    states = takebuf_string(io)
    close(io)

    io=IOBuffer()
    for v in Policy(M.variables)
        desc= " $(map(x->round(x,3),extrema(M[v.name,0])))"
        print(io,"""<span class="field-tip">
            $(v.name)
            <span class="tip-content">$desc</span>
        </span>  """)
    end
    policies = takebuf_string(io)
    close(io)


    summary = """
    <table style="width:100%">
      <tr>
        <td style="width:100%;heigh:100%;background-color:#3d3d3d;width:60%">
        <table style="width:60%;background-color:#3d3d3d">
          <tr>
            <td style="font-size:12;width:20%">States</td>
            <td style="font-size:10">$states</td>
          </tr>
          <tr>
            <td style="font-size:12">Policy</td>
            <td style="font-size:10">$policies</td>
          </tr>
        </table>
        </td><td style="width:40%"></td>
      </tr>
    </table>
    """

    #Equations

    io = IOBuffer()
    print(io,"""<table style="width:100%">""")
    for i in 1:length(M.summary.equations.args)
        print(io,"""
        <tr>
        <td style="font-size:11;background-color:#3d3d3d">
        $i
        <td>
        <td style="font-size:11;background-color:#3d3d3d">
        $(M.summary.equations.args[i])
        </td>
        <td id="fval$i" style="font-size:11;background-color:#3d3d3d">        $(signif(median(M.Fval[:,i]),2))  ($(signif(maximum(M.Fval[:,i]),2)))</td></tr>""")
    end
    print(io,"</table>")
    eqs = takebuf_string(io)
    close(io)




    # Progress meter

    progress= """
    <h3> Solve Step: <text id="solvestep"></text></h3>
    <div>Time: <text id="Telapsed"></text></div>
    <div>Iterations: <text id="iiter">0</text>/<text id="Niter">0</text></div>
    """

    # progress ="""
    # <h3>Progress (<text id="solvestep"></text>)</h3>
    # <table style="width:100%;background-color:#3d3d3d">
    # <tr><td id="time"></td><td id="iter"></td></tr>
    # </table>
    # """
    # <table style="width:100%;background-color:#3d3d3d">
    #   <tr style="font-size:12">
    #     <td style="width:8%">tid</td>
    #     <td style="width:8%">1</td>
    #     <td style="width:8%">2</td>
    #     <td style="width:8%">3</td>
    #     <td style="width:8%">4</td>
    #     <td style="width:8%">5</td>
    #     <td style="width:8%">6</td>
    #     <td style="width:8%">7</td>
    #     <td style="width:8%">8</td>
    #   </tr>
    #   <tr style="font-size:9">
    #     <td style="width:8%"></td>
    #     <td style="width:8%" id="progress1">0</td>
    #     <td style="width:8%" id="progress2">0</td>
    #     <td style="width:8%" id="progress3">0</td>
    #     <td style="width:8%" id="progress4">0</td>
    #     <td style="width:8%" id="progress5">0</td>
    #     <td style="width:8%" id="progress6">0</td>
    #     <td style="width:8%" id="progress7">0</td>
    #     <td style="width:8%" id="progress8">0</td>
    #   </tr>
    # </table>
    # """

    body!(W,"""
    <h3>Summary</h3>
    $eqs
    $summary
    $progress
    """),
    push!(M.summary.displays,W)
end


function displaysolve(M,iiter,Niter,MaxError,SError,t)
    if !active(M.summary.displays[end])
        initdisplay(M)
    else
        # @js M.summary.displays[end] document.getElementById("solvestep").innerHTML="Equations"
        @js M.summary.displays[end] document.getElementById("iiter").innerHTML=$iiter
        @js M.summary.displays[end] document.getElementById("Niter").innerHTML=$Niter
        @js M.summary.displays[end] document.getElementById("Telapsed").innerHTML=$t
        for i = 1:length(M.summary.equations.args)
            s = """$(signif(MaxError[i],2))  ($(signif(SError[i],2)))"""
            targ = "fval$i"
            @js M.summary.displays[end] document.getElementById($targ).innerHTML=$s
        end
    end
end

# function updatedisplay(M::Model,step,i,n,crit,t)
#     if !active(M.summary.displays[end])
#         initdisplay(M)
#     else
#         @js M.summary.displays[end] document.getElementById("solvestep").innerHTML=$step
#         iter = "$i/$n"
#         @js M.summary.displays[end] document.getElementById("iter").innerHTML=$iter
#         time = "$(signif(t,2))s"
#         @js M.summary.displays[end] document.getElementById("time").innerHTML=$time
#         for i = 1:length(M.summary.equations.args)
#             s = """$(signif(median(M.Fval[:,i]),2))  ($(signif(maximum(M.Fval[:,i]),2)))"""
#             targ = "fval$i"
#             @js M.summary.displays[end] document.getElementById($targ).innerHTML=$s
#         end
#     end
# end


#
# for i = 1:8
#     pn="progress$i"
# @js W document.getElementById($pn).innerHTML=$(i*10)
# end
#
# @js W document.getElementById("solvestep").innerHTML="a"
# @js M.summary.displays[end] document.getElementById("fval1").innerHTML
