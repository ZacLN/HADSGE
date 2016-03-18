function getE(x,list=Expr[])
    if isa(x,Expr)
        push!(list,x)
        for i = 1:length(x.args)
            list = getE(x.args[i],list)
        end
    end
    return list
end

function cntnode(x,cnt=0)
    if isa(x,Expr) && x.head!=:ref
        for i = 1:length(x.args)
            cnt += cntnode(x.args[i],cnt)
        end
    else
        return 1
    end
    return (cnt)
end

function exprreport(E::Expr,reps=2,nl = 100)
    es=getE(E)
    filter!(x->cntnode(x)>nl,es)
    h=[(sum(es.==e),e) for e in unique(es)]
    # cnts = [cntnode(e) for e in es]
    filter!(x->x[1]>reps,h)
    length(h)==0 && (return HTML("No expression repeated more than 5 times."))
    out = IOBuffer()
    for i = 1:length(h)
        println(out,Hiccup.div("#.boxed",string(h[i][2])))
        print(out,h[i][1])
    end
    s = takebuf_string(out)
    close(out)
    HTML(s)
end
