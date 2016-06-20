using Blink,WebSockets,HttpServer


function gui()
    wsh = WebSocketHandler() do req,client
        while true
            msg = read(client)
            println(String(msg))
        end
    end
    server = Server(wsh)
    hadsgewsid=rand(1:9999)
    @async run(server,hadsgewsid)

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
    body!(W,"""
    <form>
    <h2>Equations</h2>
    <textarea id="equations" style="color:white;background-color:#4d4d4d;width:80%;height:50%">
    equations1
    </textarea>
    </form>
    <button id="submit" onclick="hadsgews.send('00')">submit</button>
    """)
    @js W resizeTo(683,400)
    @js W moveTo(0,5)

    Blink.addscript!(W,"""
    hadsgews = new WebSocket("ws://localhost:$hadsgewsid")
    """)

    return W,server
end


w,s=gui()


Blink.addscript!(w,"""
function addEq()
{
    var eqs = document.getElementById('eqs');
    var eq = document.createElement("textarea");
    eq.id = "eq"+eqs.childElementCount;
    eqs.appendChild(eq);
}
function addRow()
{
    var eqs = document.getElementById('eqs').children[0];
    var eq = document.createElement("tr");
    var n = eqs.childElementCount+1;
    eq.innerHTML = '<td style="width:70%"><textarea id="eq'+n+'"style="width:100%">Eq.'+n+'</textarea>        </td>        <td style="width:5%">  <input size="5" type="text" name="polname'+n+'"><br></td>        <td style="width:2.5%"><input size="2" type="text" name="pollb'+n+'"><br></td>        <td style="width:2.5%"><input size="2" type="text" name="polub'+n+'"><br></td>        <td style="width:20%"><input size="10" type="text" name="polinit'+n+'"><br></td>'
    eqs.appendChild(eq);
}

function remLastEq()
{
    var eqs = document.getElementById('eqs');
    eqs.removeChild(eqs.lastChild);
}
""")

body!(w,"""
<h2>Equations </h2>
<table id="eqs" >
<tr>
    <td style="width:70%">
        <textarea id="eq1" style="width:100%">Eq.1</textarea>
    </td>
    <td style="width:5%">  <input size="5" type="text" name="polname1"><br></td>
    <td style="width:2.5%"><input size="2" type="text" name="pollb1"><br></td>
    <td style="width:2.5%"><input size="2" type="text" name="polub1"><br></td>
    <td style="width:20%"><input size="10" type="text" name="polinit1"><br></td>
</tr>
</table>

<button id="addeq" onclick="addEq()">add equation</button>
<button id="addr" onclick="addRow()">+</button>
<button id="addeq" onclick="remLastEq()">X</button>
<button id="submit" onclick="hadsgews.send('00')">submit</button>
""")

@js w document.getElementById("equations").value



if false
    close(w)
    close(s)
end
