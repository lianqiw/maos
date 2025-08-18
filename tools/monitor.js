"use strict"; //every variable need to be defined
var port=':80';
var pcol="ws://";
const wssGlobal={};//for checking connection. state is not dependable.
const drawHost={};
const iconName={1:"▶️", 2:"🕒", 3:"▶️", 4: "🕒", 11:"✅",12:"❌",13:"❌",14:"❌",15:"❌"};
const draw={
  start:0,
  data:1,
  heartbeat:2,
  points:3,
  style:4,
  circle:5,
  limit:6,
  fig:7,
  name:8,
  title:9,
  xlabel:10,
  ylabel:11,
  zlim:12,
  legend:13,
  xylog:14,
  figfn:15,
  pause:16,
  resume:17,
  final:18,
  float:19,
  frame:20,
  single:21,
  udpport:22,
  init:23,
  pid:24,
  zlog:25,
  path:26,
  exename:27,
  end:100,
  entry:9999
};
var byteFloat=4;

function get_hostname() {
	/*
	* We open the websocket encrypted if this page came on an
	* https:// url itself, otherwise unencrypted
	*/
	console.log(window.location.href);
	var u = window.location.href;
	if (u.substring(0, 5) === "https") {
		pcol = "wss://";
		u = u.substring(8);
	} else {
		if (u.substring(0, 4) === "http"){
			u = u.substring(7);
		}
	}
  u = u.split('/')[0];
  const p=u.indexOf(':');
  if(p!=-1){
    port=u.substring(p);
  }
	return u;
}
function split_hostname(hostname){
  let host=hostname.split('://').pop().split(':')[0];
  if (/^(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)$/.test(host)) {  
    return host;
  }else{
    return host.split('.')[0];
  }
}
function linspace(start, stop, length) {
    if (stop === undefined) {
      stop = start;
      start = 0;
    }
    if(length===undefined) length=Math.round(stop-start);
    const step=(length==1)?1:(stop-start)/(length-1);
    const arr = Array(length);
    for (let i = 0; i < length; i++) {
      arr[i] = start + i * step;
    }

    return arr;
}
function arange(start, stop, step) {
    if (stop === undefined) {
      stop = start;
      start = 0;
    }
    if (step === undefined) {
      step = 1;
    }

    const length = Math.max(Math.ceil((stop - start) / step), 0);
    const arr = Array(length);

    for (let i = 0; i < length; i++) {
      arr[i] = start + i * step;
    }

    return arr;
}
function mkRGB(val){
  const r=(val>>16)&0xFF;
  const g=(val>>8 )&0xFF;
  const b=val&0xFF;
  return `rgb(${r},${g},${b})`;
}
function mkCircle(cir){
  const cx=cir[0];
  const cy=cir[1];
  const radius=cir[2];
  const color=cir[3];
  const np=128;
  const x=new Float32Array(np+1);
  const y=new Float32Array(np+1);
  const dt=Math.PI*2./np;
  for(let i=0; i<np; i++){
    const angle=dt*i;
    x[i]=Math.cos(angle)*radius+cx;
    y[i]=Math.sin(angle)*radius+cy;
  }
  x[np]=x[0];
  y[np]=y[0];
  return {x,y,mode:'lines',color:mkRGB(color),showlegend:false};
}
function mkCircles(cirs){
  const ncir=cirs.length>>2;
  const traces=[];
  for(let i=0; i<ncir; i++){
    traces.push(mkCircle(cirs.slice(i*4,i*4+4)));
  }
  return traces;
}
function reshapeTo2D(p, nx, ny){
  const res=[]
  for(let iy=0; iy<ny; iy++){
    res.push(p.slice(iy*nx, (iy+1)*nx));
  }
  //console.log(`reshapeTo2D nx=${nx} ny=${ny} res[0].length=${res[0].length} res.length=${res.length}`)
  return res;
}
function App() {
  const {useState, useEffect, useRef}=React;
  //Notice: useState is async. It does NOT update the variable in the current rendering cycle
  //When the state update depends on the previous state value, it is crucial to use the functional update form of the setting function.
  //closure captures state variable might retain the value it had when the effect was initially created.
  //useState variable should only be accessed in set or useEffect with state listed in [].
  const [data, setData] = useState([]);//save received jobs
  const [wss, setWss] = useState({});//save websocket information
  const [active, setActive]=useState('');//current active host
  const [text, setText]=useState('');//text for input
  const [topnb, setTopnb]=useState([]);//for fig
  const [botnb, setBotnb]=useState([]);//for fig name
  const [topActive, setTopActive]=useState({});//active fig
  const [botActive, setBotActive]=useState({});//active name
  const [drawActive, setDrawActive]=useState({});//active drawData
  const [figCount, setFigCount]=useState(0);//total number of figures
  const [figSeq, setFigSeq]=useState(0);//figure sequence number (index of drawData). triggers plotly update.
  const chartRef = useRef(null);//for plotly. References a DOM object.

    //useEffect(() => {//Effect function runs after React updates the DOM.
  function connect(hostname){
    if(hostname.length==0) return false;
    const host=split_hostname(hostname);//host:port->shortname
    if(hostname.indexOf(':')==-1){
      hostname+=port;
    }
    let ws;
    if(!wssGlobal[host]){
      try{
        ws = new WebSocket(pcol+hostname+"/xxx", "maos-monitor-protocol");	/* + "/xxx" bit is for IE10 workaround */
        ws.binaryType = 'arraybuffer';
      }catch(error){
        console.log(error);
      }
    }
    if(!ws) return false;

    ws.onopen = () => {
      wssGlobal[host]=ws;
      console.log(wssGlobal);
      setWss((oldVal)=>{
        return {...oldVal, [host]:ws};});
      setActive((oldVal)=>oldVal?oldVal:host);
      console.log(host, "WebSocket connected", ws);
      if(localStorage.hosts.indexOf(hostname)==-1){
        localStorage.hosts+=";"+hostname;
      }
      resetdrawHost(host);
    };

    ws.onmessage = (event) => {
      var j;
      if(typeof event.data=== "string") {//instanceof does not work for string
        if (event.data.indexOf('$')!=-1){
          j = event.data.split('$'); //newer format to avoid conflict with actual data.
        }else{
          j = event.data.split(';');
        }
        if(j.length>0){
          for(let ij=0; ij<j.length-1; ij++){
            const i=j[ij].split('&');
            var newdata;
            if(i.length==3){
              let path="";
              const ip=i[2].indexOf(" ");
              if(ip!=-1){
                path=i[2].substring(0,ip);
                i[2]=i[2].substring(ip);
              }
              const io=i[2].lastIndexOf("-O");
              let dirout="";
              if(io!=-1){
                dirout=i[2].substring(io+2).trim().split(' ')[0];
                i[2]=i[2].replace(/[ \t]+-o[ \t]+[^ \t]+/g, " ").trim();
              }
              newdata={PID:i[0], Host:host, Path:path, "Dir Out":dirout, Args:i[2], status:0};
            }else if(i.length==14){
              if(i[3]!=14){//14: remove
                i[11]=parseInt(i[11]);
                i[12]=parseInt(i[12]);
                let icon="⏩";
                if(i[12]>0 || i[3]!=11){
                  icon=iconName[i[3]];
                }
                newdata={PID:i[2], Host:host, status:i[3], Time:i[4], High: i[5], Low:i[6], iseed:i[7], nseed:i[8], isim:i[9], nsim:i[10], 
                  rest:i[11], tot:i[12], Step:i[13], icon:icon};
              }
            }else{
              continue;//invalid data
            }
            setData((oldValue)=>{
              const ind=oldValue.findIndex((v=>v.PID==i[0]&&v.Host===host));
              if(newdata===undefined && ind!=-1){//remove
                return oldValue.filter((v,i)=>i!=ind);
              }else if(ind==-1){//append
                return [...oldValue, newdata];
              }else{//update
                return oldValue.map((v,i)=>(i===ind)?{...v, ...newdata}:v);
              }
            });
          }
        }
      }else if(event.data instanceof Blob){//Blob is read only
        console.log("Got blob data with bytes ", event.data.size);
      }else if(event.data instanceof ArrayBuffer){//ArrayBuffer support writing.
        procBuffer(event.data); //separate into function for testing
      }else{
        console.error("Invalid data:", event.data);
      }
    };

    ws.onclose = () => {
      wssGlobal[host]=false;
      console.log(wssGlobal);
      setWss((oldVal)=>({...oldVal, [host]:false}));
      setData((oldVal)=>oldVal.filter((v)=>v.Host!=host));//remove jobs of the host
      console.warn(hostname,"disconnected");
      resetdrawHost(host)
    };

    ws.onerror = (err) => {
      console.error(hostname,"WebSocket error", err);
    };
    return true;
  
    function procBuffer(wsBuf){
      //console.log("Got ArrayBuffer data with bytes ", wsBuf.byteLength);
      const dataView = new DataView(wsBuf);
      var buf={pos:0,cmd:0};
      var drawData={};
      function getKey(){
        Object.keys(draw).find(key=>draw[key]===buf.cmd)
      }
      function getInt(){
        let v=dataView.getInt32(buf.pos, true); buf.pos+=4;
        return v;
      }
      function getByteArray(e=-1){
        const nx=e>-1?e:getInt();
        let v;
        try{
          v=new Uint8Array(wsBuf, buf.pos, nx); buf.pos+=nx;
        }catch(err){
          console.log(err, getKey(), wsBuf.byteLength, buf.pos, nx, e);
        }
        return v;
      }
      function getIntArray(e=-1){
        const nx=e>-1?e:getInt();
        const newBuffer=wsBuf.slice(buf.pos, buf.pos+nx*4); buf.pos+=nx*4;
        const v=new Int32Array(newBuffer); 
        return v;
      }
      function getFloatArray(e=-1){
        const nx=e>-1?e:getInt();
        const newBuffer=wsBuf.slice(buf.pos, buf.pos+nx*byteFloat);buf.pos+=nx*byteFloat;
        let v;
        try{
          if(byteFloat==4){//use slice to avoid byteoffset alignment error.
            v=new Float32Array(newBuffer);
          }else if(byteFloat==8){
            v=new Float64Array(newBuffer); 
          }
        }catch(error){
          console.log(getKey(buf.cmd), "newFloatArray fails with", error, 
            newBuffer, newBuffer.byteLength, buf.pos, nx, byteFloat, wsBuf.byteLength);
        }
        return v;
      }
      function getStr(){
        return String.fromCharCode.apply(null, getByteArray().slice(0,-1));//drop terminating \0
      }
      let len;
      while(buf.pos<wsBuf.byteLength){
        buf.cmd=getInt()
        switch(buf.cmd){
        case draw.entry:
          len=getInt();//length
          continue;
          break;
        case draw.frame:
          buf.pos+=16;//skip 4 ints
          break;
        case draw.start:
          break;
        case draw.data:{
          let nx=getInt();
          let ny=getInt();
          drawData['data']=reshapeTo2D(getFloatArray(nx*ny), nx, ny);
          }break;
        case draw.heartbeat:
          break;
        case draw.points:{
          let nx=getInt();
          let ny=getInt();
          let square=getInt();
          if(!drawData['pts']){
            drawData['pts']=[];
          }
          drawData['pts'].push({p:getFloatArray(nx*ny),nx:nx,ny:ny,square:square});
        }break;
        case draw.style:
          drawData['style']=getIntArray();
        break;
        case draw.circle:{
          const nc=getInt()
          drawData['circle']=getFloatArray(4*nc);
        }break;
        case draw.limit:
          drawData['limit']=getFloatArray(4);
          drawData['limitManual']=1;
          break;
        case draw.fig:
          drawData['fig']=getStr();
          break;
        case draw.name:
          drawData['name']=getStr();
          break;
        case draw.title:
          drawData['title']=getStr();
          break;
        case draw.xlabel:
          drawData['xlabel']=getStr();
          break;
        case draw.ylabel:
          drawData['ylabel']=getStr();
          break;
        case draw.zlim:
          drawData['zlim']=getFloatArray(2);
          drawData['zlimManual']=1;
          drawData['zlogLast']=0;
          break;
        case draw.legend:
          drawData['legend']=[]
          for(let il=0; il<drawData['pts'].length; il++){
            drawData['legend'].push(getStr());
          }
          break;
        case draw.xylog:
          drawData['xylog']=getByteArray(2);
          break;
        case draw.final:
          //console.log("draw.final");
          break;
        case draw.float:
          byteFloat=getInt();
          //console.log("byteFloat is", byteFloat);
          break;
        case draw.udpport:
          getInt();//ignored
          console.log("UDP is not supported");
          break;
        case draw.init:
          //console.log("draw.init");
          break;
        case draw.pid:
          drawHost[host]['pid']=getInt();
          console.log("drawHost[host] is", drawHost[host])
          break;
        case draw.zlog:
          drawData['zlog']=getInt();
          break;
        case draw.path:
          drawHost[host]['path']=getStr();
          console.log("drawHost[host] is", drawHost[host])
          break;
        case draw.exename:
          drawHost[host]['exename']=getStr();
          console.log("drawHost[host] is", drawHost[host])
          break;
        case draw.end:
          //console.log(`drawData (${wsBuf.byteLength} bytes):`, drawData)
          if(!drawHost[host]['drawData'][drawData['fig']]){
            drawHost[host]['drawData'][drawData['fig']]={}
            setFigCount(old=>old+1);
          }
          drawHost[host]['drawData'][drawData['fig']][drawData['name']]=drawData;
          setFigSeq(old=>old+1);
          //start making plot
          break;
        default:
          console.warn("unknown command", buf.cmd);
          buf.pos+=len;//skip its payload
      }//switch
      }//while
    }//function procbuffer
  }//function connect
  useEffect(()=>{//run once upon mount
    const hostname=get_hostname();//connect current host
    if(!localStorage.hosts){
      localStorage.hosts=hostname;
    }else if(localStorage.hosts.indexOf(hostname)==-1){
      localStorage.hosts+=";"+hostname;
    }
    console.log(localStorage.hosts);
    const hosts=[...new Set(localStorage.hosts.split(";"))];//unique entries
    localStorage.hosts=hosts.filter((hostname)=>connect(hostname)).join(";");//keep only valid hosts
  },[]);
  function resetdrawHost(host){
      if(drawHost[host]){
        Object.keys(drawHost[host]).forEach(key=>delete drawHost[host][key]);
        drawHost[host]['drawData']={};
      }else{
        drawHost[host]={'drawData':{}};
      }
      setTopnb([]);
  }
  useEffect(()=>{
    if(drawHost[active]&&drawHost[active]['drawData']){
      setTopnb(Object.keys(drawHost[active]['drawData']).sort());
    }else{
      setTopnb([]);
    }
  },[active,figCount]);//number of figure changes

  useEffect(()=>{
    try{
      if(topnb.length>0){
        if(topActive[active]==='' || topnb.indexOf(topActive[active])==-1){
          setTopActive(oldval=>({...oldval, [active]:topnb[0]}));
        }
        if(drawHost[active]&&drawHost[active]['drawData'] && topActive[active]!==''){
          setBotnb(Object.keys(drawHost[active]['drawData'][topActive[active]]));//length may change
        }
      }else{
        setTopActive(oldval=>({...oldval, [active]:''}));
      }
    }catch(err){
      console.log(err, active, drawHost[active], topnb, botnb, topActive, botActive);
    }
  },[topnb]);//executes whenever variables changes

  useEffect(()=>{
    try{
      if(drawHost[active]&&drawHost[active]['drawData'] && topActive[active]!==''){
        setBotnb(Object.keys(drawHost[active]['drawData'][topActive[active]]).sort());//length may change
      }else{
        setBotnb([]);
      }
    }catch(err){
      console.log(err, active, drawHost[active], topnb, botnb, topActive, botActive);
    }
  },[topActive]);

  useEffect(()=>{
    try{
      if(botnb.length>0){
        if(botActive[active]==='' || botnb.indexOf(botActive[active])==-1){
          setBotActive(oldval=>({...oldval, [active]:botnb[0]}))
        }
      }else{
        setBotActive(oldval=>({...oldval, [active]:''}))
      }
    }catch(err){
      console.log(err, active, drawHost[active], topnb, botnb, topActive[active], botActive);
    }
  },[botnb]);

  useEffect(()=>{//trigger plotting if page switched or data is updated
    try{
      if(drawHost[active]&&drawHost[active]['drawData'] && topActive[active]!=='' && botActive[active]!==''){
        setDrawActive(drawHost[active]['drawData'][topActive[active]][botActive[active]]);
      }else{
        Plotly.purge(chartRef.current, 0);
      }
    }catch(err){
      console.log(err, active, drawHost[active], topnb, botnb, topActive[active], botActive);
    }
  },[topActive,botActive,figSeq]);

  useEffect(()=>{//let server know the current active page
    try{
      if(drawHost[active]&&drawHost[active]['drawData'] && topActive[active]!=='' && botActive[active]!==''){
        cmdHostPid(active, drawHost[active]['pid'], `DRAW_FIGFN&${topActive[active]}&${botActive[active]}`);
      }
    }catch(err){
      console.log(err, active, drawHost[active], topnb, botnb, topActive[active], botActive);
    }
  },[topActive,botActive]);

  function cmdHostPid(host, PID, cmd){
    if(wssGlobal[host]){
      wssGlobal[host].send(PID+"&"+cmd+";");
    }
  }
  function cmdHost(cmd){
    if(cmd==="kill_all" && !window.confirm("Kill all jobs on " + active +"?")){
      return;
    }
    console.log(cmd, active)
    data.map((v)=>{
      if(active.length==0 || v.Host===active){
        if(cmd==="clear_finished"){
          if(v.status==11){
            cmdHostPid(v.Host, v.PID, "REMOVE");
          }
        }else if(cmd==="clear_skipped"){
          if(v.status==11 && v.tot==0){
            cmdHostPid(v.Host, v.PID, "REMOVE");
          }
        }else if(cmd==="clear_crashed"){
          if(v.status>11){
            cmdHostPid(v.Host, v.PID, "REMOVE");
          }
        }else if(cmd==="kill_all"){
          if(v.status<11){
            cmdHostPid(v.Host, v.PID, "KILL");
          }
        }
      }
    })
  }
  useEffect(()=>{//whenever drawActive changes
    const drawData=drawActive;
    if(!drawData){
      return;
    }
    
    let square=0;
    let ratio=1;
    const layout={title:drawData['title'], width:800, height:600, 
            xaxis:{showgrid:true, showline:true,mirror:true,ticks:"inside",zeroline:false},
            yaxis:{showgrid:true, showline:true,mirror:true,ticks:"inside",zeroline:false},
          }
    let traces=[];
    try{
      if(drawData['data']){
        square=1;
        const nx=drawData['data'][0].length;
        const ny=drawData['data'].length;
        ratio=(nx&&ny)?(ny/nx):1;
        let x,y;
        if(drawData['limit']){
          x=linspace(drawData['limit'][0],drawData['limit'][1],nx)
          y=linspace(drawData['limit'][2],drawData['limit'][3],ny)
        }
        traces.push({x, y, z:drawData['data'], type:'heatmap', colorscale:'Jet', showlegend:false});
      }else if(drawData['pts']){
        const showlegend=drawData['pts'].length>1?1:0;
        traces=traces.concat(drawData['pts'].map((pts, ipts)=>{
        if(pts.ny==2){
          if(pts.square) square=1;
          return {x:pts.p.slice(0,pts.nx), y:pts.p.slice(pts.nx,2*pts.nx), 
            mode:"markers", marker:{symbol:'x'},  name:drawData['legend']?drawData['legend'][ipts]:"", showlegend}
        }else{
          return {y:pts.p, mode:"lines+markers", name:drawData['legend']?drawData['legend'][ipts]:"", showlegend}
        }}));
      }
      if(drawData['circle']){
        traces=traces.concat(mkCircles(drawData['circle']));       
      }
      if(square){
        layout.yaxis.scaleanchor='x';
        layout.yaxis.scaleratio=ratio;
        layout.height=layout.width*ratio;
      }
      if(traces.length){
          Plotly.newPlot(chartRef.current, traces, layout, {responsive:true});
      }else{
        Plotly.purge(chartRef.current, 0);
      }
    }catch(err){
      console.log(err, drawData, traces, layout);
    }
  },[drawActive])
  function Progress({text, frac}){
    return(
      <td>
      <div className="progress" style={{"width":"160px"}} key="progress">
        <div className="progressbar" style={{"width":frac}}></div>
        <span className="progressbar">{text}</span>
      </div>
      </td>
    );
  }
  const columns=["Time", "Host", "PID", "Path", "Dir Out", "Low", "High", "Step"];
  const cn={Time:"",Host:"",PID:"",Path:"tdpath","Dir Out":"tdout",Low:"",High:"",Step:""};
  return (
    <div>
      <div className="inline">
        <img src="icon-monitor.png" alt="icon"></img>
        <ul className="tab_hosts">
        <li onClick={()=>{cmdHost("clear_finished")}}>Clear Finished</li>
        <li onClick={()=>{cmdHost("clear_skipped")}}>Clear Skipped</li>
        <li onClick={()=>{cmdHost("clear_crashed")}}>Clear Crashed</li>
        <li onClick={()=>{cmdHost("kill_all")}}>Kill All</li>
        <li><input value={text} onChange={e=>setText(e.target.value)}></input>
          <button onClick={()=>{setText(''); if(text.length){connect(text);}}}>Connect</button>
        </li> 
        {Object.keys(wss).map((host)=>(
        <li key={host} className={active===host?"active":""}>
          <span onClick={()=>{if(!wss[host])connect(host); setActive(host); setFigSeq(old=>old+1)}}>{(wss[host]?"🟢":"🔴")+host}</span>
          <span onClick={()=>{
            if(wss[host]){
              wss[host].close();
            }else{//already closed, remove
              setWss((oldVal)=>{const {[host]:_, ...rest } = oldVal; return oldVal;})
              localStorage.hosts.replace(host,""); 
            }}}>✖️</span>
        </li>
      ))}
      </ul></div>
      {columns.length > 0 ? (
        <table className="monitor">
          <thead>
            <tr>{columns.map((col) => <th key={col} className={cn[col]}>{col}</th>)}
            <th key="progress">Progress</th>
            <th key="icon"></th>
            </tr>
          </thead>
          <tbody>
            {data.filter((row)=>row.status>0 && (active.length==0 || row.Host===active)).map((row, i) => (
              <tr key={row.PID}>
                {columns.map((col) => <td key={col} className={cn[col]}>{row[col]}</td>)}
                <Progress
                  text={row.tot==0?"":row.iseed+'/'+row.nseed+' '+row.isim+'/'+row.nsim+' '+row.rest+'/'+row.tot}
                  frac={Math.round(100*(1-(row.tot==0?1:row.rest/row.tot)))+'%'}>
                </Progress>                  
                <td onClick={()=>{
                    if(row.status < 11){
                      if(window.confirm("Kill the job " + row.PID +"?")){
                        cmdHostPid(row.Host, row.PID, "KILL");
                      }
                    }else{
                      cmdHostPid(row.Host, row.PID, "REMOVE");
                    }
                  }
                }onContextMenu={(event)=>{
                    event.preventDefault();
                    if(row.status < 11){
                      //if(window.confirm("Plot the job " + row.PID +"?")){
                      if(!drawHost[row.Host] || !drawHost[row.Host]['pid']){
                        cmdHostPid(row.Host, row.PID, "DRAW");
                      }else if(drawHost[row.Host]['pid'] !== row.PID){
                        cmdHostPid(row.Host, row.PID, "DRAW");
                      }else{
                        console.log("Already plotting\n");
                      }
                    }
                  }
                }>{row.icon}</td>
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        <p>No data yet</p>
      )}
      <ul className="inline tab_hosts">
        {topnb.map((fig)=>(<li key={fig} className={topActive[active]===fig?"active":""}onClick={()=>{setTopActive(oldval=>({...oldval, [active]:fig}));}}
        >{fig}</li>))}
      </ul>
      <table border="0"><tbody><tr valign="top"><td>
      <ul className="tab_hosts">
        {botnb.map((fn)=>(<li key={fn} className={botActive[active]===fn?"active":""}onClick={()=>{setBotActive(oldval=>({...oldval, [active]:fn}));}}
        >{fn}</li>))}
      </ul></td><td>
      <div id="chart" ref={chartRef}></div>
      </td></tr></tbody></table></div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
