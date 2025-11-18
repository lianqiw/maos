"use strict"; //every variable need to be defined
function lpfDrawdata(newData, oldData, lpf){
  try{
    if('data' in oldData && 'data' in newData && lpf!=1){
      const ny=newData['data'].length;
      if(ny>0 && ny==oldData['data'].length){
        const nx=newData['data'][0].length;
        if(nx==oldData['data'][0].length){
          for (let iy = 0; iy < ny; iy++) {
            for (let ix=0; ix<nx; ix++){
              newData['data'][iy][ix]=newData['data'][iy][ix]*lpf+oldData['data'][iy][ix]*(1.-lpf);
            }
          }
        }
      }
    }/*else if('pts' in oldData && 'pts' in newData){
      const nx=newData['pts'].length;
      if(nx==oldData['pts'].length){
        for (let ix=0; ix<nx; ix++){
          newData['pts'][ix]=newData['pts'][ix]*lpf+oldData['pts'][ix]*(1-lpf);
        }
      }
    }*/
  }catch(err){
    console.log(now(), err);
  }
}
const DrawDaemon = React.memo(({ drawInfo, jobActive, updateDrawInfo}) => {
  //Notice that the state variables do CHANGE at every re-rendering. It just does not change within a rendering
  const { useState, useEffect, useRef } = React;
  const [topActive, setTopActive] = useState({});//active fig
  const [botActive, setBotActive] = useState({});//active name
  const [figSeq, setFigSeq] = useState(0);//figure sequence number (index of drawData). triggers plotly update.
  const [cumStart, setCumStart]=useState(0);//time step to start cumulative plot.
  const [cumInput, setCumInput]=useState(0);//capture input context
  const [lpfInput, setLpfInput]=useState(1);//low pass filter to slow down update
  const lpf=useRef(1); //use reference avoid closure capture stale data
  const figName=useRef('');
  const [cumPlot, setCumPlot]=useState(false);//plot cumulative plot
  const [wss, setWss]=useState({});
  const [pause, setPause]=useState({});
  //Use useRef for data that is not directly related to rendering.
  //const pause = useRef({});//pause plotting
  const jobRef = useRef({});//job data. 
  const wssRef = useRef({});//mutable ref object that persists for the entire lifecycle of the component.
  const chartRef = useRef(null);//for plotly. References a DOM object.
  const cumInputRef=useRef(null);//for input of cumStart
  const lpfInputRef=useRef(null);//for input of lpf
  const layout={autosize: true, margin: {t:40, b:40, l:50, r:20}};
  function connect(hostname, job) {
    if (!hostname || hostname.length == 0) return false;
    const [host, pid] = job.split(':', 2);
    if (hostname.indexOf(':') == -1) {
      hostname += port;
    }
    let ws;
    if (!wssRef.current[job]) {
      //console.log(now(),`Drawdaemon is connecting for ${job} at ${pcol}${hostname}`);
      try {
        ws = new WebSocket(pcol + hostname + "/xxx", "maos-drawdaemon-protocol");	/* + "/xxx" bit is for IE10 workaround */
      } catch (error) {
        console.log(now(),error);
      }
    } else {
      console.log(now(),`Drawdaemon is already active for ${host} ${pid}`);
    }
    if (!ws) return false;
    ws.onopen = () => {
      ws.binaryType = 'arraybuffer';
      setWss((oldVal) => ({ ...oldVal, [job]: true }));
      wssRef.current[job] = ws;
      setBotActive((oldVal) => ({ ...oldVal, [job]: {} }))//Initialize botActive to empty object to avoid undefined error
      setTopActive((oldVal) => ({ ...oldVal, [job]: {} }))//initialize to empty
      jobRef.current[job] = { 'drawData': {} };//empty plots
      console.log(now(),`Drawdaemon connected for ${job}`);
      ws.send(pid + "&" + "DRAW" + ";");
    };

    ws.onclose = (event) => {//we keep plots when connection is closed.
      setWss((oldVal) => ({ ...oldVal, [job]: false }));
      delete wssRef.current[job];
      console.log(now(), `Drawdaemon disconnected from ${job}`);
    };
    ws.onerror = (err) => {
      console.error(now(), hostname, "Drawdaemon connection error", err);
    };
    ws.onmessage = (event) => {
      if (typeof event.data === "string") {//instanceof does not work for string
        console.log(now(), `Got unexpected text data ${event.data}`);
      } else if (event.data instanceof Blob) {//Blob is read only
        console.log(now(),"Got unexpected blob data with bytes ", event.data.size);
      } else if (event.data instanceof ArrayBuffer) {//ArrayBuffer support writing.
        let drawData = procBuffer(event.data); //separate into function for testing
        if ('pid' in drawData) {
          jobRef.current[job]['pid'] = drawData['pid'];
          updateDrawInfo([job, host+':'+drawData['pid']])
        }
        if ('path' in drawData) {
          jobRef.current[job]['path'] = drawData['path'];
        }
        if ('exename' in drawData) {
          jobRef.current[job]['exename'] = drawData['exename']
        }
        if ('final' in drawData) {
          jobRef.current[job]['pid'] = -1;//draw is no longer active
        }
        const fig = drawData['fig'];
        const name = drawData['name'];
        if (fig && name) {
          try {
            if (!jobRef.current[job]['active']) {
              jobRef.current[job]['active'] = fig;
              //Initialize topActive to first figure 
              setTopActive((oldVal) => ({ ...oldVal, [job]: fig }))
            }
            if (!jobRef.current[job]['drawData'][fig]) {
              jobRef.current[job]['drawData'][fig] = {}
              //Initialize botActive to first figure of the group
              setBotActive((oldVal) => ({ ...oldVal, [job]: { ...oldVal[job], [fig]: name } }))
            }
            if (lpf.current!=1 && name in jobRef.current[job]['drawData'][fig]){
              lpfDrawdata(drawData, jobRef.current[job]['drawData'][fig][name], lpf.current);
            }
            jobRef.current[job]['drawData'][fig][name] = drawData;
          } catch (err) {
            console.log(now(),{ err, drawData, jobRef, job });
          }
        }
        setFigSeq(old => old + 1);
      } else {
        console.error(now(), "Invalid data:", event.data);
      }
    };

    return true;
  }//function connect

  useEffect(() => {//run ones upon mount and state change
    Object.keys(drawInfo).map((job) => {
      if (drawInfo[job]) {//host:port->hostname
        if (!wssRef.current[job]) {
          const pid = job.split(':').at(-1);
          if(pid<=0 || !jobRef.current[job]){
            connect(drawInfo[job].hostname, job);
          }//else: do not auto-reconnect existing job plots
        }
      } else {
        if (wssRef.current[job]) {//close page
          console.warn(`Close Drawdaemon connection to ${job}`)
          wssRef.current[job].close();//drop connection
        }
        jobRef.current[job] = { 'drawData': {} };//remove plots
        Plotly.purge(chartRef.current, 0);//delete plot
      }
    })
  }, [drawInfo]);
  
  useEffect(() => {//trigger plotting if page switched or data is updated
    if (chartRef.current && !pause[jobActive]) {
      try {
        let clear = 1;
        if (jobActive !== '' && jobRef.current[jobActive] && jobRef.current[jobActive]['drawData']
          && jobRef.current[jobActive]['drawData'][topActive[jobActive]]
          && botActive[jobActive][topActive[jobActive]]) {
          const drawData = (jobRef.current[jobActive]['drawData'][topActive[jobActive]][botActive[jobActive][topActive[jobActive]]]);
          if (drawData) {
            const id=drawData.fig+'_'+drawData.name;
            if(figName.current!==id){
              figName.current=id;
            }
            layout.uirevision=figName.current;//preserve zoom for the plot
            drawData.cumStart=cumStart;
            drawData.cumPlot=cumPlot;
            const traces  = makeTraces(drawData, layout);
            if (traces.length) {
              Plotly.react(chartRef.current, traces, layout, { responsive: true });
              clear = 0;
            }
          }
        }
        if (clear) {
          Plotly.purge(chartRef.current, 0);
        }
      } catch (err) {
        Plotly.purge(chartRef.current, 0);
        console.log(now(),err, { jobActive, jobRef, topActive, botActive });
      }
    }
    //Use ResizeObserver to handle container resize
    const resizeObserver = new ResizeObserver(() => {
      Plotly.Plots.resize(chartRef.current);
    });
    resizeObserver.observe(chartRef.current);

    return () => resizeObserver.disconnect();
  }, [jobActive, topActive, botActive, figSeq, cumStart, cumPlot, pause]);

  useEffect(() => {//let server know the current jobActive page.
    try {
      if (wssRef.current[jobActive] && jobRef.current[jobActive] && jobRef.current[jobActive]['pid'] > -1
        && topActive[jobActive] && botActive[jobActive][topActive[jobActive]]) {
        wssRef.current[jobActive].send(`${jobRef.current[jobActive]['pid']}&DRAW_FIGFN&${topActive[jobActive]}&${botActive[jobActive][topActive[jobActive]]}`);
      }
    } catch (err) {
      console.log(now(),{ err, jobActive, wssRef, jobRef, topActive, botActive });
    }
  }, [jobActive, topActive, botActive]);
  try {
    return jobRef.current[jobActive] && jobRef.current[jobActive]['drawData'] && Object.keys(jobRef.current[jobActive]['drawData']).length ? (<div>
      <ul className="inline tab_hosts"> {/*draw top-notebook (horizontal)*/}
        {Object.keys(jobRef.current[jobActive]['drawData']).sort().map((fig) => (
          <li key={fig} className={topActive[jobActive] === fig ? "active" : ""}
            onClick={() => { setTopActive(oldVal => ({ ...oldVal, [jobActive]: fig })); setPause(oldVal=>({...oldVal, [jobActive]:false}));}}
          >{fig}</li>))}  
        <li title="Pause or Resume ploting" onClick={() => { setPause(oldVal=>({...oldVal, [jobActive]:oldVal[jobActive]?false:true}));}}>{pause[jobActive] ? "▶️" : "⏸️"}</li>
        <li title="Stop receiving more data for plotting" onClick={() => { if(wssRef.current[jobActive]) wssRef.current[jobActive].close() }}>{wss[jobActive]?"⏹️" : "🔴"}</li>
        <li title="Cumulative ploting" className={cumPlot?"active":""} onClick={()=>{if(!cumInput) {setCumInput(0.1); setCumStart(0.1);}setCumPlot(oldVal=>!oldVal)}}>🎢</li>
        <Menu label={(<span style={{padding:'0.3em'}}>Menu</span>)} child={<ul className="menu-list" >
          <li title="Set cumulative plotting starting index"><span>Cum Start</span><div className="spring-spacer"></div>
            <form onSubmit={(e)=>{e.preventDefault(); setCumStart(parseFloat(cumInput.length?cumInput:"0")); setCumPlot(true);}}>
            <input ref={cumInputRef} style={{width:'4em'}} value={cumInput} type="number" step="0.1" min="0"
            onClick={(e)=>{if(cumInputRef.current) cumInputRef.current.select(); e.stopPropagation();}} 
            onChange={e=>{setCumInput(e.target.value);}}></input></form></li>
          <li title="Set image update low pass filter"><span>Update LPF </span><div className="spring-spacer"></div>
            <form onSubmit={(e)=>{e.preventDefault(); lpf.current=(parseFloat(lpfInput.length?lpfInput:"0.01"));}}>
            <input ref={lpfInputRef} style={{width:'4em'}} value={lpfInput} type="number" step="0.001" min="0.001" max="1"
            onClick={()=>{if(lpfInputRef.current) lpfInputRef.current.select();}} 
            onChange={e=>{setLpfInput(e.target.value);}}></input></form></li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'png',filename: figName.current});}}>💾 Save as {figName.current}.png</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'jpeg',filename: figName.current});}}>💾 Save as {figName.current}.jpg</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'svg',filename: figName.current});}}>💾 Save as {figName.current}.svg</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'webp',filename: figName.current});}}>💾 Save as {figName.current}.webp</li>
        </ul>}></Menu>
      </ul>
      <div className="layout">
        <ul className="tab_hosts" style={{width:'8em'}}>{/*draw sub-notebook (vertical)*/}
          {jobRef.current[jobActive]['drawData'][topActive[jobActive]] &&
            Object.keys(jobRef.current[jobActive]['drawData'][topActive[jobActive]]).sort().map((name) => (
              <li key={name} className={botActive[jobActive][topActive[jobActive]] === name ? "active" : ""}
                onClick={() => { setBotActive(oldVal => ({ ...oldVal, [jobActive]: { ...oldVal[jobActive], [topActive[jobActive]]: name } })); }}
              >{name}</li>))}
        </ul>
        <div id="chart" ref={chartRef}></div></div>
    </div>) : (<div id="chart" ref={chartRef}></div>)
  } catch (err) {
    console.log(now(),{ err, jobRef, topActive, botActive, jobActive });
  }
})

window.DrawDaemon = DrawDaemon;
