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
  const [cumInput, setCumInput]=useState(0);//capture input value without changing cumStart instantenously
  const [lpfInput, setLpfInput]=useState(1);//low pass filter to slow down update
  const [climLpfInput, setClimLpfInput]=useState(0.01)
  const lpf=useRef(1); //use reference avoid closure capture stale data
  const figName=useRef('');
  const [cumPlot, setCumPlot]=useState(false);//plot cumulative plot
  const [xRange, setXRange]=useState(0);//x axis range. Only applicable when the x axis lable contains Step
  const [wss, setWss]=useState({});
  const [pause, setPause]=useState({});
  //Use useRef for data that is not directly related to rendering.
  const jobRef = useRef({});//job data. 
  const wssRef = useRef({});//mutable ref object that persists for the entire lifecycle of the component.
  const chartRef = useRef(null);//for plotly. References a DOM object.
  const layout=useRef({});//remember the current layout
  const [chartLoaded, setChartLoaded]=useState(0)
  const gif=useRef(null)
  const gifActive=useRef(0)
  const chartSize=useRef({})
  //Do not use state variables in function. Stale value will be captured
  //When update state vars, must use callback style to avoid stale reference
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
        console.log(now(), `Unable to create websocket for ${host} error is ${error}`);
        return false;
      }
    } else {
      console.log(now(),`Drawdaemon is already active for ${host} ${pid}`);
      return true;//already setup
    }
    ws.onopen = () => {
      console.log(now(),`Drawdaemon connecting for ${job}`);
      ws.binaryType = 'arraybuffer';
      setWss((oldVal) => ({ ...oldVal, [job]: true }));
      wssRef.current[job] = ws;
      setBotActive((oldVal) => ({ ...oldVal, [job]: {} }))//Initialize botActive to empty object to avoid undefined error
      setTopActive((oldVal) => ({ ...oldVal, [job]: {} }))//initialize to empty
      jobRef.current[job] = { 'drawData': {}, 'layout':{}};//empty plots
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
          if(jobRef.current[job]['pid']==-1){//new session
            jobRef.current[job]['active'] = undefined //reset active page
          }
          jobRef.current[job]['pid'] = drawData['pid'];
          if(drawInfo[job]['jobname']!==host+':'+drawData['pid']){
            updateDrawInfo([job, host+':'+drawData['pid']])
          }
        }
        if ('path' in drawData) {
          jobRef.current[job]['path'] = drawData['path'];
        }
        if ('exename' in drawData) {
          jobRef.current[job]['exename'] = drawData['exename']
        }
        if ('final' in drawData) {
          jobRef.current[job]['pid'] = -1;//pid is no longer active
          const session=drawData['session']
          //Remove drawData with session less than this session
          jobRef.current[job]['drawData']=Object.fromEntries(
            Object.entries(jobRef.current[job]['drawData']).map(([key, item])=>
            [key, Object.fromEntries(//filter each fn
                Object.entries(item).filter(([key2, item2])=>
                  item2 && item2.session==session))
            ]).filter(([key, item])=>Object.keys(item).length>0)
          )
          finishGif()
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
            if (lpf.current<1 && lpf.current>0 && name in jobRef.current[job]['drawData'][fig]){
              lpfDrawdata(drawData, jobRef.current[job]['drawData'][fig][name], lpf.current);
            }
            jobRef.current[job]['drawData'][fig][name] = drawData;
            if(!jobRef.current[job]['layout'][fig]){
              jobRef.current[job]['layout'][fig]={};
            }
            if(!jobRef.current[job]['layout'][fig][name]){
              jobRef.current[job]['layout'][fig][name]={autosize: true, margin: {t:40, b:40, l:50, r:120, autoexpand:true}, 
                uirevision: Math.random(), title:{}, xaxis:{autorange: true, range:[-1,1]}, yaxis:{autorange: true, range:[-1,1]}}
            }
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
  //Persist chart user initiated zoom. 
  useEffect(() => {
    if (!chartLoaded || !chartRef.current) {
      return;
    }
    const chart = chartRef.current;
    const handleRelayout = (event) => {
      if (event['xaxis.range[0]'] !== undefined &&
          event['xaxis.range[1]'] !== undefined) {
        layout.current.xaxis.range = [
            event['xaxis.range[0]'],
            event['xaxis.range[1]']
        ];
      }
      if (event['yaxis.range[0]'] !== undefined &&
          event['yaxis.range[1]'] !== undefined) {
        layout.current.yaxis.range = [
            event['yaxis.range[0]'],
            event['yaxis.range[1]']
        ];
      }
    };
    chart.on('plotly_relayout', handleRelayout);
    //console.log(now(), 'registered handlerRelayout');
    return () => {
      chart.removeListener('plotly_relayout', handleRelayout);
    };
  }, [chartLoaded]);
  useEffect(() => {//run ones upon mount and state change
    Object.keys(drawInfo).map((job) => {
      if (drawInfo[job]) {//host:port->hostname
        if (!wssRef.current[job]) {        
          connect(drawInfo[job].hostname, job);
        }
      } else { // close plot
        if (wssRef.current[job]) {//close page
          console.warn(`Close Drawdaemon connection to ${job}`)
          wssRef.current[job].close();//drop connection
        }
        jobRef.current[job] = { 'drawData': {} , 'layout': {}};//remove plots
        Plotly.purge(chartRef.current, 0);//delete plot
      }
    })
  }, [drawInfo]);
  function startGif(filename){
    if(!chartRef.current) return;
    if(gifActive.current) finishGif()
    chartSize.current=chartRef.current.getBoundingClientRect();
    const {width, height}=chartSize.current
    if(width<=0 || height<=0) return
    //console.log(now(), 'Create new GIF instance')
    gifActive.current=1
    let g=new GIF({
      workers: 4, quality: 10, dither: false, repeat: 0,
      width: Math.round(width*window.devicePixelRatio),
      height: Math.round(height*window.devicePixelRatio),
      workerScript: '/gif.worker.js'
    })
    gif.current=g

    g.on('start', () => console.log(now(), `GIF started processing ${g.frames.length} frames for ${filename}.gif`));
    /*let lastStep=-1
    g.on('progress', (p) => {
      const step = Math.floor(p * 10);
      if (step > lastStep) {
          lastStep = step;
          console.log(now(), `GIF processed ${Math.min(step * 10, 100)}% of ${g.frames.length} frames for ${filename}.gif`);
      }
    });*/
    g.on('finished', blob =>{
      console.log(now(), `GIF finished processing ${g.frames.length} frames for ${filename}.gif`)
      const url=URL.createObjectURL(blob)
      const a=document.createElement('a')
      a.href=url
      a.download=filename+'.gif'
      a.click()
      setTimeout(()=>URL.revokeObjectURL(url), 1000)
    })
  }
  function finishGif(){
    if(gifActive.current){
      gif.current.render()
      gifActive.current=0
    }
  }
  function addCounter(layout) {
    const counter = {
        xref: 'paper',
        yref: 'paper',
        x: 1,
        y: 1,
        xanchor: 'right',
        yanchor: 'top',
        text: `${gif.current.frames.length}`,
        showarrow: false,
        font: { size: 15, color: '#000' },
    };

    layout.annotations = [
        ...(layout.annotations || []).filter(a => !a.name || a.name !== 'frameCounter'),
        { ...counter, name: 'frameCounter' }
    ];
  }

  function removeCounter(layout) {
    layout.annotations = (layout.annotations || [])
        .filter(a => a.name !== 'frameCounter');
  }
  useEffect(() => {//trigger plotting if page switched or data is updated
    if (chartRef.current && !pause[jobActive]) {
      try {
        let clear = 1;
        if (jobActive !== '' && jobRef.current[jobActive] && jobRef.current[jobActive]['drawData']
          && jobRef.current[jobActive]['drawData'][topActive[jobActive]]
          && botActive[jobActive][topActive[jobActive]]) {
          const drawData = (jobRef.current[jobActive]['drawData'][topActive[jobActive]][botActive[jobActive][topActive[jobActive]]]);
          layout.current = (jobRef.current[jobActive]['layout'][topActive[jobActive]][botActive[jobActive][topActive[jobActive]]]);
          if (drawData) {
            const id=(drawData.fig+'_'+drawData.name).replace(/[:\s]+/g, "_");
            if(figName.current!==id){
              figName.current=id;
            }
            if(drawData.cumPlot!==cumPlot || drawData.cumStart!==cumStart){
              drawData.cumPlot=cumPlot;
              drawData.cumStart=cumStart;
              //Plotly.purge(chartRef.current, 0);
            }
            drawData.climLpf=climLpfInput
            if(drawData.xlabel && drawData.xlabel.includes('Step') && xRange>0){
              const i=drawData.xlabel.indexOf('(')
              if(i>=0) drawData.xlabel=drawData.xlabel.slice(0, i)
              drawData.xlabel+='(normalized)'
              drawData.xRange=xRange
            }
            const traces  = makeTraces(drawData, layout.current);
            if (traces.length) {
              if(gifActive.current){
                addCounter(layout.current)
              }else{
                removeCounter(layout.current)
              }
              Plotly.react(chartRef.current, traces, layout.current, { responsive: true })
              .then(()=>{
                  setChartLoaded(1)
                  return gifActive.current?Plotly.toImage(chartRef.current, 
                    {format: 'png', 
                    width: chartSize.current.width, 
                    height: chartSize.current.height,
                    scale: window.devicePixelRatio}):null
                }).then((png)=>{
                  if(png){
                    const img=new Image()
                    img.src=png
                    img.decode().then(()=>{gif.current.addFrame(img, {delay: 100, copy:true})})
                  }
                })
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
  }, [figSeq, cumStart, cumPlot, pause, xRange]);

  useEffect(() => {//switch page
    try {
      if (wssRef.current[jobActive] && jobRef.current[jobActive] && jobRef.current[jobActive]['pid'] > -1
        && topActive[jobActive] && botActive[jobActive][topActive[jobActive]]) {
        wssRef.current[jobActive].send(`${jobRef.current[jobActive]['pid']}&DRAW_FIGFN&${topActive[jobActive]}&${botActive[jobActive][topActive[jobActive]]}`);//let server know the current jobActive page.
      }
      finishGif()
      setFigSeq(old => old + 1);//trigger plot update
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
        <li title="Stop receiving more data for plotting" onClick={() => { if(wssRef.current[jobActive]) wssRef.current[jobActive].close() }}>{wss[jobActive]?"⏹️" : "🛑"}</li>
        <li title="Cumulative ploting" className={cumPlot?"active":""} onClick={()=>{if(!cumInput) {setCumInput(0.1); setCumStart(0.1);}setCumPlot(oldVal=>!oldVal)}}>🎢</li>
        <Menu label={(<span style={{padding:'0.3em'}}>Options</span>)} child={<ul className="menu-list" >
          <FloatInput title='Set x axis range' name='Set x axis range' value={xRange} min="0"
            onChange={setXRange} onApply={setXRange}/>
          <FloatInput title='Set cumulative plotting starting index' name='Set cumu start' value={cumInput} step="0.1" min="0" 
            onChange={setCumInput} onApply={(value)=>{setCumStart(value); setCumPlot(true);}}/>
          <FloatInput title="Set frame update low pass filter" name="Set frame update LPF" value={lpfInput} step="0.001" min="0" max="1"
            onChange={setLpfInput} onApply={(value)=>{lpf.current=value}}/>
          <FloatInput title="Set clim update low pass filter" name="Set clim update LPF" value={climLpfInput} step="0.001" min="0" max="1"
            onChange={setClimLpfInput}/>
        </ul>}/>
        <Menu label={(<span style={{padding:'0.3em'}}>💾 Save</span>)} child={<ul className="menu-list" >
          {(gifActive.current)&&
           (<li onClick={()=>{finishGif()}}>💾 Stop Recording Animation</li>)||
           (<li onClick={()=>{startGif(figName.current)}}>🎞️ Record Animation as {figName.current}.gif</li>)
          }
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'png',filename: figName.current});}}>🖼️ Save as {figName.current}.png</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'jpeg',filename: figName.current});}}>🏞️ Save as {figName.current}.jpg</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'svg',filename: figName.current});}}>📐 Save as {figName.current}.svg</li>
          <li onClick={()=>{Plotly.downloadImage(chartRef.current, {format: 'webp',filename: figName.current});}}>🅿️ Save as {figName.current}.webp</li>
        </ul>}/>
      </ul>
      <div className="layout">
        <ul className="tab_hosts" style={{width:'8em', overflowY: 'auto'}}>{/*draw sub-notebook (vertical)*/}
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
