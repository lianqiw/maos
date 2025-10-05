"use strict"; //every variable need to be defined

const DrawDaemon = React.memo(({ drawInfo, jobActive }) => {
  //Notice that the state variables do CHANGE at every re-rendering. It just does not change within a rendering
  const { useState, useEffect, useRef } = React;
  const [topActive, setTopActive] = useState({});//active fig
  const [botActive, setBotActive] = useState({});//active name
  const [figSeq, setFigSeq] = useState(0);//figure sequence number (index of drawData). triggers plotly update.
  //Use useRef for data that is not directly related to rendering.
  const pause = useRef({});//pause plotting
  const jobRef = useRef({});//job data. 
  const wssRef = useRef({});//mutable ref object that persists for the entire lifecycle of the component.
  const chartRef = useRef(null);//for plotly. References a DOM object.

  function connect(hostname, pid) {
    if (!hostname || hostname.length == 0) return false;
    const host = split_hostname(hostname);//host:port->shortname
    const job = host + ':' + pid;
    if (hostname.indexOf(':') == -1) {
      hostname += port;
    }
    let ws;
    if (!wssRef.current[job]) {
      console.log(`Drawdaemon is connecting for ${job} at ${pcol}${hostname}`);
      try {
        ws = new WebSocket(pcol + hostname + "/xxx", "maos-drawdaemon-protocol");	/* + "/xxx" bit is for IE10 workaround */
        ws.binaryType = 'arraybuffer';
      } catch (error) {
        console.log(error);
      }
    } else {
      console.log(`Drawdaemon is already active for ${host} ${pid}`);
    }
    if (!ws) return false;
    ws.onopen = () => {
      wssRef.current[job] = ws;
      setBotActive((oldVal) => ({ ...oldVal, [job]: {} }))//Initialize botActive to empty object to avoid undefined error
      setTopActive((oldVal) => ({ ...oldVal, [job]: {} }))//initialize to empty
      jobRef.current[job] = { 'drawData': {} };//empty plots
      console.log(`Draw WebSocket connected for ${job}`);
      ws.send(pid + "&" + "DRAW" + ";");
    };

    ws.onclose = (event) => {//we keep plots when connection is closed.
      delete wssRef.current[job];
      console.log(`Draw WebSocket disconnected for ${job}`);
    };

    ws.onmessage = (event) => {
      if (typeof event.data === "string") {//instanceof does not work for string
        console.log("Got unexpected text data with bytes ", event.data.size);
      } else if (event.data instanceof Blob) {//Blob is read only
        console.log("Got unexpected blob data with bytes ", event.data.size);
      } else if (event.data instanceof ArrayBuffer) {//ArrayBuffer support writing.
        let drawData = procBuffer(event.data); //separate into function for testing
        if ('pid' in drawData) {
          jobRef.current[job]['pid'] = drawData['pid'];
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

            jobRef.current[job]['drawData'][fig][name] = drawData;
          } catch (err) {
            console.log({ err, drawData, jobRef, job });
          }
        }
        setFigSeq(old => old + 1);
      } else {
        console.error("Invalid data:", event.data);
      }
    };


    ws.onerror = (err) => {
      console.error(hostname, "WebSocket error", err);
    };
    return true;
  }//function connect

  useEffect(() => {//run ones upon mount and state change
    Object.keys(drawInfo).map((job) => {
      const hostname = drawInfo[job];//host:port->shortname
      if (hostname) {
        if (!wssRef.current[job]) {
          const pid = job.split(':').at(-1);
          connect(hostname, pid);
        }
      } else if (wssRef.current[job]) {//close page
        console.warn(`Close Draw Websocket connection to ${job}`)
        wssRef.current[job].close();//drop connection
        jobRef.current[job] = { 'drawData': {} };//remove plots
        Plotly.purge(chartRef.current, 0);//delete plot
      }
    })
  }, [drawInfo]);

  useEffect(() => {//trigger plotting if page switched or data is updated
    if (chartRef.current && !pause.current[jobActive]) {
      try {
        let clear = 1;
        if (jobActive !== '' && jobRef.current[jobActive] && jobRef.current[jobActive]['drawData']
          && jobRef.current[jobActive]['drawData'][topActive[jobActive]]
          && botActive[jobActive][topActive[jobActive]]) {
          const drawData = (jobRef.current[jobActive]['drawData'][topActive[jobActive]][botActive[jobActive][topActive[jobActive]]]);
          if (drawData) {
            const { traces, layout } = makeTraces(drawData);
            if (traces.length) {
              Plotly.newPlot(chartRef.current, traces, layout, { responsive: true });
              clear = 0;
            }
          }
        }
        if (clear) {
          Plotly.purge(chartRef.current, 0);
        }
      } catch (err) {
        Plotly.purge(chartRef.current, 0);
        console.log(err, { jobActive, jobRef, topActive, botActive });
      }
    }
  }, [jobActive, topActive, botActive, figSeq]);

  useEffect(() => {//let server know the current jobActive page.
    try {
      if (wssRef.current[jobActive] && jobRef.current[jobActive] && jobRef.current[jobActive]['pid'] > -1
        && topActive[jobActive] && botActive[jobActive][topActive[jobActive]]) {
        wssRef.current[jobActive].send(`${jobRef.current[jobActive]['pid']}&DRAW_FIGFN&${topActive[jobActive]}&${botActive[jobActive][topActive[jobActive]]}`);
      }
    } catch (err) {
      console.log({ err, jobActive, wssRef, jobRef, topActive, botActive });
    }
  }, [jobActive, topActive, botActive]);
  try {
    return jobRef.current[jobActive] && jobRef.current[jobActive]['drawData'] ? (<div>
      <ul className="inline tab_hosts"> {/*draw top-notebook (horizontal)*/}
        {Object.keys(jobRef.current[jobActive]['drawData']).sort().map((fig) => (
          <li key={fig} className={topActive[jobActive] === fig ? "active" : ""}
            onClick={() => { setTopActive(oldVal => ({ ...oldVal, [jobActive]: fig })); pause.current[jobActive] = false; }}
          >{fig}</li>))}
        <li title="Pause or Resume ploting" onClick={() => { pause.current[jobActive] = !pause.current[jobActive] }}>{pause.current[jobActive] ? "▶️" : "⏸️"}</li>
        <li title="Stop receiving more data for plotting" onClick={() => { wssRef.current[jobActive].close() }}>⏹️</li>
      </ul>
      <table border="0"><tbody><tr valign="top">
        <td><ul className="tab_hosts">{/*draw sub-notebook (vertical)*/}
          {jobRef.current[jobActive]['drawData'][topActive[jobActive]] &&
            Object.keys(jobRef.current[jobActive]['drawData'][topActive[jobActive]]).sort().map((name) => (
              <li key={name} className={botActive[jobActive][topActive[jobActive]] === name ? "active" : ""}
                onClick={() => { setBotActive(oldVal => ({ ...oldVal, [jobActive]: { ...oldVal[jobActive], [topActive[jobActive]]: name } })); pause.current[jobActive] = false; }}
              >{name}</li>))}
        </ul></td>
        <td>
          <div id="chart" ref={chartRef}></div>
        </td>
      </tr></tbody></table>
    </div>) : (<div id="chart" ref={chartRef}></div>)
  } catch (err) {
    console.log({ err, jobRef, topActive, botActive, jobActive });
  }
})

window.DrawDaemon = DrawDaemon;
