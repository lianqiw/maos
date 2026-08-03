"use strict"; //every variable need to be defined
window.port = ':80';
window.pcol = "ws://";
const iconName = {0:"❓", 1: "▶️", 2: "🕒", 3: "🆕", 4: "🕒", 11: "✅", 12: "❌", 13: "❌", 14: "❌", 15: "❌" };
//Notice: 
//Always define Components outside of App.
//Components defined within App() will always be rerendered when App() renders. Move it outside to avoid.
//Components must start with Capital letter
//Notice that the state variables do CHANGE at every re-rendering. It just does not change within a rendering call.
//Use useRef to avoid capturing a stale state variable in closure.
//State variable can only be used in 2 places: 1) in rendering jsx 2) in useEffect with variable as depency. In all other places, use current from useRef
//import {get_hostname, split_hostname} from "utils.js"; //ES Module imports
//import {DrawDaemon} from "drawdaemon.js";
window.now = () => `[${new Date().toLocaleString()}]`;
/*
function saveHost(hostname) {//add a host to localStorage.hosts
  if (!hostname) return;
  const hosts = [...new Set(localStorage.hosts.split(";"))];
  if (hosts.indexOf(hostname) == -1) {
    hosts.push(hostname);
    localStorage.hosts = hosts.filter((v) => v.length > 0).join(";");
    console.log(now(), "localStorage.hosts:", localStorage.hosts);
  }
}*/
function sec2str(seconds){
  seconds = Math.floor(seconds);

  const hours = Math.floor(seconds / 3600);
  const minutes = Math.floor((seconds % 3600) / 60);
  const secs = seconds % 60;
  
  if(hours>0){
    return String(hours).padStart(2, "0")+'h'+String(minutes).padStart(2, "0")
  }else{
    return String(minutes).padStart(2, "0")+':'+String(secs).padStart(2, "0")
  }
}
function step2str(step){
  const stepk=Math.floor(step/1000)
  if(stepk*1000==step){
    return stepk.toLocaleString()+'k'
  }else{
    return step.toLocaleString()
  }
}
function App() {
  const { useState, useEffect, useRef } = React;
  const [job, setJob] = useState([]);//save received jobs
  const [wss, setWss] = useState({});//save websocket connection status. true: connected, false: disconnected, undefined: removed
  const [active, setActive] = useState('');//current active host
  const [text, setText] = useState('');//text for input
  const [drawInfo, setDrawInfo] = useState([]);//set for draw
  //useRef are persistent and is not a new variable for every render.
  const wssRef = useRef({});//save wss information in mutatble reference.
  //window.wssRef=wssRef;//for debugging in js console
  const wssReconnectRef = useRef({});//to track reconnection attempts
  const jobRef = useRef(job);
  const activeRef = useRef(active);
  //const fullName = useRef({});//full hostname.
  const [hosts, setHosts]=useState(()=>{
    const hostname = get_hostname();
    const host = split_hostname(hostname);
    const savedHosts = localStorage.getItem("savedHosts");
    return savedHosts?JSON.parse(savedHosts):{
      [host]:hostname
    }
  })
  const rempx = parseFloat(getComputedStyle(document.documentElement).fontSize);
  const [widths, setWidths] = useState(()=>{
    const savedWidth = localStorage.getItem("savedWidth");
    const pageWidth=document.documentElement.clientWidth;
    return savedWidth?JSON.parse(savedWidth):{//for resizeable column
    "Time": rempx*7.5,
    "Host": rempx*5,
    "PID": rempx*4,
    "Start Dir": Math.max(pageWidth*0.1, rempx*10),
    "Out Dir": Math.max(pageWidth*0.1, rempx*10),
    "Low": rempx*4,
    "High": rempx*4,
    "Step": rempx*3,
    " ":rempx*1.2,
    "Progress":rempx*10
  }});
  //useEffect(()=>{//Effect function runs after React updates the DOM.
  function connect(host) {
    if (!host || host.length == 0) return false;
    var hostname=hosts[host]
    if(!hostname) return false;
    if (hostname.indexOf(':') == -1) {
      hostname += port;
    }
    let ws;
    if (wssRef.current[host] && wssRef.current[host].readyState === WebSocket.OPEN) {
      return true;//no need to setup event handlers again
    } else {
      try {
        ws = new WebSocket(pcol + hostname + "/xxx", "maos-monitor-protocol");	/* + "/xxx" bit is for IE10 workaround */
      } catch (err) {
        console.log(now(), `Monitor cannot create WebSocket to ${host}:`, err);
        return false;
      }
    }

    ws.onopen = () => {
      console.log(now(), `Monitor connected to ${host}`);
      ws.binaryType = 'arraybuffer';
      //saveHost(hostname);//save to localStorage only after successful connection.
      //fullName.current[host] = hostname;//for translation from host to hostname
      wssRef.current[host] = ws;//save the ws object only after successful connection.
      setWss((oldVal) => ({ ...oldVal, [host]: true }));
      wssReconnectRef.current[host] = 1;//reset reconnction counter
      setActive((oldVal) => oldVal ? oldVal : host);
      setJob((oldVal) => oldVal.filter((v) => v.Host != host));//remove jobs of the host
    };
    ws.onclose = () => {
      console.log(now(), `Monitor disconnected from ${host}`);
      delete wssRef.current[host];//websocket cannot be reused.
      setWss((oldVal) => ({ ...oldVal, [host]: oldVal[host] === undefined ? undefined : false }))//undefined: remove. we keep ws for reconnction
      if (wssReconnectRef.current[host] > 0 && wssReconnectRef.current[host] < 10) { //try to reconnect only 10 times max
        console.log(now(), `Monitor trying reconnect to ${host} after ${wssReconnectRef.current[host]} seconds`);
        wssReconnectRef.current[host] *= 2;//exponential backoff
        setTimeout(() => {
          connect(hostname);
        }, wssReconnectRef.current[host] * 1000);//in ms 
      }
    };
    ws.onerror = (err) => {
      console.log(now(), `Monitor connection error from ${host}`);
    };
    ws.onmessage = (event) => {
      if (typeof event.data === "string") {//instanceof does not work for string
        var j;
        if (event.data.indexOf('$') != -1) {
          j = event.data.split('$'); //newer format to avoid conflict with actual data.
        } else {
          j = event.data.split(';');
        }
        for (let ij = 0; ij < j.length; ij++) {
          if(j[ij].length==0) continue;
          const i = j[ij].split('&');
          var newdata;
          if (i.length>=2 && i[1]=='DRAW'){
            const job = `${host}:${i[0]}`;
            setDrawInfo((oldInfo) => ({ ...oldInfo, [job]: {hostname:hosts[host], jobname:job} }))
            setActive(job);
            //console.log(now(), "Answer drawing for "+job)
            continue;//skip setJob
          }else if (i.length == 3 && i[1]==='PATH') {//0: pid, 1: PATH, 2: path+args
            let startdir = "";
            const ip = i[2].indexOf(" ");
            if (ip != -1) {
              startdir = i[2].substring(0, ip)+" ";
              i[2] = i[2].substring(ip);
            }
            const io = i[2].lastIndexOf("-o");
            let dirout = "";
            if (io != -1) {
              dirout = " -o "+i[2].substring(io + 2).trim().split(' ')[0]+" "
              i[2] = i[2].replace(/[ \t]+-o[ \t]+[^ \t]+/g, " ").trim();
            }
            //newdata = { PID: i[0], Host: host, "Start Dir": startdir, "Arguments": i[2], "Out Dir": dirout, status: 0 };
            newdata = { PID: i[0], Host: host, "Arguments":startdir+dirout+i[2], status: 0 };
          } else if (i.length == 14 && i[1]==='STATUS') {//0:pid, 1:STATUS, 2:pidnew, 3:status, 4:start time,5: errhi, 6:errlo, 7:iseed, 8:nseed, 9:isim, 10:nsim, 11:rest, 12:tot, 13:step timing
            i[3] = parseInt(i[3]);//status. 1: running, 2: wait, 3: started, 4: queued. 11: finished. 12: crashed. 13: to kill; 14: remove; 15: killed;
            if (i[3] != 14) {//14: remove
              i[2] = parseInt(i[2]);//pid
              i[11] = parseInt(i[11]);//remaining time
              i[12] = parseInt(i[12]);//total time
              const icon = (i[12] > 0 || i[3] != 11) ? iconName[i[3]] : "⏩";
              const frac = Math.round(100 * (1 - (i[12] == 0 ? 1 : i[11] / i[12]))) + '%';
              let prog = i[12] == 0 ? "" : i[7] + '/' + i[8] + ' ' + step2str(i[9]) + '/' + step2str(i[10])+' '+sec2str(i[3]==1?i[11]:i[12])
              if (i[5] === '0.00') i[5] = '';
              if (i[6] === '0.00') i[6] = '';
              if (i[13] === '0.000') i[13] = '';
              newdata = { PID: i[2], Host: host, status: i[3]?i[3]:12, Time: i[4], High: i[5], Low: i[6], Step: i[13], prog: prog, frac: frac, icon: icon };
            }
          } else {
            continue;//invalid data
          }
          try {
            setJob((oldValue) => {
              const ind = oldValue.findIndex((v) => v.PID == i[0] && v.Host === host);
              if (newdata === undefined){//remove
                if(ind != -1) {//found
                  return oldValue.filter((v, i) => i != ind);
                }else{
                  return oldValue;//not found
                }
              } else if (ind == -1) {//append
                return [...oldValue, newdata];
              } else {//update
                return oldValue.map((v, i) => (i === ind) ? { ...v, ...newdata } : v);
              }
            });
          } catch (err) {
            console.log(now(), { err, job: jobRef.current, newdata })
          }
        }
      } else if (event.data instanceof Blob) {//Blob is read only
        console.log(now(), "Got unexpected blob data with bytes ", event.data.size);
      } else if (event.data instanceof ArrayBuffer) {//ArrayBuffer support writing.
        console.log(now(), "Got unexpected arraybuffer data with bytes ", event.data.size);
      } else {
        console.error("Invalid data:", event.data);
      }
    };
    return true;
  }//function connect
  function removeHost(host) {//remove a host from localStorage.hosts
    if (!host) return;
    //const hostname = fullName.current[host];
    if (wssRef.current[host]) wssRef.current[host].close();
    wssReconnectRef.current[host] = 0;//prevent auto reconnection
    setWss((oldVal) => ({ ...oldVal, [host]: undefined }));//hide from list
    //const hosts = [...new Set(localStorage.hosts.split(";"))];//unique entries
    //localStorage.hosts = hosts.filter((v) => v != hostname && v.length > 0).join(";");
    //console.log(now(), "localStorage.hosts:", localStorage.hosts);
    setHosts((oldVal)=>({...oldVal, [host]: undefined}))
  }
  useEffect(() => {//run once upon mount
    /*const hostname = get_hostname();//connect current host
    if (!localStorage.hosts) localStorage.hosts = "";
    let hosts = [...new Set(localStorage.hosts.split(";"))]
    if (hosts.indexOf(hostname) == -1) {
      hosts.push(hostname);
    }
    hosts.filter((v) => v.length > 0).forEach((v) => connect(v));//connect to other hosts
    */
    (async () => {
      const { ZstdSimple, ZstdStream } = await zstdCodec.ZstdInit();
      window.ZstdSimple=ZstdSimple
      window.ZstdStream=ZstdStream
    // Now use ZstdSimple.compress(...) or ZstdSimple.decompress(...)
    })();
  }, []);

  useEffect(() => {
    jobRef.current = job;//update reference
  }, [job])

  useEffect(() => {
    activeRef.current = active;//update reference
  }, [active])
  useEffect(()=>{
    localStorage.setItem("savedWidth", JSON.stringify(widths))
  }, [widths])
  useEffect(()=>{
    Object.keys(hosts).map(host=>{
      if(hosts[host] && !wssRef.current[host]){
        connect(host)
      }
    })
    localStorage.setItem("savedHosts", JSON.stringify(hosts))
  }, [hosts])
  function cmdHostPid(host, PID, cmd) {
    if (wssRef.current[host] && wssRef.current[host].readyState === WebSocket.OPEN) {
      if (cmd === "DRAW") {
        const job = `${host}:${PID}`;
        setDrawInfo((oldInfo) => ({ ...oldInfo, [job]: {hostname: hosts[host], jobname:job} }))
        setActive(job);
        //console.log(now(), "Initiate drawing for "+job)
      } else if (cmd === "KILL_ASK") {
        if (window.confirm("Kill the job " + PID + "?")) cmdHostPid(host, PID, "KILL");
      } else {
        wssRef.current[host].send(PID + "&" + cmd + ";");
        console.log(now(), "Sending " + cmd + " to " + host + " for " + PID);
      }
    }
  }
  function cmdHost(host, cmd) {
    //Sending command to the host will change jobRef. So we create the entire
    //command string together for all the jobs
    const hosts = host === "" ? Object.keys(wssRef.current) : [host];
    const cmdstr = hosts.reduce((bcc, h) => ({
      ...bcc, [h]:
        jobRef.current.reduce((acc, v) => {
          if (v.Host === h) {
            if (acc.length>4000) {
              return acc;//prevent overflow
            } else if (cmd === "clear_all" && v.status >= 11) {
              return acc + `${v.PID}&REMOVE;`;
            } else if (cmd === "clear_finished" && v.status == 11) {
              return acc + `${v.PID}&REMOVE;`;
            } else if (cmd === "clear_skipped" && v.status == 11 && v.prog === "") {
              return acc + `${v.PID}&REMOVE;`;
            } else if (cmd === "clear_crashed" && v.status > 11) {
              return acc + `${v.PID}&REMOVE;`;
            } else if (cmd === "kill_all" && v.status < 11) {
              return acc + `${v.PID}&KILL;`;
            }
          }
          return acc;
        }, "")
    }), {});
    Object.keys(cmdstr).forEach((host) => {
      const command = cmdstr[host];
      if (command && command.length > 0 && wssRef.current[host] && wssRef.current[host].readyState === WebSocket.OPEN) {
        if (cmd !== "kill_all" || window.confirm("Kill all jobs on " + (host === "" ? "all" : host) + "?")) {
          wssRef.current[host].send(command);
          console.log(now(), "Sending " + command + " to " + host);
        }
      }
    });
  }

  const columns = ["Time", "Host", "PID", "Arguments", "Low", "High", "Step", " ", "Progress"];
  try{
    return (
      <div>
        <ul className="inline tab_hosts">
          <Menu label={(<span><img src="icon-monitor.png" alt="icon"></img><span>Menu</span></span>)} child={
            <ul className="menu-list" >
              <li onClick={() => { cmdHost(active, "clear_all"); }}>☑️ Clear All Jobs on {active === "" || wss[active] === undefined ? "all" : active}</li>
              <li onClick={() => { cmdHost(active, "clear_finished"); }}>✅ Clear Finished Jobs on {active === "" || wss[active] === undefined ? "all" : active}</li>
              <li onClick={() => { cmdHost(active, "clear_skipped"); }}>⏩ Clear Skipped Jobs on {active === "" || wss[active] === undefined ? "all" : active}</li>
              <li onClick={() => { cmdHost(active, "clear_crashed"); }}>❌ Clear Crashed Jobs on {active === "" || wss[active] === undefined ? "all" : active}</li>
              <li onClick={() => { cmdHost(active, "kill_all"); }}>🟥 Kill All Jobs on {active === "" || wss[active] === undefined ? "all" : active}</li>
              {job.filter((row) => (row.status == 1 || row.status == 3)).map((row) =>
                (<li key={row.Host + row.PID} onClick={() => { cmdHostPid(row.Host, row.PID, "DRAW") }}>▶️ Plot {row.PID} at {row.Host}</li>))}
            </ul>
          }></Menu>
          <li className={active === "" ? "active" : ""} onClick={() => setActive("")}>
            <span>▶️ Active</span>
          </li>
          {Object.keys(hosts).filter((v) => hosts[v] != undefined).map((host) => (//List of hosts
            <li key={host} className={active === host ? "active" : ""}>
              <span title="Connect to host" onClick={() => { setActive(host); if (!wss[host]) connect(host); }}>{wss[host] ? "🟢" : "🔴"}</span>
              <span title="Switch to host" onClick={() => { setActive(host); }}>{host}</span>
              <span title="Remove Host" onClick={() => { setActive(""); removeHost(host); }}>&nbsp;⛌</span>
            </li>
          ))}
          <li><form onSubmit={(e) => { e.preventDefault(); setText(''); if (text.length) { setHosts((oldVal)=>({...oldVal, [split_hostname(text)]:text}))}}}>
            <input style={{ width: '10em' }} id="add_host" value={text} onChange={e => setText(e.target.value)}></input>
            <button type="submit">Connect</button></form>
          </li>
          {Object.keys(drawInfo).filter((job) => (drawInfo[job])).map((job) => (
            <li key={job} className={active === job ? "active" : ""}>
              <span onClick={() => { setActive(job) }}>{drawInfo[job].jobname}</span>
              <span onClick={() => { setDrawInfo((oldInfo) => ({ ...oldInfo, [job]: undefined })); setActive(""); }}>&nbsp;⛌</span>
            </li>))}
        </ul>
        {!active.includes(':') && (
          <table className="monitor">
            <colgroup>
            {columns.map((cl) => <col key={cl} title={cl} style={{width:widths[cl]}}></col>)}
            </colgroup>
            <thead>
              <tr>{columns.map((col) => col in widths?<ResizableColumn
              width={widths[col]} title={col} key={col}
              onWidthChange={(w) =>
                setWidths({ ...widths, [col]: w })
              }/>:<th key={col}>{col}</th>)}
              </tr>
            </thead>
            <tbody>
              {job.filter((row) => ((active.length == 0 && row.status < 10) || row.Host === active)).map((row, i) => (
                <tr key={row.PID}>
                  {columns.slice(0,-2).map((col) => <td key={col} title={row[col]}>{row[col]}</td>)}
                  <td className="jobIcon" key="icon" onClick={() => { cmdHostPid(row.Host, row.PID, row.status > 10 ? "REMOVE" : "KILL_ASK") }}>{row.icon}</td>
                  <Progress key="progress" text={row.prog} frac={row.frac}></Progress>
                </tr>
              ))}
            </tbody>
          </table>
        )}
        <DrawDaemon drawInfo={drawInfo} jobActive={active} 
        updateDrawInfo={(info)=>{setActive(info[0]); setDrawInfo((oldInfo) => ({ ...oldInfo, [info[0]]: {...oldInfo[info[0]], jobname: info[1]}}))}}></DrawDaemon>
      </div>
    );
  }catch(err){
    console.log(Now(), "Failed to create jsx:", {wss, job, drawInfo} )
  }
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
