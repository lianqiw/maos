"use strict"; //every variable need to be defined
window.port = ':80';
window.pcol = "ws://";
const iconName = { 1: "▶️", 2: "🕒", 3: "🆕", 4: "🕒", 11: "✅", 12: "❌", 13: "❌", 14: "❌", 15: "❌" };
//Notice: 
//Always define Components outside of App.
//Components defined within App() will always be rerendered when App() renders. Move it outside to avoid.
//Components must start with Capital letter
//Notice that the state variables do CHANGE at every re-rendering. It just does not change within a rendering call.
//Use useRef to avoid capturing a stale state variable in closure.
//State variable can only be used in 2 places: 1) in rendering jsx 2) in useEffect with variable as depency. In all other places, use current from useRef
//import {get_hostname, split_hostname} from "utils.js"; //ES Module imports
//import {DrawDaemon} from "drawdaemon.js";
function removeHost(hostname) {
  if (!hostname) return;
  const hosts = [...new Set(localStorage.hosts.split(";"))];//unique entries
  localStorage.hosts = hosts.filter((v) => v != hostname && v.length > 0).join(";");
}
function App() {
  const { useState, useEffect, useRef } = React;
  const [job, setJob] = useState([]);//save received jobs
  const [wss, setWss] = useState({});//save websocket information
  const [active, setActive] = useState('');//current active host
  const [text, setText] = useState('');//text for input
  const [drawInfo, setDrawInfo] = useState([]);//set for draw
  //useRef are persistent and is not a new variable for every render.
  const wssRef = useRef({});//for immediate testing, do not wait for state update
  const jobRef = useRef(job);
  const activeRef = useRef(active);
  const fullName = useRef({});//full hostname.
  //useEffect(()=>{//Effect function runs after React updates the DOM.
  function connect(hostname) {
    if (hostname.length == 0) return false;
    const host = split_hostname(hostname);//host:port->shortname
    if (hostname.indexOf(':') == -1) {
      hostname += port;
    }
    let ws;
    if (!wssRef.current[host]) {
      try {
        ws = new WebSocket(pcol + hostname + "/xxx", "maos-monitor-protocol");	/* + "/xxx" bit is for IE10 workaround */
        ws.binaryType = 'arraybuffer';
      } catch (error) {
        console.log(error);
        return false;
      }
    }
    if (!ws) return false;

    ws.onopen = () => {
      fullName.current[host] = hostname;
      wssRef.current[host] = ws;
      setWss((oldVal) => ({ ...oldVal, [host]: ws }));
      setActive((oldVal) => oldVal ? oldVal : host);
      console.log(host, "Monitor WebSocket connected", { ...ws });
      if (localStorage.hosts.indexOf(hostname) == -1) {
        localStorage.hosts += ";" + hostname;
      }
      console.log({ ...wssRef.current });
    };
    ws.onclose = () => {
      console.log({ ...wssRef.current });
      delete wssRef.current[host];
      setWss((oldVal) => (oldVal[host] ? { ...oldVal, [host]: false } : { ...oldVal, [host]: undefined }))//false: disconnected. undefined: removed.
      setJob((oldVal) => oldVal.filter((v) => v.Host != host));//remove jobs of the host
      console.warn(host, "Monitor WebSocket disconnected", ws);
    };
    ws.onerror = (err) => {
      removeHost(hostname);
      console.error(hostname, "WebSocket error", err);
    };
    ws.onmessage = (event) => {
      if (typeof event.data === "string") {//instanceof does not work for string
        var j;
        if (event.data.indexOf('$') != -1) {
          j = event.data.split('$'); //newer format to avoid conflict with actual data.
        } else {
          j = event.data.split(';');
        }
        if (j.length > 0) {
          for (let ij = 0; ij < j.length - 1; ij++) {
            const i = j[ij].split('&');
            var newdata;
            if (i.length == 3) {
              let path = "";
              const ip = i[2].indexOf(" ");
              if (ip != -1) {
                path = i[2].substring(0, ip);
                i[2] = i[2].substring(ip);
              }
              const io = i[2].lastIndexOf("-O");
              let dirout = "";
              if (io != -1) {
                dirout = i[2].substring(io + 2).trim().split(' ')[0];
                i[2] = i[2].replace(/[ \t]+-o[ \t]+[^ \t]+/g, " ").trim();
              }
              newdata = { PID: i[0], Host: host, Path: path, "Dir Out": dirout, Args: i[2], status: 0 };
            } else if (i.length == 14) {
              i[3] = parseInt(i[3]);//status. 1: running, 2: wait, 3: started, 4: queued. 11: finished. 12: crashed. 13: to kill; 14: remove; 15: killed;
              if (i[3] != 14) {//14: remove
                i[2] = parseInt(i[2]);//pid
                i[11] = parseInt(i[11]);//remaining time
                i[12] = parseInt(i[12]);//total time
                const icon = (i[12] > 0 || i[3] != 11) ? iconName[i[3]] : "⏩";
                const prog = i[12] == 0 ? "" : i[7] + '/' + i[8] + ' ' + i[9] + '/' + i[10] + ' ' + i[11] + '/' + i[12];
                const frac = Math.round(100 * (1 - (i[12] == 0 ? 1 : i[11] / i[12]))) + '%'
                newdata = { PID: i[2], Host: host, status: i[3], Time: i[4], High: i[5], Low: i[6], prog: prog, frac: frac, icon: icon };
              }
            } else {
              continue;//invalid data
            }
            setJob((oldValue) => {
              const ind = oldValue.findIndex((v) => v.PID == i[0] && v.Host === host);
              if (newdata === undefined && ind != -1) {//remove
                return oldValue.filter((v, i) => i != ind);
              } else if (ind == -1) {//append
                return [...oldValue, newdata];
              } else {//update
                return oldValue.map((v, i) => (i === ind) ? { ...v, ...newdata } : v);
              }
            });
          }
        }
      } else if (event.data instanceof Blob) {//Blob is read only
        console.log("Got unexpected blob data with bytes ", event.data.size);
      } else if (event.data instanceof ArrayBuffer) {//ArrayBuffer support writing.
        console.log("Got unexpected arraybuffer data with bytes ", event.data.size);
      } else {
        console.error("Invalid data:", event.data);
      }
    };


    return true;
  }//function connect
  useEffect(() => {//run once upon mount
    const hostname = get_hostname();//connect current host
    let hosts_old;
    if(localStorage.hosts){
      hosts_old = localStorage.hosts.split(";")
    }else{
      hosts_old=[];
    }
    hosts_old.push(hostname);
    const hosts = [...new Set(hosts_old)];//unique entries
    console.log(hosts);
    localStorage.hosts = hosts.filter((hostname) => connect(hostname)).join(";");//keep only valid hosts
    console.log({ storage: localStorage.hosts })
  }, []);

  useEffect(() => {
    jobRef.current = job;//update reference
  }, [job])

  useEffect(() => {
    activeRef.current = active;//update reference
  }, [active])

  function cmdHostPid(host, PID, cmd) {
    if (wssRef.current[host]) {
      if (cmd === "DRAW") {
        const job = `${host}:${PID}`;
        setDrawInfo((oldInfo) => ({ ...oldInfo, [job]: fullName.current[host] }))
        setActive(job);
      } else if (cmd === "KILL_ASK") {
        if (window.confirm("Kill the job " + PID + "?")) cmdHostPid(host, PID, "KILL");
      } else {
        wssRef.current[host].send(PID + "&" + cmd + ";");
        console.log("Sending " + cmd + " to " + host + " for " + PID);
      }
    }
  }
  function cmdHost(host, cmd) {
    let confirm = 0; //0: not confirmed. 1: confirmed.
    let terminate = 0;//1: skip the mapping
    console.log(cmd, host)
    jobRef.current.map((v) => {
      if ((host === "" || v.Host === host) && !terminate) {
        if (cmd === "clear_finished" && v.status == 11) {
          cmdHostPid(v.Host, v.PID, "REMOVE");
        } else if (cmd === "clear_skipped" && v.status == 11 && v.tot == 0) {
          cmdHostPid(v.Host, v.PID, "REMOVE");
        } else if (cmd === "clear_crashed" && v.status > 11) {
          cmdHostPid(v.Host, v.PID, "REMOVE");
        } else if (cmd === "kill_all" && v.status < 11) {
          if (!confirm) {
            confirm = window.confirm("Kill all jobs on " + (host===""?"all":host) + "?");
            if (!confirm) {//only ask once
              terminate=1;
            }
          }
          if(confirm){
            cmdHostPid(v.Host, v.PID, "KILL");
          }
        } else if (cmd === "plot_job" && (v.status == 1 || v.status == 3)) {
          cmdHostPid(v.Host, v.PID, "DRAW");
        }
      }
    })
  }

  const columns = ["Time", "Host", "PID", "Path", "Dir Out", "Low", "High", "Step"];
  const cn = { Time: "", Host: "", PID: "", Path: "tdpath", "Dir Out": "tdout", Low: "", High: "", Step: "" };
  return (
    <div>
      <ul className="inline tab_hosts">
        <Menu label={(<div><img src="icon-monitor.png" alt="icon"></img><span>Menu</span></div>)} child={
          (active===""||wss[active])&&(//Menubar
          <ul className="menu-list" >
            <li onClick={() => { cmdHost(active, "clear_finished"); }}>✅ Clear Finished Jobs on {active === "" ? "all" : active}</li>
            <li onClick={() => { cmdHost(active, "clear_skipped"); }}>⏩ Clear Skipped Jobs on {active === "" ? "all" : active}</li>
            <li onClick={() => { cmdHost(active, "clear_crashed"); }}>❌ Clear Crashed Jobs on {active === "" ? "all" : active}</li>
            <li onClick={() => { cmdHost(active, "kill_all"); }}>🟥 Kill All Jobs on {active === "" ? "all" : active}</li>
            {job.filter((row) => (row.status == 1 || row.status == 3)).map((row) =>
              (<li key={row.Host + row.PID} onClick={() => { cmdHostPid(row.Host, row.PID, "DRAW") }}>▶️ Plot {row.PID} at {row.Host}</li>))}
          </ul>
        )}></Menu>
        <li className={active === "" ? "active" : ""} onClick={() => setActive("")}>
          <span>▶️ Active</span>
        </li>
        {Object.keys(wss).filter((v) => wss[v] !== undefined).map((host) => (//List of hosts
          <li key={host} className={active === host ? "active" : ""}>
            <span title="Switch to host" onClick={() => { setActive(host); if (!wss[host]) connect(fullName.current[host]); }}>{wss[host] ? "🟢" : "🔴"}{host}</span>
            <span title="Remove Host" onClick={() => {
              if (wss[host]) wss[host].close(); setWss((oldVal) => ({ ...oldVal, [host]: undefined })); console.log(wss);
              removeHost(fullName.current[host]); setActive("");
            }}> ⛌</span>
          </li>
        ))}
        <li><input id="add_host" value={text} onChange={e => setText(e.target.value)}></input>
          <button onClick={() => { setText(''); if (text.length) { connect(text); } }}>Connect</button>
        </li>
        {Object.keys(drawInfo).filter((info) => (drawInfo[info])).map((info) => (
          <li key={info} className={active === info ? "active" : ""}>
            <span onClick={() => { setActive(info) }}>{info}</span>
            <span onClick={() => { setDrawInfo((oldInfo) => ({ ...oldInfo, [info]: undefined })); setActive(""); }}>⛌</span>
          </li>))}
      </ul>
      {!active.includes(':') && job.filter((row) => ((active.length == 0 && row.status < 10) || row.Host === active)).length > 0 ? (
        <table className="monitor">
          <thead>
            <tr>{columns.map((col) => <th key={col} className={cn[col]}>{col}</th>)}
              <th key="icon"></th>
              <th key="progress">Progress</th>
            </tr>
          </thead>
          <tbody>
            {job.filter((row) => ((active.length == 0 && row.status < 10) || row.Host === active)).map((row, i) => (
              <tr key={row.PID}>
                {columns.map((col) => <td key={col} className={cn[col]}>{row[col]}</td>)}
                <td onClick={() => { cmdHostPid(row.Host, row.PID, row.status > 10 ? "REMOVE" : "KILL_ASK") }}>{row.icon}</td>
                <Progress text={row.prog} frac={row.frac}></Progress>
              </tr>
            ))}
          </tbody>
        </table>
      ) : null
      }
      <DrawDaemon drawInfo={drawInfo} jobActive={active}></DrawDaemon>
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
