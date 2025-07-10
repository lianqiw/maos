function get_hostname() {
	/*
	* We open the websocket encrypted if this page came on an
	* https:// url itself, otherwise unencrypted
	*/
	console.log(window.location.href);
	var u = window.location.href;
	var pcol = ""
	if (u.substring(0, 5) == "https") {
		pcol = "wss://";
		u = u.substring(8);
	} else {
		pcol = "ws://";
		if (u.substring(0, 4) == "http"){
			u = u.substring(7);
		}
	}

	u = u.split('/')[0];

	/* + "/xxx" bit is for IE10 workaround */
	return pcol+u+"/xxx";
}
function App() {
  const [data, setData] = React.useState([]);
  const [connected, setConnected] = React.useState(false);

  React.useEffect(() => {
    const ws = new WebSocket(get_hostname(), "maos-monitor-protocol");

    ws.onopen = () => {
      setConnected(true);
      console.log("WebSocket connected");
    };

    ws.onmessage = (event) => {
	  var j;
      try {
		if (event.data.indexOf('$')!=-1){
			j = event.data.split('$'); //newer format to avoid conflict with actual data.
		}else{
			j = event.data.split(';');
		}
		if(j.length>0){
			const i=j[0].split('&');
        	setData(i);
		}
      } catch (err) {
        console.error("Invalid data:", event.data, err);
      }
    };

    ws.onclose = () => {
      setConnected(false);
      console.warn("WebSocket disconnected");
    };

    ws.onerror = (err) => {
      console.error("WebSocket error", err);
    };

    return () => ws.close();
  }, []);

  const columns = data.length ? Object.keys(data[0]) : [];

  return (
    <div style={{ padding: 20 }}>
      <h2>Monitor (WebSocket)</h2>
      <p>Status: {connected ? "🟢 Connected" : "🔴 Disconnected"}</p>

      {columns.length > 0 ? (
        <table border="1" cellPadding="4">
          <thead>
            <tr>{columns.map((col) => <th key={col}>{col}</th>)}</tr>
          </thead>
          <tbody>
            {data.map((row, i) => (
              <tr key={i}>
                {columns.map((col) => <td key={col}>{row[col]}</td>)}
              </tr>
            ))}
          </tbody>
        </table>
      ) : (
        <p>No data yet</p>
      )}
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
