"use strict"; //every variable need to be defined
const draw = {
  start: 0,
  data: 1,
  heartbeat: 2,
  points: 3,
  style: 4,
  circle: 5,
  limit: 6,
  fig: 7,
  name: 8,
  title: 9,
  xlabel: 10,
  ylabel: 11,
  zlim: 12,
  legend: 13,
  xylog: 14,
  figfn: 15,
  pause: 16,
  resume: 17,
  final: 18,
  float: 19,
  frame: 20,
  single: 21,
  udpport: 22,
  init: 23,
  pid: 24,
  zlog: 25,
  path: 26,
  exename: 27,
  end: 100,
  entry: 9999
};
var byteFloat = 4;

function linspace(start, stop, length) {
  if (stop === undefined) {
    stop = start;
    start = 0;
  }
  if (length === undefined) length = Math.round(stop - start);
  const step = (length == 1) ? 1 : (stop - start) / (length - 1);
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
function mkRGB(val) {
  const r = (val >> 16) & 0xFF;
  const g = (val >> 8) & 0xFF;
  const b = val & 0xFF;
  return `rgb(${r},${g},${b})`;
}
function mkCircle(cir) {
  const cx = cir[0];
  const cy = cir[1];
  const radius = cir[2];
  const color = cir[3];
  const np = 128;
  const x = new Float32Array(np + 1);
  const y = new Float32Array(np + 1);
  const dt = Math.PI * 2. / np;
  for (let i = 0; i < np; i++) {
    const angle = dt * i;
    x[i] = Math.cos(angle) * radius + cx;
    y[i] = Math.sin(angle) * radius + cy;
  }
  x[np] = x[0];
  y[np] = y[0];
  return { x, y, mode: 'lines', color: mkRGB(color), showlegend: false };
}
function mkCircles(cirs) {
  const ncir = cirs.length >> 2;
  const traces = [];
  for (let i = 0; i < ncir; i++) {
    traces.push(mkCircle(cirs.slice(i * 4, i * 4 + 4)));
  }
  return traces;
}
function reshapeTo2D(p, nx, ny) {
  const res = []
  for (let iy = 0; iy < ny; iy++) {
    res.push(p.slice(iy * nx, (iy + 1) * nx));
  }
  //console.log(`reshapeTo2D nx=${nx} ny=${ny} res[0].length=${res[0].length} res.length=${res.length}`)
  return res;
}
function procBuffer(wsBuf) {
  //console.log("Got ArrayBuffer data with bytes ", wsBuf.byteLength);
  const dataView = new DataView(wsBuf);
  var buf = { pos: 0, cmd: 0 };
  var drawData = {};
  function getKey() {
    Object.keys(draw).find(key => draw[key] === buf.cmd)
  }
  function getInt() {
    let v = dataView.getInt32(buf.pos, true); buf.pos += 4;
    return v;
  }
  function getByteArray(e = -1) {
    const nx = e > -1 ? e : getInt();
    let v;
    try {
      v = new Uint8Array(wsBuf, buf.pos, nx); buf.pos += nx;
    } catch (err) {
      console.log(err, getKey(), wsBuf.byteLength, buf.pos, nx, e);
    }
    return v;
  }
  function getIntArray(e = -1) {
    const nx = e > -1 ? e : getInt();
    const newBuffer = wsBuf.slice(buf.pos, buf.pos + nx * 4); buf.pos += nx * 4;
    const v = new Int32Array(newBuffer);
    return v;
  }
  function getFloatArray(e = -1) {
    const nx = e > -1 ? e : getInt();
    const newBuffer = wsBuf.slice(buf.pos, buf.pos + nx * byteFloat); buf.pos += nx * byteFloat;
    let v;
    try {
      if (byteFloat == 4) {//use slice to avoid byteoffset alignment error.
        v = new Float32Array(newBuffer);
      } else if (byteFloat == 8) {
        v = new Float64Array(newBuffer);
      }
    } catch (error) {
      console.log(getKey(buf.cmd), "newFloatArray fails with", error,
        newBuffer, newBuffer.byteLength, buf.pos, nx, byteFloat, wsBuf.byteLength);
    }
    return v;
  }
  function getStr() {
    return String.fromCharCode.apply(null, getByteArray().slice(0, -1));//drop terminating \0
  }
  let len;
  while (buf.pos < wsBuf.byteLength) {
    buf.cmd = getInt()
    switch (buf.cmd) {
      case draw.entry:
        len = getInt();//length
        continue;
        break;
      case draw.frame:
        buf.pos += 16;//skip 4 ints
        break;
      case draw.start:
        break;
      case draw.data: {
        let nx = getInt();
        let ny = getInt();
        drawData['data'] = reshapeTo2D(getFloatArray(nx * ny), nx, ny);
      } break;
      case draw.heartbeat:
        break;
      case draw.points: {
        let nx = getInt();
        let ny = getInt();
        let square = getInt();
        if (!drawData['pts']) {
          drawData['pts'] = [];
        }
        drawData['pts'].push({ p: getFloatArray(nx * ny), nx: nx, ny: ny, square: square });
      } break;
      case draw.style:
        drawData['style'] = getIntArray();
        break;
      case draw.circle: {
        const nc = getInt()
        drawData['circle'] = getFloatArray(4 * nc);
      } break;
      case draw.limit:
        drawData['limit'] = getFloatArray(4);
        drawData['limitManual'] = 1;
        break;
      case draw.fig:
        drawData['fig'] = getStr();
        break;
      case draw.name:
        drawData['name'] = getStr();
        break;
      case draw.title:
        drawData['title'] = getStr();
        break;
      case draw.xlabel:
        drawData['xlabel'] = getStr();
        break;
      case draw.ylabel:
        drawData['ylabel'] = getStr();
        break;
      case draw.zlim:
        drawData['zlim'] = getFloatArray(2);
        drawData['zlimManual'] = 1;
        drawData['zlogLast'] = 0;
        break;
      case draw.legend:
        drawData['legend'] = []
        for (let il = 0; il < drawData['pts'].length; il++) {
          drawData['legend'].push(getStr());
        }
        break;
      case draw.xylog:
        drawData['xylog'] = getByteArray(2);
        break;
      case draw.final:
        //console.log("DRAW_FINAL");
        drawData['final']=1;
        drawData['id']=drawData['fig']+drawData['name']
        break;
      case draw.float:
        byteFloat = getInt();
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
        drawData['pid'] = getInt();
        break;
      case draw.zlog:
        drawData['zlog'] = getInt();
        break;
      case draw.path:
        drawData['path'] = getStr();
        break;
      case draw.exename:
        drawData['exename'] = getStr();
        break;
      case draw.end:
        //console.log(`drawData (${wsBuf.byteLength} bytes):`, drawData)
        //start making plot
        return drawData;
      default:
        console.warn("unknown command", buf.cmd);
        buf.pos += len;//skip its payload
    }//switch
  }//while
  return drawData;
}//function procbuffer
function makeCumu(y, cumStart, cumPlot){
  if(cumPlot){
    if(cumStart<0){
      cumStart=0; 
    }else if (cumStart<1) {
      cumStart=Math.round(cumStart*y.length);
    }else {
      cumStart=Math.round(cumStart);
    }
    let y2;
    if (byteFloat == 4) {//use slice to avoid byteoffset alignment error.
      y2 = new Float32Array(y.length-cumStart);
    } else if (byteFloat == 8) {
      y2 = new Float64Array(y.length-cumStart);
    }
    let cum=0;
    for(let i=0; i<y2.length; i++){
      cum=cum+y[i+cumStart]**2;
      y2[i]=Math.sqrt(cum/(i+1));
    }
    const x2=Array.from({length:y2.length}, (_, index)=>cumStart+index);
    return {x:x2,y:y2};
  }else{
    return {y:y};
  }
}
function makeTraces(drawData) {
  const layout = {
    title: {'text':('title' in drawData?drawData['title']:"")}, 
    xaxis: { showgrid: true, showline: true, mirror: true, ticks: "inside", zeroline: false, 
      type: ('xylog' in drawData && drawData.xylog[0]==='y'.charCodeAt(0))?'log':'linear',
      title: {'text':('xlabel' in drawData?drawData['xlabel']:"")}, 
    },
    yaxis: { showgrid: true, showline: true, mirror: true, ticks: "inside", zeroline: false,
      type: ('xylog' in drawData && drawData.xylog[1]==='y'.charCodeAt(0))?'log':'linear',
      title: {'text':('ylabel' in drawData?drawData['ylabel']:"")}, 
    },
    autosize: true,
    margin: {t:40, b:40, l:50, r:20},
    uirevision: drawData['id'],
  }
  let traces = [];
  let square = 0;
  let ratio = 1;
  if (drawData['data']) {
    const nx = drawData['data'][0].length;
    const ny = drawData['data'].length;
    if(nx && ny){
      square = 1;
      let trace={z: drawData['data'], type: 'heatmap', colorscale: 'Jet', showlegend: false };
      if (drawData['limit']) {
        const x = linspace(drawData['limit'][0], drawData['limit'][1], nx)
        const y = linspace(drawData['limit'][2], drawData['limit'][3], ny)
        ratio=(drawData['limit'][3]-drawData['limit'][2])/(drawData['limit'][1]-drawData['limit'][0]);
        trace={x,y,...trace}
      }else{
        ratio=ny/nx
      }
      traces.push(trace);
    }
  } else if (drawData['pts']) {
    const showlegend = drawData['pts'].length > 1 ? 1 : 0;
    traces = traces.concat(drawData['pts'].map((pts, ipts) => {
      const extra={name: drawData['legend'] ? drawData['legend'][ipts] : "", showlegend};
      if (pts.ny == 2) {
        if (pts.square) square = 1;
        return {
          x: pts.p.slice(0, pts.nx), y: pts.p.slice(pts.nx, 2 * pts.nx), 
          mode: square?"markers":"lines", marker: { symbol: 'x' }, ...extra
        }
      } else{
        return { ...makeCumu(pts.p, drawData.cumStart, drawData.cumPlot), mode: "lines",  ...extra}
      }
    }));
  }
  if (drawData['circle']) {
    traces = traces.concat(mkCircles(drawData['circle']));
  }
  if (square) {
    layout.xaxis.constrain='domain';
    layout.yaxis.scaleanchor = 'x';
    layout.yaxis.scaleratio = ratio;
  }
  return {traces:traces,layout:layout};
}
window.procBuffer=procBuffer;
window.makeTraces=makeTraces;

