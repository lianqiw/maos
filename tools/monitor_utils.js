"use strict"; //every variable need to be defined
function get_hostname() {
  /*
  * We open the websocket encrypted if this page came on an
  * https:// url itself, otherwise unencrypted
  */
  var u = window.location.href;
  if (u.substring(0, 5) === "https") {
    window.pcol = "wss://";
    window.hcol = "https://";
    u = u.substring(8);
  } else {
    window.pcol = "ws://";
    window.hcol = "http://"; 
    if (u.substring(0, 4) === "http") {
      u = u.substring(7);
    }
  }
  u = u.split('/')[0];
  const p = u.indexOf(':');
  if (p != -1) {
    window.port = u.substring(p);
  }
  return u;
}
function split_hostname(hostname) {
  let host = hostname.split('://').pop().split(':')[0];
  if (/^(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.(25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)$/.test(host)) {
    return host;
  } else {
    return host.split('.')[0];
  }
}
function splitText(text) {
    const result = [];
    let current = '';
    let bracketDepth = 0;

    for (let i = 0; i < text.length; i++) {
        const c = text[i];

        if (c === '[') {
            bracketDepth++;
            current += c;
        } else if (c === ']') {
            bracketDepth--;
            current += c;
        } else if (c === ' ' && bracketDepth === 0) {
            const prev = text.slice(0, i).trimEnd();
            const next = text.slice(i + 1);

            // Keep space around '='
            const aroundEquals =
                prev.endsWith('=') || next.startsWith('=');

            // Keep space after a single-letter dash option, e.g. "-c value"
            const afterDashOption = /-[A-Za-z]$/.test(prev);

            if (aroundEquals || afterDashOption) {
                current += c;
            } else {
                if (current.trim()) {
                    result.push(current.trim());
                    current = '';
                }
            }
        } else {
            current += c;
        }
    }

    if (current.trim()) {
        result.push(current.trim());
    }

    // Remove whitespace around '='
    return result.map(s => s.replace(/\s*=\s*/g, '='));
}
function Progress({ text, frac }){
    return (
      <td>
        <div className="progress" title={text}>
          <div className="progressbar" style={{ "width": frac }}></div>
          <span className="progressbar">{text}</span>
        </div>
      </td>
    );
}
function Menu({title, label, child }){
  const [isOpen, setIsOpen] = React.useState(false);
  //Ref to the entire menu container for click-outside detection
  const menuRef = React.useRef(null);
  // 2. Handle dismissal (click-outside or Escape key)
  React.useEffect(() => {
    const handleClickOutside = (e) => {
      // Check if the click was outside the menu itself
      if (menuRef.current && !menuRef.current.contains(e.target)) {
        setIsOpen(false);
      }
    };

    const handleEscape = (e) => {
      if (e.key === 'Escape') {
        setIsOpen(false);
      }
    };

    if (isOpen) {
      document.addEventListener('click', handleClickOutside);
      document.addEventListener('keydown', handleEscape);
    }
    return () => {//cleanup function that runs if dependencies is changed or if destroyed
      document.removeEventListener('click', handleClickOutside);
      document.removeEventListener('keydown', handleEscape);
    };
  }, [isOpen]);
  const handleState = () => setIsOpen(!isOpen) //avoid inline function that triggers rerender
  return (
    <div title={title} className="menu-container" >
      <div onClick={handleState} className="menu-button" ref={menuRef}>
        {label}
      </div>
      {isOpen && child}
    </div>
  );
}
function ResizableColumn({ name, width, onWidthChange }) {
  const startResize = (e) => {
    const startX = e.clientX;
    const startWidth = width;

    function onMove(e) {
      onWidthChange(Math.max(20, startWidth + e.clientX - startX));
    }

    function onUp() {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    }

    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  };

  return (
    <th title="Drag | to resize" style={{ width, position: "relative" }}>
      {name}
      <div
        onMouseDown={startResize}
        style={{
          position: "absolute",
          right: 0,
          top: 0,
          width: 6,
          height: "100%",
          cursor: "col-resize",
        }}
      />
    </th>
  );
}
function FloatInput({
  title='Input',
  name='Input',
  type="number",
  value,
  step,
  min,
  max,
  onChange,
  onApply,
}) {
  const inputRef = React.useRef(null);

  const handleSubmit = (e) => {
    e.preventDefault();
    if(onApply) {
      onApply(parseFloat(value || min))
    }
  };

  return (
    <li title={title}>
      <span>{name}</span>
      <div className="spring-spacer" />

      <form onSubmit={handleSubmit}>
        <input
          ref={inputRef}
          style={{ width: '4em' }}
          value={value}
          type={type}
          step={step}
          min={min}
          max={max}
          onClick={(e) => {
            if(inputRef.current) inputRef.current.select();
            e.stopPropagation();
          }}
          onChange={(e) => onChange(e.target.value)}
        />
      </form>
    </li>
  );
}
//Use global function instead of import to avoid error with in-browser babel transformer
window.get_hostname=get_hostname;
window.split_hostname=split_hostname;
window.splitText=splitText
window.Menu=Menu;
window.Progress=Progress;
window.ResizableColumn=ResizableColumn;
window.FloatInput=FloatInput;
