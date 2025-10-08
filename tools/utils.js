"use strict"; //every variable need to be defined
function get_hostname() {
  /*
  * We open the websocket encrypted if this page came on an
  * https:// url itself, otherwise unencrypted
  */
  var u = window.location.href;
  if (u.substring(0, 5) === "https") {
    window.pcol = "wss://";
    u = u.substring(8);
  } else {
    window.pcol = "ws://";
    if (u.substring(0, 4) === "http") {
      u = u.substring(7);
    }
  }
  u = u.split('/')[0];
  const p = u.indexOf(':');
  if (p != -1) {
    port = u.substring(p);
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
function Menu({ label, child }){
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
    <div className="menu-container" >
      <div onClick={handleState} className="menu-button" ref={menuRef}>
        {label}
      </div>
      {isOpen && child}
    </div>
  );
}
//Use global function instead of import to avoid error with in-browser babel transformer
window.get_hostname=get_hostname;
window.split_hostname=split_hostname;
window.Menu=Menu;
window.Progress=Progress;
