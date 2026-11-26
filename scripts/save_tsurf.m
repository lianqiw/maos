%Convert a tilted surface OPD array into .bin file for MAOS
%
%   Script that reads a tilted surface map (e.g. M3 print through) and its
%   parameters from FNIN and writes them to FNOUT with keywords in the
%   header, for use with ncpa.tsurf=['FNOUT'] in MAOS (see ncpa.tsurf in
%   config/maos/sim.conf and rmap_parse_keywords in lib/map.c). MAOS projects
%   it to the pupil with m3proj.
%
%   Variables used if they exist in the workspace:
%     fnin   Input file (default 'M3_PrintThrough_ZA45_AA180_YA45.bin.gz').
%            It must contain a cell with
%            a{1}=[dx, dy, ox, oy, txdeg, tydeg, ftel, fexit, fsurf] and
%            a{2}=the surface map.
%     fnout  Output file (default 'temp.bin').
%
%   Header keys written (all are required by MAOS; lengths in meters):
%     dx, dy        sampling along x and y,
%     ox, oy        coordinate of the first point (lower left corner),
%     txdeg, tydeg  x and y tilt angle of the surface in degree wrt the beam
%                   (90 is perpendicular),
%     ftel          telescope effective focal length,
%     fexit         distance between the exit pupil and the focus,
%     fsurf         distance between the surface (M3) and the focus.
%   The map is the surface (NOT the OPD) in meters, with the first index
%   along x and the second along y.
%
%   Uses the read/write mex functions if available, otherwise readbin and
%   writebin.

if ~exist('fnin','var')
    fnin='M3_PrintThrough_ZA45_AA180_YA45.bin.gz';
end
if ~exist('fnout','var')
    fnout='temp.bin';
end
if exist('read')==3 %#ok<EXIST> mex function
    a=read(fnin);
else
    a=readbin(fnin);
end
if ~iscell(a) || numel(a)<2 || numel(a{1})<9
    error('save_tsurf: %s must contain a cell {[dx, dy, ox, oy, txdeg, tydeg, ftel, fexit, fsurf], map}', fnin);
end
%%First set the OPD of the surface and other parameters
OPD=a{2};
dx=a{1}(1);
dy=a{1}(2);
ox=a{1}(3);
oy=a{1}(4);
txdeg=a{1}(5);
tydeg=a{1}(6);
ftel=a{1}(7);
fexit=a{1}(8);
fsurf=a{1}(9);
vx=0; %not used by tsurf
vy=0;

header=sprintf('ox=%20.15g\noy=%20.15g\ndx=%20.15g\ndy=%20.15g\ntxdeg=%20.15g\ntydeg=%20.15g\nftel=%20.15g\nfexit=%20.15g\nfsurf=%20.15g\n',...
               ox, oy, dx, dy, txdeg, tydeg, ftel, fexit, fsurf);

if exist('write')==3 %#ok<EXIST> mex function
    write(OPD,header,fnout);
else
    writebin(OPD,header,fnout);
end
