%Convert a surface OPD array into .bin file for MAOS
%
%   Script that reads an OPD map and its parameters from FNIN and writes them
%   to FNOUT as a map with keywords in the header, for use with
%   ncpa.surf=['FNOUT'] in MAOS (see ncpa.surf in config/maos/sim.conf and
%   map_parse_keywords in lib/map.c).
%
%   Variables used if they exist in the workspace:
%     fnin   Input file (default 'KECK'). It must contain a cell with
%            a{1}=[dx, dy, ox, oy, h, ...] (dy is not used; MAOS takes dy=dx)
%            and a{2}=OPD.
%     fnout  Output file (default 'temp.bin').
%
%   Header keys written (all in meters, m/s for vx and vy):
%     ox, oy  coordinate of OPD(1,1), the lower left corner,
%     dx      sampling (dy defaults to dx),
%     h       height of the surface conjugation from the pupil,
%     vx, vy  frozen flow speed (0 here).
%   OPD is in meters, with the first index along x and the second along y
%   (OPD(ix,iy)). It should cover about 1/2 meter beyond the meta pupil on
%   each side to avoid edge roll off in interpolation.
%
%   Uses the read/write mex functions if available, otherwise readbin and
%   writebin.

if ~exist('fnin','var')
    fnin='KECK';
end
if ~exist('fnout','var')
    fnout='temp.bin';
end
if exist('read')==3 %#ok<EXIST> mex function
    a=read(fnin);
else
    a=readbin(fnin);
end
if ~iscell(a) || numel(a)<2 || numel(a{1})<5
    error('save_surf: %s must contain a cell {[dx, dy, ox, oy, h, ...], OPD}', fnin);
end
%%First set the OPD of the surface and other parameters
OPD=a{2};
dx=a{1}(1);
ox=a{1}(3);
oy=a{1}(4);
h=a{1}(5);
vx=0;
vy=0;

header=sprintf('ox=%20.15g\noy=%20.15g\ndx=%20.15g\nh=%20.15g\nvx=%20.15g\nvy=%20.15g\n',...
               ox, oy, dx, h, vx, vy);

if exist('write')==3 %#ok<EXIST> mex function
    write(OPD,header,fnout);
else
    writebin(OPD,header,fnout);
end
