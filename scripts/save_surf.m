%Convert a surface OPD array into .bin file for MAOS

a=read('KECK');
%%First set the OPD of the surface and other parameters
OPD=a{2};
dx=a{1}(1);
ox=a{1}(3);
oy=a{1}(4);
h=a{1}(5);
vx=0;
vy=0;

header=sprintf('ox=%20g\noy=%20g\ndx=%20g\nh=%20g\nvx=%20g\nvy=%20g\n',...
               ox, oy, dx, h, vx, vy);

write(OPD,header,'temp.bin.gz');
