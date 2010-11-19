%Convert a tilted surface OPD array into .bin file for MAOS

a=read('M3_PrintThrough_ZA45_AA180_YA45.bin.gz');
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
vx=0;
vy=0;

header=sprintf('ox=%20g\noy=%20g\ndx=%20g\ndy=%20g\ntxdeg=%20g\ntydeg=%20g\nftel=%20g\nfexit=%20g\nfsurf=%20g\n',...
               ox, oy, dx, dy, txdeg, tydeg, ftel, fexit, fsurf);

write(OPD,header,'temp.bin.gz');
