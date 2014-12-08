loc=read('accphi_loc');
locin=read('accphi_locin');
hl=read('accphi_hfor');
hc=read('accphi_cub_hfor');
map2pts=read('accphi_pts1');
loc2loc=read('accphi_loc2loc');
map2loc=read('accphi_loc0');
map2stat=read('accphi_stat');
h=read('accphi_loc2h');

cubic_h=read('accphi_cub_loc2h');
cubic_loc2loc=read('accphi_cub_loc2loc');
cubic_map2loc=read('accphi_cub_map2loc');
cubic_locmap2loc=read('accphi_cub_locmap2loc');


sfigure(1)
clf
subplot(2,3,1)
draw(loc, map2loc);
title('map2loc');
subplot(2,3,2)
draw(loc, map2pts-map2loc);
title('map2pts-map2loc');
subplot(2,3,3)
draw(loc, map2stat-map2loc);
title('map2stat-map2loc');
subplot(2,3,4)
draw(loc, h-map2loc);
title('h-map2loc');
subplot(2,3,5)
draw(loc, loc2loc-map2loc);
title('loc2loc-map2loc');
subplot(2,3,6)
draw(loc, loc2loc);
title('loc2loc');

sfigure(2)
clf
subplot(2,3,1)
draw(loc, cubic_loc2loc-map2loc);
title('cubic loc2loc-linear map2loc');
subplot(2,3,2)
draw(loc, cubic_map2loc-cubic_loc2loc);
title('cubic map2loc-cubic loc2loc');
subplot(2,3,3)
draw(loc, cubic_locmap2loc-cubic_loc2loc);
title('cubic locmap2loc-cubic loc2loc');
subplot(2,3,4)
draw(loc, cubic_h-cubic_loc2loc);
title('cubic h- cubic loc2loc');
subplot(2,3,5)
draw(loc, cubic_h-h)
title('cubic h- linear h');
subplot(2,3,6)
draw(loc, cubic_loc2loc);
title('cubic loc2loc');