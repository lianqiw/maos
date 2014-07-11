%run_maos
%run the maos.mex in matlab and plot a few data
if ~exist('simu','var')
    simu=maos('setup','', 20);
end
parms=simu.parms;
recon=simu.recon;
aper=simu.aper;
powfs=simu.powfs;

im=1;


%TOMO output
sfigure(im); im=im+1; clf;
draw(recon.xloc, simu.opdr)

%DM Fitting output
sfigure(im); im=im+1; clf;
draw(recon.aloc, simu.dmfit)

%DM commands applied
sfigure(im); im=im+1; clf;
draw(recon.aloc, simu.dmreal)

%Gradients
sfigure(im); im=im+1; clf;
saloc=powfs(1).saloc;
gcl=simu.gradlastcl{1};
gol=simu.gradlastol{1};
subplot(1,2,1)
quiver(saloc(:,1), saloc(:,2), gcl(1:end/2), gcl(end/2+1:end));
axis square equal tight
title('CL grads');
subplot(1,2,2)
quiver(saloc(:,1), saloc(:,2), gol(1:end/2), gol(end/2+1:end));
axis square equal tight
title('PSOL grads');

%Subaperture images
sfigure(im); im=im+1; clf;
plot_saimage(saloc, 0.5, simu.ints{1}, 2);
