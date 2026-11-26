%run_maos
%run the maos.mex in matlab and plot a few data
%
%   Script that sets up MAOS inside Matlab with the maos mex function (built
%   from mex/maos.c), runs NSTEP time steps, and plots the tomography output,
%   the DM commands, the gradients and the subaperture images of WFS 1.
%
%   Variables used if they exist in the workspace:
%     simu   Result of a previous run. If it exists, the simulation is not
%            run again; clear it to rerun.
%     conf   Option string passed to maos('setup', conf), e.g.
%            '-c scao_lgs.conf -o out -g-1 sim.end=100'. Default is ''
%            (default.conf; nothing is saved without -o).
%     nstep  Number of time steps to run. Default is 20.
%
%   Requirements: the maos mex function, and the plotting helpers draw(loc,
%   opd) and plot_saimage(saloc, dsa, ints, scale), which are not part of
%   MAOS. sfigure is used if available, otherwise figure.
%
%   MAOS keeps its state between calls. To set up again (e.g. with a
%   different conf), call maos('reset') and clear simu first.
%
%   The interface of the maos mex function:
%     maos('setup', '-o dirout -n N -c scao_ngs.conf -g0')
%     maos('reset')
%     simu=maos('sim', nstep)
%     simu=maos('get','sim')
%     parms=maos('get','parms')

if ~exist('simu','var')
    if ~exist('conf','var')
        conf='';
    end
    if ~exist('nstep','var')
        nstep=20;
    end
    if exist('maos')~=3 %#ok<EXIST>
        error('run_maos: the maos mex function is not found. Build it in mex/ and add it to the path.');
    end
    %'setup' does not run any time step (its 3rd argument is not used) and
    %there is no simulation data yet, so run the steps with 'sim'.
    maos('setup', conf);
    simu=maos('sim', nstep);
end
if ~isstruct(simu) || ~isfield(simu, 'parms')
    error('run_maos: maos returned no simulation data. Check the maos output; clear simu and call maos(''reset'') to retry.');
end
for fun={'draw', 'plot_saimage'}
    if ~exist(fun{1}) %#ok<EXIST>
        error('run_maos: the plotting function %s is not found on the path.', fun{1});
    end
end
if exist('sfigure') %#ok<EXIST>
    newfig=@sfigure; %does not steal focus
else
    newfig=@figure;
end
parms=simu.parms;
recon=simu.recon;
aper=simu.aper;
powfs=simu.powfs;

im=1;


%TOMO output
newfig(im); im=im+1; clf;
draw(recon.xloc, simu.opdr)

%DM Fitting output
newfig(im); im=im+1; clf;
draw(recon.aloc, simu.dmrecon) %sim_t has no dmfit; dmrecon is the fitting output

%DM commands applied
newfig(im); im=im+1; clf;
draw(recon.aloc, simu.dmreal)

%Gradients
%Assumes WFS 1 belongs to powfs 1. Gradients are ordered as all x then all y.
newfig(im); im=im+1; clf;
saloc=powfs(1).saloc;
gcl=simu.gradlastcl{1};
gol=simu.gradlastol{1};
subplot(1,2,1)
quiver(saloc(:,1), saloc(:,2), gcl(1:end/2), gcl(end/2+1:end));
axis square equal tight
title('CL grads');
subplot(1,2,2)
if ~isempty(gol) %PSOL gradients exist only for some reconstructors
    quiver(saloc(:,1), saloc(:,2), gol(1:end/2), gol(end/2+1:end));
    axis square equal tight
end
title('PSOL grads');

%Subaperture images (only with physical optics WFS)
if ~isempty(simu.ints)
    newfig(im); im=im+1; clf;
    plot_saimage(saloc, 0.5, simu.ints{1}, 2);
end
