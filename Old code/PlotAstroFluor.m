function astroFluorFig = PlotAstroFluor( T, Fastro, varargin )

IP = inputParser;
addRequired( IP, 'T', @iscell )
addRequired( IP, 'Fastro', @isstruct )
addOptional( IP, 'Fback', {}, @iscell )
addOptional( IP, 'velocity', {}, @iscell )
addParameter( IP, 'unit', 's', @ischar);
addParameter( IP, 'title', '', @ischar);
parse( IP, T, Fastro, varargin{:} ); 
Fback = IP.Results.Fback; 
velocity = IP.Results.velocity; 
timeUnit = IP.Results.unit;
titleString = IP.Results.title;

Nmovie = numel(Fastro);
Ncell = size(Fastro(1).cell, 2 );
cellColor = distinguishable_colors(Ncell);
astroFluorFig = figure('WindowState','max');
if ~isempty(velocity)
    sp(2) = subplot(3,1,3); sp(1) = subplot(3,1,1:2); 
else
    sp(1) = subplot(1,1,1);
end
% Plot fluorescence
for m = 1:Nmovie
    if ~isempty(Fback)
        plot( T{m}, Fback{m}, 'k'); hold on;
    end
    for c = 1:Ncell
        plot( T{m}, Fastro(m).cell(:,c), 'color', cellColor(c,:)); hold on;
    end
end
xlim([-Inf,Inf]); 
ylabel('Fluorescence');
title( sprintf('%s', titleString), 'interpreter','none' );
% Plot wheel velocity
if ~any( cellfun( @isempty, velocity ) ) %~isempty(velocity)
    subplot(sp(2));
    for m = 1:Nmovie
        plot( T{m}, velocity{m}, 'k'); hold on;
    end
    ylabel('Wheel Velocity (cm/s)');
end
xlabel(sprintf('Time (%s)', timeUnit));
linkaxes(sp,'x');

%{
if iscell(T)

else
    h(1) = plot( T, Fback, 'k' ); hold on;
    h(2) = plot( T, Fastro.cell, 'b' ); 
    h(3) = plot( T, Fastro.soma, 'r' ); 
    h(4) = plot( T, Fastro.proc, 'g' ); 
    legend(h, {'Background','Cell','Soma','Processes'});
end
%}

end