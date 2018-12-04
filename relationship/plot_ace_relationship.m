function [ ] = plot_ace_relationship( tanstruct_relation_in )
%A function to plot the relationship between two gases measured by ace. The
%function uses a scatter plot.
%
% *INPUT*
%           tanstruct_relation_in: STRUCTURE - the data to be plotted. The
%           structure is created by running 'get_ace_gas_relationship.m'.
%
% *OUTPUT*
%           makes a plot 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
tanrel = tanstruct_relation_in;
norbit = length(tanrel.occultation);
nalt = length(tanrel.altitude_km(:,1));
%% reshape the data into single columns
column_size = [norbit.*nalt, 1];
vmrcolx = reshape(tanrel.vmr_x, column_size);
% vmrerrorcolx = reshape(tanrel.vmr_error_x, column_size);
vmrcoly = reshape(tanrel.vmr_y, column_size);
% vmrerrorcoly = reshape(tanrel.vmr_error_y, column_size);
if length(tanrel.altitude_km(1,:)) == 1
    altcol = repmat(tanrel.altitude_km,[norbit,1]);
else
    altcol = reshape(tanrel.altitude_km, column_size);
end
if isfield(tanrel,'lat')
    latcol = reshape(tanrel.lat, column_size);   
end

%% get the line of best fit from slope and intercept
fittedX = linspace(min(vmrcolx), max(vmrcolx), 200);
fittedY = polyval([tanrel.slope, tanrel.intercept], fittedX);

%% plot the data
lw = 2;
fs = 16;
fs_text = 10;
ms = 2;
% l2 = sprintf('m = %0.2e%s%0.2f, c = %0.2e%s%0.2f', tanrel.slope, char(177), tanrel.slope_error, tanrel.intercept, char(177), tanrel.intercept_error);
% l1 = sprintf('r = %0.2f', round2(tanrel.correlation, 0.01));
l1 = sprintf('r = %0.2f\nm = %0.2e%s%0.2e\nc = %0.2e%s%0.2e', round2(tanrel.correlation, 0.01), tanrel.slope, char(177), tanrel.slope_error, tanrel.intercept, char(177), tanrel.intercept_error);

h1 = scatter(vmrcolx, vmrcoly, ms, latcol, 'filled'); hold on; %#ok<NASGU>
h2 = plot(fittedX, fittedY, 'k-', 'LineWidth', lw); %#ok<NASGU>
% labels and text
text(0.5*min(vmrcolx), 0.7*max(vmrcoly), l1, 'Fontsize', fs_text)
title(sprintf('%s vs. %s', tanrel.gas_y, tanrel.gas_x));
% legend(h2, l2, 'Location', 'best')
% legend(h2,l2)
xlabel('VMR')
ylabel('VMR')
%colobar for altitude
c = colorbar;
% title(c,'altitude [km]','Position', [12 -23 0])
title(c,'latitude [degN]','Position', [25 344 0])
set(gca,'FontSize',fs)
xlim([min(vmrcolx)*1.1, max(vmrcolx)*1.1])
ylim([min(vmrcoly)*1.1, max(vmrcoly)*1.1])
%
end

