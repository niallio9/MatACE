function [ stats_out  ] = plot_ace_comparison( tanstruct1, varargin )
%A funcion to compare a number of data sets with one other.

% *INPUT*
%           gas1: STRING - the name of the gas with which you want to
%           comapre the others
%
%           varargin: STRING - The other inputs are the data STRINGs of the
%           names of the gases with which to compare data1. e.g., 'gas2,
%           ..., gasN'.
%
% *OUTPUT*
%           stat_out: STRUCTURE - the statistics of the data comparison.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% matdirectory = 'C:\Users\ryann\ACE\matdata\';
%%STANDARD
% filein_pre = 'ACE_v3p6_';
% filein_post = '.mat';
figi = randi(100);

tanstructN = varargin;
ngas = length(tanstructN);
stats_out = cell(1, ngas);
% do_plot = varargin{end};
% if ~isscalar(do_plot)
%     error('the last element of varargin must be a scalar value of 0 or 1')
% elseif do_plot ~= 0 && do_plot ~= 1
%     error('the last element of varargin must be a scalar value of 0 or 1')
% end

gas1 = tanstruct1; gas1 = apply_ace_flags(gas1);
gasname1 = gas1.gas;
zfix = gas1.altitude_km;
figpos = [238          55        1094         611];
legend_text = cell(2*ngas, 1);

for i = 1 : ngas
    gas2 = tanstructN{i}; gas2 = apply_ace_flags(gas2);
    gasname2 = gas2.gas;
    [gas1_match, gas2] = match_ace_data(gas1, gas2);
    stats = compare_data_sets_2d( gas1_match.vmr, gas1_match.vmr_error, gas2.vmr, gas2.vmr_error );
    stats_out{i} = stats;
    legend_text{2*i - 1} = gasname2;
    legend_text{2*i} = gasname1;
    
    % %% Profile Plots
    %
    lw = 2;
    fs = 16;
    numplots = 6;
    %
    figure(figi), colorOrder = get(gca, 'ColorOrder'); c = colorOrder(i, :);
    figure(figi), set(gcf,'Position', figpos)
    ax(1) = subplot(1, numplots, 1:2); plot(stats.dataN_mean, zfix, 'LineWidth', lw,  'color', c), hold on, %c = p.Color;
    subplot(1, numplots, 1:2); plot(stats.data1_mean, zfix, 'LineStyle', '--', 'LineWidth', lw,  'color', c)
    set(gca,'Fontsize',fs)
    xlabel('VMR','Fontsize',fs)
    ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','XMinorTick','on','color','none')
    grid minor
    if i == ngas
        legend(legend_text)
    end

    ax(2) = subplot(1, numplots, 3); errorbar(stats.difference_mean, zfix, stats.difference_std, 'horizontal', 'LineWidth', lw, 'color', c), hold on
    set(gca,'Fontsize',fs)
    xlabel('\Delta VMR','Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    grid minor
    
    ax(3) = subplot(1, numplots, 4); errorbar(stats.pdifference_mean, zfix, stats.pdifference_std, 'horizontal', 'LineWidth', lw, 'color', c), hold on
    set(gca,'Fontsize',fs)
    xlabel('\Delta VMR [%]','Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    grid minor
    xlim([-100 100])

    ax(4) = subplot(1, numplots, 5); plot(stats.correlation, zfix, 'Linewidth', lw, 'color', c), hold on
    xlabel('correlation','Fontsize',fs)
    set(gca,'Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
%     linkaxes(ax,'y')
    grid minor
    xlim([0 1])
    
    ax(5) = subplot(1, numplots, 6);  errorbar(stats.regression, zfix, stats.regression_error, 'horizontal', 'LineWidth', lw, 'color', c), hold on
    xlabel('slope','Fontsize',fs)
    set(gca,'Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    linkaxes(ax,'y')
    grid minor
    xlim([0 2])
    legend()
    
end
%
end

