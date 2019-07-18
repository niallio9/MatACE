function [ stats_out  ] = plot_ace_comparison_by_gas(lat_bounds, time_subset, gasname1, varargin )
%A funcion to compare a number of data sets with one other.

% *INPUT*
%           lat_bounds: VECTOR - the start and end latitudes within which
%           you want to subset the data. Can be left as an empty vector.
%
%           time_subset: STRING - to seasonally subset the data. Options
%                        are 'DJF, 'MAM', 'JJA', 'SON'.
%                        SCALAR - The subset the data by month of the year.
%                        Data will be used if it occurs in that month.
%                        Options are 1 to 12. Can be left as an empty
%                        vector. 
%                        VECTOR - The start and end years within which you
%                        want to subset the data. Can be left as an empty
%                        vector. 
%                        
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

matdirectory = 'C:\Users\ryann\ACE\matdata\';
%% USER DEFINED
iplot_alt = flip(15: 3: 27); % the indices of the altitudes at which you want to plot the timeseries
nalt = length(iplot_alt);

%%STANDARD
filein_pre = 'ACE_v3p6_';
filein_post = '.mat';
figi = randi(100);
figii = figi + 1;
gasnameN = varargin;
ngas = length(gasnameN);
stats_out = cell(1, ngas);

%%
filein_1 = strcat(matdirectory,filein_pre, gasname1, filein_post);
gas1 = load(filein_1); gas1 = gas1.tanstruct;
if length(gas1.gas) > 6
    if strcmp(gas1.gas(end - 6: end), 'maestro')
        gas1 = apply_ace_filter_mad3(gas1);
    else
        gas1 = apply_ace_flags(gas1);
    end
else
    gas1 = apply_ace_flags(gas1);
end
% if you want to subset the data by lat or time
if ~isempty(lat_bounds)
    gas1 = subset_ace_by_lat_tangent(gas1, lat_bounds(1), lat_bounds(2));
end
if ~isempty(time_subset)
    if ischar(time_subset)
        gas1 = subset_ace_by_3month(gas1, time_subset);
    elseif isscalar(time_subset)
        gas1 = subset_ace_by_month(gas1, time_subset);
    else
        gas1 = subset_ace_by_date(gas1, [time_subset(1), 0, 0], [time_subset(2), 12, 31]);
    end
end
zfix = gas1.altitude_km;
figpos = [238          55        1094         611];
% legend_text = {};
legend_text = cell(2*ngas, 1);

%% profile plots
for i = 1 : ngas
    filein_i = strcat(matdirectory,filein_pre, gasnameN{i}, filein_post)
    gas2 = load(filein_i); gas2 = gas2.tanstruct;
    if length(gas2.gas) > 6
        if strcmp(gas2.gas(end - 6: end), 'maestro')
            gas2 = apply_ace_filter_mad3(gas2);
        else
            gas2 = apply_ace_flags(gas2);
        end
    else
        gas2 = apply_ace_flags(gas2);
    end
    [gas1_match, gas2] = match_ace_data(gas1, gas2);
    fprintf('there are %i coinciding profiles\n', length(gas1_match.occultation))
    stats = compare_data_sets_2d( gas1_match.vmr, gas1_match.vmr_error, gas2.vmr, gas2.vmr_error );
    stats_out{i} = stats;
%     legend_text = [legend_text; gasnameN{i} ; gasname1]
    legend_text{2*i - 1} = gasnameN{i};
    legend_text{ 2*i} = gasname1;
    %
    lw = 2;
    fs = 16;
    numplots = 6;
    %
    figure(figi), colorOrder = get(gca, 'ColorOrder'); c1 = colorOrder(2*i - 1, :); c2 = colorOrder(2*i, :);
    figure(figi), set(gcf,'Position', figpos)
    ax(1) = subplot(1, numplots, 1:2); errorbar(stats.dataN_mean, zfix, stats.dataN_std, 'horizontal', 'LineWidth', lw,  'color', c1), hold on, %c = p.Color;
    subplot(1, numplots, 1:2); errorbar(stats.data1_mean, zfix, stats.data1_std, 'horizontal', 'LineWidth', lw,  'color', c2)
% % %     ax(1) = subplot(1, numplots, 1:2); errorbar(stats.dataN_mean, zfix, stats.dataN_mean_error, 'horizontal', 'LineWidth', lw,  'color', c1), hold on, %c = p.Color;
% % %     subplot(1, numplots, 1:2); errorbar(stats.data1_mean, zfix, stats.data1_mean_error, 'horizontal', 'LineWidth', lw,  'color', c2)   
    set(gca,'Fontsize',fs)
    xlabel('VMR','Fontsize',fs)
    ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','XMinorTick','on','color','none')
    grid minor
    if i == ngas
        legend(legend_text)
    end

    ax(2) = subplot(1, numplots, 3); errorbar(stats.difference_mean, zfix, stats.difference_std, 'horizontal', 'LineWidth', lw, 'color', c1), hold on
    set(gca,'Fontsize',fs)
    xlabel('\Delta VMR','Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    grid minor
    
    ax(3) = subplot(1, numplots, 4); errorbar(stats.pdifference_mean, zfix, stats.pdifference_std, 'horizontal', 'LineWidth', lw, 'color', c1), hold on
    set(gca,'Fontsize',fs)
    xlabel('\Delta VMR [%]','Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    grid minor
    xlim([-100 100])

    ax(4) = subplot(1, numplots, 5); plot(stats.correlation, zfix, 'Linewidth', lw, 'color', c1), hold on
    xlabel('correlation','Fontsize',fs)
    set(gca,'Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
%     linkaxes(ax,'y')
    grid minor
    xlim([0 1])
    
    ax(5) = subplot(1, numplots, 6);  errorbar(stats.regression, zfix, stats.regression_error, 'horizontal', 'LineWidth', lw, 'color', c1), hold on
    xlabel('slope','Fontsize',fs)
    set(gca,'Fontsize',fs)
    %ylabel('altitude [km]','Fontsize',fs)
    set(gca, 'YGrid', 'on', 'XGrid', 'on','YMinorTick','on','yticklabel',[],'XMinorTick','on','color','none')
    linkaxes(ax,'y')
    grid minor
    xlim([0 2])
%     legend() 
    
    
    %% timeseries plots
    sdates = mjd2datenum(gas1_match.date_mjd);
    axx = nan(nalt,1);
    figure(figii), colorOrder = get(gca, 'ColorOrder'); c1 = colorOrder(2*i - 1, :); c2 = colorOrder(2*i, :);
    figure(figii), set(gcf,'Position', figpos)
    for n = 1: nalt
        figure(figii), axx(n) = subplot(nalt,1,n); hold(axx(n),'on')
%         errorbar(axx(n),sdates, gas1_match.vmr(iplot_alt(n),:), gas1_match.vmr_error(iplot_alt(n),:), 'vertical', '.', 'color', c2 )
%         errorbar(axx(n),sdates, gas2.vmr(iplot_alt(n),:), gas2.vmr_error(iplot_alt(n),:), 'vertical','.', 'color', c1 )
        plot(axx(n),sdates, gas2.vmr(iplot_alt(n),:), '.', 'color', c1 )
        plot(axx(n),sdates, gas1_match.vmr(iplot_alt(n),:), '.', 'color', c2 )
        ytext = max(gas1_match.vmr(iplot_alt(n),:));
        text(sdates(1), ytext, sprintf('%0.1f km)', gas1_match.altitude_km(iplot_alt(n)) ));
        if i == ngas && n == 1
            legend(legend_text)
        end
%         nanstd(gas2.vmr(iplot_alt(n),:))
%         nanstd(gas1_match.vmr(iplot_alt(n),:))
        
        ylabel('VMR')
%         if n == 1
%             title(sprintf('Cly family zonal VMR, %i-%i%c latitude', latmin, latmax, char(176)))
%             %             legend('ClO', 'HOCl', 'HCl', 'ClONO2', 'Cly', 'Location', 'eastoutside' )
%         end
        if n == nalt
            %                 n
            dynamicDateTicks([],'x','mm')
            %                 dynamicDateTicks
            xlabel('year')
            linkaxes(axx,'x')
        else
            %                 datetick('x','mmyy')
            set(gca,'xticklabel',{[]})
        end
    end
    
end




%
end

