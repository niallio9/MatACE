function [ ] = plot_ace_climatology_bygas_serialmonth( gas_in, year_in, do_log, do_obs_limit )
%A funcion to compare the climatology made by Jaho, with the current
%version of the climatology. The assumption is that the two versions are
%made on the same latitude and altitude grid, and ar for the same gas.

% *INPUT*
%           gasname: STRING - the name of the gas for which you want to
%           compare the climatologies.
%
%           filename_oldclim: STRING - the netcdf file that contains the
%           old version of the climatology data.
%
% *OUTPUT*
%           vmrzon_dif: CELL OF ARRAYS - the differnce of the new and old
%           zonal vmrs for each month of the year.
%
%           vmrzonvar_dif: CELL OF ARRAYS - the differnce of the new and
%           old standard deviation of the zonal vmrs for each month of the
%           year.
%
%           obscount_dif: CELL OF ARRAYS - the differnce of the new and old
%           observation counts for each month of the year.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *INPUT*
%           climstruct_in: STRUCTURE - contains the gas specific ACE
%           climatology data. This structure can be created with
%           'make_ace_climatology.m' or with
%           'make_ace_climatology_multiple.m'.
%
% *OUTPUT*
%           oldclim_in: netcdf - contains the gas specific climatology data
%           of the previous version of the climatology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
% newclim_dir = 'C:\Users\ryann\ACE\climdata_testing\';
% newclim_dir = 'C:\Users\ryann\ACE\climdata\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
newclim_dir = 'C:\Users\ryann\ACE\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\MAESTRO\climdata\';

vmrzon = nan(48,36,12);
if nargin > 2
    dolog = do_log;
else
    dolog = 0;
end

if nargin > 3
    do_obs_lim = do_obs_limit;
else
    do_obs_lim = 1;
end
gasname = gas_in;
yearin = year_in;
newfile_pre = strcat('ACEFTS_CLIM_v3_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_12.mat
% newfile_pre = strcat('ACEMAESTRO_CLIM_v3_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_12.mat

newfile_post = '.mat';
% newmonthnames = {'DJF', 'MAM', 'JJA', 'SON'};
months = 1:12;



%% loop through the months and fill in the cells of arrays
for i = 1:12 % do all months
    filenewi = strcat( newclim_dir, gasname,'/','serial_month','/', newfile_pre, sprintf('%i_%0.2i',yearin, months(i)), newfile_post)

    if exist(filenewi,'file') ~= 2
        fprintf('There is no file for %i_%0.2i. Moving on...\n', yearin, months(i))
        vmrzon(:,:,i) = nan(48,36); % this is hardcoded here, for now.
    else
        clim = load(filenewi); clim = clim.climstruct; % the variable is called climstruct in the new data
        if do_obs_lim == 1
            clim = reduce_climstruct_data_by_obs_nr(clim,5);
        end
        vmrzon(:,:,i) = clim.vmr_zonal;
    end
end

vmrzon = vmrzon*1e6; % put into ppm
%% If plotting in log
if dolog == 1
    cmax = 2000;
    c = [ 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]; % for tick marks
    cstart = 4;
    c = c(cstart:end);
    c = round2(cmax * c, 1);
%     c = [2.5 3 3.5 4 4.5 5 5.5 6 7 8 20 50 5000];
    vmrzon(vmrzon <= 0) = nan;%c(1)*cmax;
    vmrzon = log10(vmrzon);
    cgrid = log10(c); % for caxis
    ymax = 300;
    ymin = 60;
%     c = [0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 1]; % for tick marks
%     cstart = 4;
%     c = c(cstart:end);
%     c = cmax * c;
%     vmrzon(vmrzon <= 0) = nan;%c(1)*cmax;
%     vmrzon = log10(vmrzon);
%     cgrid = log10(c); % for caxis
else
    vmrzon_short = vmrzon(6:33,:);
    cmax = max(vmrzon_short(:))
    cgrid = linspace(0, cmax, 10);
    c = cgrid; 
%     ymax = 500;
%     ymin = 0.1;
    ymax = 300;
    ymin = 0.1;
end

% Do for the zonal vmr
figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     figure(figi), set(gcf,'Position', [358,61,722,532])
figure(figi), set(gcf,'Position', [97,49,852,630])
suptitle(sprintf('%s zonal VMR, %i',gasname, yearin))
hold on
for i = 1:12
    subplot(4,3,i), plot_on_ace_clim_latlev_contourf( vmrzon(:,:,i), sprintf('%0.2i',months(i)), cgrid )
%     subplot(4,3,i), plot_on_ace_clim_latlev_pcolor( vmrzon(:,:,i), sprintf('%0.2i',months(i)) )
    ylim([10^0,10^3])
    caxis([cgrid(1),cgrid(end)])
    colormap(parula(length(cgrid) - 1))
    ylim([ymin, ymax])
    if i ~= 1 && i ~= 4 && i ~= 7 && i ~= 10
        ylabel('')
        set(gca,'yticklabel',{[]})
    end
    if i ~= 10 && i ~= 11 && i ~= 12
        xlabel('')
        set(gca,'xticklabel',{[]})
    end
    if i == 12
        c = colorbar('YTick', cgrid, 'Yticklabel', c, 'location','eastoutside','position',[0.923 0.309 0.016 0.35]);
        title(c,'VMR','Position', [8 -23 0])
    end
end
%     figi = figi + 1;
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     suptitle(sprintf('percent difference in %s zonal VMR',gasname))
%     hold on
%     for i = 1:4
%         subplot(2,2,i), plot_on_ace_clim_latlev( vmrzon_dif_percent{i}, newmonthnames{i} )
%         cmax11 = 100;
% %         cmax11 = max(max(abs(vmrzon_dif_percent{i})));
%         caxis([-cmax11,cmax11]);
%         if i == 4
%             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             title(c,'%','Position', [8 -23 0])
%         end
%     end
%     % Do for the zonal vmr variance
%     figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     suptitle(sprintf('difference in %s zonal VMR variance',gasname))
%     hold on
%     for i = 1:4
%         subplot(2,2,i), plot_on_ace_clim_latlev( vmrzonvar_dif{i}, newmonthnames{i}(1:3) )
%         cmax2 = max(max(abs(vmrzonvar_dif{i})));
%         caxis([-cmax2,cmax2]);
%         if i == 4
%             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             title(c,'VMR','Position', [8 -23 0])
%         end
%     end
%     % Do for the observation count
%
%     figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     suptitle(sprintf('difference in %s observation count',gasname))
%     hold on
%     for i = 1:4
%         subplot(2,2,i), plot_on_ace_clim_latlev( obscount_dif{i}, newmonthnames{i} )
%         cmax3 = max(max(abs(obscount_dif{i})));
%         caxis([-cmax3,cmax3]);
%         if i == 4
%             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             title(c,'VMR','Position', [8 -23 0])
%         end
%     end
%
%     % Plot a mean profile for reference
%     figi = randi(100);
%     for i = 1:4
%         figure(figi), plot(vmr_mean{i}(:,1), vmr_mean{i}(:,2), 'LineWidth', 2)
%         hold on
%         ytickspace = [10^-4 10^-2 1 10^2];
%         yticks(ytickspace);
%         xlabel('VMR')
%         ylabel('altitude [hPa]')
%         set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on');
%         set(gca,'FontSize',16)
%     end
%     legend(string(newmonthnames))

%
end

