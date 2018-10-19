function [ ] = plot_ace_climatology_bygas_3monthly( gas_in, do_plot )
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
newclim_dir = 'C:\Users\ryann\ACE\climdata\';

vmrzon = cell(1,4);
if nargin > 1
    yplot = do_plot;
else
    yplot = 1;
end
gasname = gas_in;
newfile_pre = strcat('ACEFTS_CLIM_v3_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_12.mat
newfile_post = '.mat';
newmonthnames = {'DJF', 'MAM', 'JJA', 'SON'};


%% loop through the months and fill in the cells of arrays
for i = 1:4 % do all 3-month periods
    filenewi = strcat( newclim_dir, gasname,'/', newfile_pre, newmonthnames{i}, newfile_post);
    clim = load(filenewi); clim = clim.climstruct; % the variable is called climstruct in the new data
    vmrzon{i} = clim.vmr_zonal;
end
%% Make the plots if you want
if yplot == 1
    % Do for the zonal vmr
    figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
    figure(figi), set(gcf,'Position', [358,61,722,532])
    suptitle(sprintf('%s zonal VMR, 2004 - 2017',gasname))
    hold on
    cmax_all = nan(1,4);
    for i = 1:4
        subplot(2,2,i), plot_on_ace_clim_latlev( vmrzon{i}, newmonthnames{i} )
        %         cmax1 = 1e-6; % may need to be gas specific
        cmax_all(i) = max(max(abs(vmrzon{i})));
        if i == 4
            c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
            cmax1 = max(cmax_all);
            if ~isnan(cmax1)
                caxis([0,cmax1]);
            end
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
    
end
%
end

