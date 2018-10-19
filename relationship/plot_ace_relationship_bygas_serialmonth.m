function [ ] = plot_ace_relationship_bygas_serialmonth( gas_in, years_in, do_plot )
%A funcion to plot the relationships between ace gases for a given year.

% *INPUT*
%           gas_in: STRING - the name of the gas for which you want to
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
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\\MAESTRO\climdata\';
% newclim_dir = '/Users/niall/Dropbox/climatology/nryan/reldata/';
% newclim_dir = 'C:\Users\ryann\Dropbox\climatology\nryan\reldata\';
newclim_dir = 'C:\Users\ryann\ACE\reldata\';

if nargin > 2
    yplot = do_plot;
else
    yplot = 1;
end
gasname = gas_in;
yearsin = years_in;
newfile_pre = strcat('ACEFTS_REL_v3_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_12.mat
% newfile_pre = strcat('ACEMAESTRO_CLIM_v1_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_12.mat

newfile_post = '.mat';
% newmonthnames = {'DJF', 'MAM', 'JJA', 'SON'};
months = 1:12;
slopej = nan(2,length(yearsin)*12);
sdatej = nan(1,length(yearsin)*12);
j = 0;

%% loop through the months and make the plots if you want
if yplot == 1
    % Do for the zonal vmr
    fs = 10;
    figi = randi(100);
    for y = 1:length(yearsin)
        figi = figi+1;
        %     figure(figi), set(gcf,'Position', [5,12,1096,704])
        %     figure(figi), set(gcf,'Position', [358,61,722,532])
%         figure(figi), set(gcf,'Position', [97,49,852,630]) %
        figure(figi), set(gcf,'Position', [97    49   713   898])
        suptitle(sprintf('%s relationship, %i',gasname, yearsin(y)))
        hold on
        cmin1 = 20;
        cmax1 = 70; % the maximum ace altitude used
        for i = 1:12
            j = j+1;
            sdatej(j) = datenum([yearsin(y),i,1]);
            filenewi = strcat( newclim_dir, gasname,'/','serial_month','/', newfile_pre, sprintf('%i_%0.2i',yearsin(y), months(i)), newfile_post);
            if exist(filenewi,'file') ~= 2
                fprintf('There is no file for %i_%0.2i. Moving on...\n', yearsin(y), months(i))
            else
                fprintf('\nloading %s\n',filenewi)
                reltan = load(filenewi); reltan = reltan.relstruct; % the variable is called climstruct in the new data
                slopej(1,j) = reltan.slope;
                slopej(2,j) = reltan.slope_error;
                
                subplot(4,3,i), plot_ace_relationship(reltan)
                caxis([cmin1,cmax1]);
                colorbar('delete')
                title ''
                set(gca,'FontSize',fs)
                if i ~= 1 && i ~= 4 && i ~= 7 && i ~= 10
                    ylabel('')
                    %                 set(gca,'yticklabel',{[]})
                end
                if i ~= 10 && i ~= 11 && i ~= 12
                    xlabel('')
                    %                 set(gca,'xticklabel',{[]})
                end
                if i == 12
                    c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]); % [0.923 0.309 0.016 0.35]
                    title(c,'altitude [km]','Position', [8 -23 0])
                end
            end
        end
    end
    %% plot the slopes
%     nansum(slopej(1,:))
    if nansum(slopej(1,:)) ~= 0
        figure(99)
        errorbar(sdatej,slopej(1,:),slopej(2,:),'Linewidth', 2), hold on;
        dynamicDateTicks([],[],'mm/yy')
        title(sprintf('slope of %s' ,gasname))
        xlabel('date')
        ylabel('slope [vmr/vmr]')
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

