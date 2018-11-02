function [ ] = plot_ace_relationship_bygas_year( gas_in, years_in, do_plot )
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
slopej = nan(2,length(yearsin));
sdatej = nan(1,length(yearsin));

%% loop through the months and make the plots if you want
if yplot == 1
    % Do for the zonal vmr
    fs = 10;
    figi = randi(100);
    %     figure(figi), set(gcf,'Position', [5,12,1096,704])
    %     figure(figi), set(gcf,'Position', [358,61,722,532])
    %         figure(figi), set(gcf,'Position', [97,49,852,630]) %
%     figure(figi), set(gcf,'Position', [97    49   969   898])
    figure(figi), set(gcf,'Position', [97  -211   969   895])
    suptitle(sprintf('%s relationship',gasname))
    hold on
    cmin1 = -90;
    cmax1 = -50; % the maximum ace altitude used
%     cmin1 = -90;
%     cmax1 = -40; %
    for i = 1:length(yearsin)
        sdatej(i) = datenum([yearsin(i),1,1]);
        filenewi = strcat( newclim_dir, gasname,'/', newfile_pre, sprintf('%i',yearsin(i)), newfile_post);
        if exist(filenewi,'file') ~= 2
            fprintf('There is no file for %i. Moving on...\n', yearsin(i))
        else
            fprintf('\nloading %s\n',filenewi)
            reltan = load(filenewi); reltan = reltan.relstruct; % the variable is called climstruct in the new data
            slopej(1,i) = reltan.slope;
            slopej(2,i) = reltan.slope_error;
            
            subplot(4,4,i), plot_ace_relationship(reltan)
            caxis([cmin1,cmax1]);
            colorbar('delete')
            title(sprintf('%i', years_in(i)))
            set(gca,'FontSize',fs)
            if i ~= 1 && i ~= 5 && i ~= 9 && i ~= 13
                ylabel('')
                %                 set(gca,'yticklabel',{[]})
            end
            if i ~= 13 && i ~= 14 && i ~= 15 && i ~= 16
                xlabel('')
                %                 set(gca,'xticklabel',{[]})
            end
            if i == 14
                c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]); % [0.923 0.309 0.016 0.35]
                title(c,'altitude [km]','Position', [8 -23 0])
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
end
%
end

