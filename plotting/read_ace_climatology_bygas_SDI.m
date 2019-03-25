function [ ] = read_ace_climatology_bygas_SDI( gas_in, years_in, do_plot )
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
%USER DEFINED
% climdir = '/Users/niall/Dropbox/climatology/nryan/climdata/'; % edit this to your directory that contains the ACE netcdf data
climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata_netcdf/';
% climdir = 'F:\ACE\climdata\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
% climdir = '/Users/niall/ACE/climdata_netcdf';
% climdir = 'C:\Users\ryann\ACE\climdata_netcdf';
% climdir = 'C:\Users\ryann\ACE\MAESTRO\climdata_netcdf';

if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end

if nargin > 2
    yplot = do_plot;
else
    yplot = 1;
end

gasname = gas_in;
yearin = years_in;
lyears = length(yearsin);

newfile = sprintf('SPARC_DI_T2Mz_%s_%i_ACEFTS_v3.6_i01.nc', gasname, yearin); % SPARC_DI_T2Mz_C2H2_2004_ACEFTS_v3.6_i01.nc
% newfile = sprintf('SPARC_DI_T2Mz_%s_%i_ACEMAESTRO_v3.13_i01.nc', gasname, yearin); % SPARC_DI_T2Mz_C2H2_2004_ACEFTS_v3.6_i01.nc
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

%% loop through the months and fill in the cells of arrays

filenewi = fullfile(climdir,gasname,newfile);
vmrzoni = ncread(filenewi,gasname);
% whos
% vmrzoni = vmrzoni(:,6:33,:);
vmrzoni(vmrzoni == -999) = nan;

%% Make the plots if you want
cmax_all = nan(1,12);
if yplot == 1
%     Do for the zonal vmr
    figi = randi(100);
    figure(figi), set(gcf,'Position', [5,12,1096,704])
    figure(figi), set(gcf,'Position', [71,-270,1112,954])
    
    suptitle(sprintf('%s zonal VMR for %i',gasname, yearin))
    hold on
    %     cmax1 = max(vmrzoni(:));
    cmax1 = 12e-9;
    for i = 1:12
        subplot(4,3,i), plot_on_ace_SDI_latlev( vmrzoni(:,:,i), monthnames{i}(1:3) );
        %         caxis([0,cmax1])
        cmax_all(i) = max(max(abs(vmrzoni(:,:,i))));
        caxis([0,cmax1]);
        if i == 12
            c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             cmax1 = max(cmax_all);
            if ~isnan(cmax1)
                caxis([0,cmax1]);
            end
            title(c,'VMR','Position', [8 -23 0]);
        end
    end
    ip = 0;
    figure(figi+1), set(gcf,'Position', [358,61,722,532])
    suptitle(sprintf('%s zonal VMR for %i',gasname, yearin))
%     cmax1 = max(vmrzoni(:));
    for i = 2:3:11
        ip = ip+1;
        subplot(2,2,ip), plot_on_ace_SDI_latlev( vmrzoni(:,:,i), monthnames{i}(1:3) );
%         caxis([0,cmax1])
        cmax_all(ip) = max(max(abs(vmrzoni(:,:,i))));
        if i == 11
            c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
            cmax1 = max(cmax_all);
            if ~isnan(cmax1)
                caxis([0,cmax1]);
            end
            title(c,'VMR','Position', [8 -23 0]);
        end
    end
%     % Do for the zonal vmr variance
%     figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     suptitle(sprintf('difference in %s zonal VMR variance',gasname))
%     hold on
%     for i = 1:12
%         subplot(3,4,i), plot_on_ace_SDI_latlev( vmrzonvar_dif(:,:,i), monthnames{i}(1:3) )
%         cmax11 = max(max(abs(vmrzonvar_dif(:,:,i))));
%         if ~isnan(cmax11)
%             caxis([-cmax11,cmax11]);
%         end
%         colorbar
%         %         if i == 12
%         %             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%         %             title(c,'VMR','Position', [8 -23 0])
%         %         end
%     end
%     % Do for the observation count
%     figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     suptitle(sprintf('difference in %s observation count',gasname))
%     hold on
%     for i = 1:12
%         subplot(3,4,i), plot_on_ace_SDI_latlev( obscount_dif(:,:,i), monthnames{i}(1:3) )
%         cmax3 = max(max(abs(obscount_dif(:,:,i))));
%         if ~isnan(cmax3)
%             caxis([-cmax3,cmax3]);
%         end
%         colorbar
%         %         if i == 12
%         %             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%         %             title(c,'VMR','Position', [8 -23 0])
%         %         end
%     end
%     
%     % Plot a mean profile for reference
%     figi = randi(100);
%     
%     figure(figi), plot_on_ace_SDI_latlev( vmr_mean, 'mean over year' )
%     cmax4 = max(max(vmr_mean));
%     if ~isnan(cmax4)
%         caxis([-cmax4,cmax4]);
%     end
%     colorbar
end

end

