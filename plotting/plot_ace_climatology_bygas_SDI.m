function [ ] = plot_ace_climatology_bygas_SDI( gas_in, year_in, do_log )
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
% climdir = '/Volumes/Seagate Backup Plus Drive/ACE/climdata_netcdf/';
% climdir = 'F:\ACE\climdata\';
% climdir = '/net/deluge/pb_1/users/nryan/ACE/climdata/';
% climdir = '/Users/niall/ACE/climdata_netcdf';
climdir = 'C:\Users\ryann\ACE\ACE_climdata_netcdf_SDI_20190409';
% climdir = 'C:\Users\ryann\ACE\MAESTRO\climdata_netcdf';

if ~isdir(climdir)
    fprintf('\nIt doesn''t look like ''%s'' exists...\n',climdir)
    error('The directory containing the .mat climatology data couldn''t be found')
end
if nargin > 2
    dolog = do_log;
else
    dolog = 0;
end
gasname = gas_in;
yearin = year_in;

newfile = sprintf('SPARC_DI_T2Mz_%s_%i_ACEFTS_v3.6_i01.nc', gasname, yearin); % SPARC_DI_T2Mz_C2H2_2004_ACEFTS_v3.6_i01.nc
% newfile = sprintf('SPARC_DI_T2Mz_%s_%i_ACEMAESTRO_v3.13_i01.nc', gasname, yearin); % SPARC_DI_T2Mz_C2H2_2004_ACEFTS_v3.6_i01.nc
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};
months = 1:12;

%% get the name of the gas for the input folder name
if length(gasname) > 10 && strcmp(gasname(end-8:end-6), 'sap') % for '_sap_s10am' and '_sap_s10pm', etc.
    gasname_short = gasname(1:end-10);
elseif length(gasname) > 7 && strcmp(gasname(end-5:end), 'sap_am') % for '_sap_am' and '_sap_pm'
    gasname_short = gasname(1:end-7);
elseif length(gasname) > 7 && strcmp(gasname(end-5:end), 'sap_pm') % for '_sap_am' and '_sap_pm'
    gasname_short = gasname(1:end-7);
elseif length(gasname) > 4 && strcmp(gasname(end-2:end), 'sap') % for '_sap'
    gasname_short = gasname(1:end-4);
elseif length(gasname) > 3 && strcmp(gasname(end-1:end), 'am') % for '_am'
    gasname_short = gasname(1:end-3);
elseif length(gasname) > 3 && strcmp(gasname(end-1:end), 'pm') % for '_pm'
    gasname_short = gasname(1:end-3);
else
    gasname_short = gasname;
end

%% loop through the months and fill in the cells of arrays

filenewi = fullfile(climdir,gasname_short,newfile)
vmrzon = ncread(filenewi,gasname_short);
% whos
% vmrzoni = vmrzoni(:,6:33,:);
vmrzon(vmrzon == -999) = nan;
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
    ymax = 500;
    ymin = 60;
%     c = [0.00001 0.00002 0.00005 0.0001 0.0002 0.0005 0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 1]; % for tick marks
%     cstart = 4;
%     c = c(cstart:end);
%     c = cmax * c;
%     vmrzon(vmrzon <= 0) = nan;%c(1)*cmax;
%     vmrzon = log10(vmrzon);
%     cgrid = log10(c); % for caxis
else
    cmax = max(vmrzon(:))
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
    subplot(4,3,i), plot_on_ace_SDI_latlev_contourf( vmrzon(:,:,i)', sprintf('%0.2i',months(i)), cgrid )
%     subplot(4,3,i), plot_on_ace_SDI_latlev_pcolor( vmrzon(:,:,i)', sprintf('%0.2i',months(i)) )
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




% %% Make the plots if you want
% cmax_all = nan(1,12);
% if yplot == 1
%     %     Do for the zonal vmr
%     figi = randi(100);
%     figure(figi), set(gcf,'Position', [5,12,1096,704])
%     figure(figi), set(gcf,'Position', [71,-270,1112,954])
%     
%     suptitle(sprintf('%s zonal VMR for %i',gasname, yearin))
%     hold on
%     %     cmax1 = max(vmrzoni(:));
%     cmax1 = 12e-9;
%     for i = 1:12
%         subplot(4,3,i), plot_on_ace_SDI_latlev( vmrzoni(:,:,i), monthnames{i}(1:3) );
%         %         caxis([0,cmax1])
%         cmax_all(i) = max(max(abs(vmrzoni(:,:,i))));
%         caxis([0,cmax1]);
%         if i == 12
%             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             %             cmax1 = max(cmax_all);
%             if ~isnan(cmax1)
%                 caxis([0,cmax1]);
%             end
%             title(c,'VMR','Position', [8 -23 0]);
%         end
%     end
%     ip = 0;
%     figure(figi+1), set(gcf,'Position', [358,61,722,532])
%     suptitle(sprintf('%s zonal VMR for %i',gasname, yearin))
%     %     cmax1 = max(vmrzoni(:));
%     for i = 2:3:11
%         ip = ip+1;
%         subplot(2,2,ip), plot_on_ace_SDI_latlev( vmrzoni(:,:,i), monthnames{i}(1:3) );
%         %         caxis([0,cmax1])
%         cmax_all(ip) = max(max(abs(vmrzoni(:,:,i))));
%         if i == 11
%             c = colorbar('location','eastoutside','position',[0.923 0.309 0.016 0.35]);
%             cmax1 = max(cmax_all);
%             if ~isnan(cmax1)
%                 caxis([0,cmax1]);
%             end
%             title(c,'VMR','Position', [8 -23 0]);
%         end
%     end
    
    
    
    
    
    
    
    
    
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
% end

end

