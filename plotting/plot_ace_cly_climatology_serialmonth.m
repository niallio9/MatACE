function [ vmrzon_clo ] = plot_ace_cly_climatology_serialmonth( years_in, lat_minmax, do_plot )
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


%% Define which species you will be plotting
%  datasource = 'instrument';
% datasource = 'model';
datasource = 'both';

switch datasource
    case 'instrument'
        clo = {'ClOmls_sap'};
        hocl = {'HOClmls_sap'};
        hcl = {'HCl'};
        clono2 = {'ClONO2_sap'};
%         cly = 'ClOy';
    case 'model'
        clo = {'ClOcmam'};
        hocl = {'HOClcmam'};
        hcl = {'HClcmam'};
        clono2 = {'ClONO2cmam'};
%         cly = 'Clycmam' % uses more than the 4 gases listed above
    case 'both'
        clo = {'ClOmlspratlatnegfixampmvortex_sap','ClOcmam'};
        hocl = {'HOClmlspratlatnegfixampmvortex_sap','HOClcmam'};
        hcl = {'HCl','HClcmam'};
        clono2 = {'ClONO2','ClONO2cmam'};
%         clo = {'ClOcmam','ClOmlsprat10_sap'};
%         hocl = {'HOClcmam', 'HOClmls_sap'};
%         hcl = {'HClcmam','HCl'};
%         clono2 = {'ClONO2cmam','ClONO2'};
end

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
newclim_dir = 'C:\Users\ryann\ACE\climdata_testing\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\\MAESTRO\climdata\';

latmin = lat_minmax(1);
latmax = lat_minmax(2);
if mod(latmin,5) ~= 0 || mod(latmax,5) ~= 0
    error('you must choose latitude limits that are a multiple of 5 degrees, within [-90, 90]')
end
%         20km, 30km, 40km, 50km
%         50hPa, 10hPa, 2hPa, .5hPa
ipplot = [17, 21, 25, 29]; % the indexes of the pressure levels on which to plot
ipplot = flip(ipplot); % flip around the pressure levels for plotting
yearsin = years_in;
lyears = length(yearsin);

% vmrzon_clo = nan(48,36,12*lyears);
% vmrzon_hocl = nan(48,36,12*lyears);
% vmrzon_hcl = nan(48,36,12*lyears);
% vmrzon_clono2 = nan(48,36,12*lyears);

if nargin > 2
    yplot = do_plot;
else
    yplot = 1;
end

file_pre = 'ACEFTS_CLIM_v3_lat_'; % ACEFTS_CLIM_v3_O3_2008_12.mat
% newfile_pre = strcat('ACEMAESTRO_CLIM_v1_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_2008_12.mat

file_post = '.mat';
% newmonthnames = {'DJF', 'MAM', 'JJA', 'SON'};
months = 1:12;
sdates = nan(1,12*lyears);

%% loop through the years and months and fill in the cells of arrays
for n = 1:length(clo)
    vmrzon_clo = nan(48,36,12*lyears);
    vmrzon_hocl = nan(48,36,12*lyears);
    vmrzon_hcl = nan(48,36,12*lyears);
    vmrzon_clono2 = nan(48,36,12*lyears);
    vmrzon_clo_error = nan(48,36,12*lyears);
    vmrzon_hocl_error = nan(48,36,12*lyears);
    vmrzon_hcl_error = nan(48,36,12*lyears);
    vmrzon_clono2_error = nan(48,36,12*lyears);
    vmrzon_clo_obscount = nan(48,36,12*lyears);
    vmrzon_hocl_obscount = nan(48,36,12*lyears);
    vmrzon_hcl_obscount = nan(48,36,12*lyears);
    vmrzon_clono2_obscount = nan(48,36,12*lyears);
    ij = 0;
    n
    for j = 1:lyears
%         j
        for i = 1:12 % do all months
            ij = ij+1;
%             ij
            sdates(ij) = datenum(yearsin(j),i,1); % the 15th of each month for each year;
            
            file_cloi = strcat( newclim_dir, clo{n},'/','serial_month','/', file_pre, clo{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_hocli = strcat( newclim_dir, hocl{n},'/','serial_month','/', file_pre, hocl{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_hcli = strcat( newclim_dir, hcl{n},'/','serial_month','/', file_pre, hcl{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_clono2i = strcat( newclim_dir, clono2{n},'/','serial_month','/', file_pre, clono2{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            
            if exist(file_cloi,'file') ~= 2 || isempty(clo{n})
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', clo{n}, yearsin(j), months(i))
                vmrzon_clo(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_clo_error(:,:,ij) = nan(48,36);
                vmrzon_clo_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(file_cloi); clim = clim.climstruct; % the variable is called climstruct in the new data
                vmrzon_clo(:,:,ij) = clim.vmr_zonal;
                vmrzon_clo_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_clo_obscount(:,:,ij) = clim.obs_count;
            end
            if exist(file_hocli,'file') ~= 2 || isempty(hocl{n})
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', hocl{n}, yearsin(j), months(i))
                vmrzon_hocl(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_hocl_error(:,:,ij) = nan(48,36);
                vmrzon_hocl_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(file_hocli); clim = clim.climstruct; % the variable is called climstruct in the new data
                vmrzon_hocl(:,:,ij) = clim.vmr_zonal;
                vmrzon_hocl_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_hocl_obscount(:,:,ij) = clim.obs_count;
            end
            if exist(file_hcli,'file') ~= 2 || isempty(hcl{n})
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', hcl{n}, yearsin(j), months(i))
                vmrzon_hcl(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_hcl_error(:,:,ij) = nan(48,36);
                vmrzon_hcl_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(file_hcli); clim = clim.climstruct; % the variable is called climstruct in the new data
                vmrzon_hcl(:,:,ij) = clim.vmr_zonal;
                vmrzon_hcl_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_hcl_obscount(:,:,ij) = clim.obs_count;
            end
            if exist(file_clono2i,'file') ~= 2 || isempty(clono2{n})
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', clono2{n}, yearsin(j), months(i))
                vmrzon_clono2(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_clono2_error(:,:,ij) = nan(48,36);
                vmrzon_clono2_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(file_clono2i); clim = clim.climstruct; % the variable is called climstruct in the new data
                vmrzon_clono2(:,:,ij) = clim.vmr_zonal;
                vmrzon_clono2_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_clono2_obscount(:,:,ij) = clim.obs_count;
            end
        end
    end
    pace = clim.pressure_hPa;
    zace = clim.altitude_km_mean;
    izace_30up = find(zace>30);
    
    %% replace nans with zeros in the vmr data so that we can ignore points
    % get locations of all-nan profiles
    vmrzon_hocl_nan2zero = vmrzon_hocl;
    vmrzon_clono2_nan2zero = vmrzon_clono2;
    [Jnanprofile_hocl, Knanprofile_hocl] = find(squeeze(nansum(vmrzon_hocl_nan2zero,1) == 0));
    for i = 1:length(Jnanprofile_hocl)
       vmrzon_hocl_nan2zero(:,Jnanprofile_hocl(i),Knanprofile_hocl(i)) = 999; % change all-nan profiles to all-999 profiles 
    end
    [Jnanprofile_clono2, Knanprofile_clono2] = find(squeeze(nansum(vmrzon_clono2_nan2zero,1) == 0));
    for i = 1:length(Jnanprofile_clono2)
       vmrzon_clono2_nan2zero(:,Jnanprofile_clono2(i),Knanprofile_clono2(i)) = 999; % change all-nan profiles to all-999 profiles 
    end
    %hocl: assume missing values values are negligible
    vmrzon_hocl_nan2zero(isnan(vmrzon_hocl_nan2zero)) = 0;
    vmrzon_hocl_error_nan2zero = vmrzon_hocl_error;
    vmrzon_hocl_error_nan2zero(isnan(vmrzon_hocl_error_nan2zero)) = 0;
    %clono2: upper scaled a priori sucks so we cant use it. assume missing
    %values above 30km are negligible.
    dummy = vmrzon_clono2_nan2zero(izace_30up,:,:);
    dummy(isnan(dummy)) = 0;
    vmrzon_clono2_nan2zero(izace_30up,:,:) = dummy;
    % restore the 999 values to nans.
    vmrzon_hocl_nan2zero(vmrzon_hocl_nan2zero == 999) = nan;
    vmrzon_clono2_nan2zero(vmrzon_clono2_nan2zero == 999) = nan;
    
    %% make Cly
    %decide which types of data to use in the sum
    cly_clo = vmrzon_clo; % ignore missing data except for when there is no data in a profile at all
    cly_hocl = vmrzon_hocl; % ignore missing data except for when there is no data in a profile at all
    cly_hcl = vmrzon_hcl; % don't ignore missing data because it is always relevent
    cly_clono2 = vmrzon_clono2; % ignore missing data except for when there is no data in a profile at all
    
    cly_clo_edit = vmrzon_clo; % ignore missing data except for when there is no data in a profile at all
    cly_hocl_edit = vmrzon_hocl_nan2zero; % ignore missing data except for when there is no data in a profile at all
    cly_hcl_edit = vmrzon_hcl; % don't ignore missing data because it is always relevent
    cly_clono2_edit = vmrzon_clono2_nan2zero; % ignore missing data except for when there is no data in a profile at all
    
    cly_clo_error = vmrzon_clo_error;
    cly_hocl_error = vmrzon_hocl;
    cly_hcl_error = vmrzon_hcl_error;
    cly_clono2_error = vmrzon_clono2;
    
    vmrzon_cly = cly_clo + cly_hocl + cly_hcl + cly_clono2;
    vmrzon_cly_error = cly_clo_error + cly_hocl_error + cly_hcl_error + cly_clono2_error;
    
    vmrzon_cly_edit = cly_clo_edit + cly_hcl_edit + cly_clono2_edit;
    vmrzon_cly_error_edit = cly_clo_error + cly_hcl_error + cly_clono2_error;
    
    %% subset according to the chosen lat limits
    % average over some latitude bins if needed
    lat_bounds = clim.lat_bounds;
    if latmax < latmin
        latmin_old = latmin;
        latmax_old = latmax;
        latmin = latmax_old;
        latmax = latmin_old;
    end
    ilatmin = find(lat_bounds == latmin)
    ilatmax = find(lat_bounds == latmax) - 1
    
    vmrzon_clo = squeeze(nanmean(vmrzon_clo(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hocl = squeeze(nanmean(vmrzon_hocl(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hcl = squeeze(nanmean(vmrzon_hcl(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clono2 = squeeze(nanmean(vmrzon_clono2(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly = squeeze(nanmean(vmrzon_cly(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_edit = squeeze(nanmean(vmrzon_cly_edit(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    
    vmrzon_clo_error = squeeze(nanmean(vmrzon_clo_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hocl_error = squeeze(nanmean(vmrzon_hocl_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hcl_error = squeeze(nanmean(vmrzon_hcl_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clono2_error = squeeze(nanmean(vmrzon_clono2_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_error = squeeze(nanmean(vmrzon_cly_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_error_edit = squeeze(nanmean(vmrzon_cly_error_edit(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    
    %% Make the plots if you want
    lw = 1;
    ms = 5;
    % whos
    %% line plot
    nalt = length(ipplot);
    ytext = 4.4; % the position for the text in the plots, in ppbv
    if yplot == 1
        % Do for the zonal vmr
        %     figi = randi(100);
        figi = 19;
        %     figure(figi), set(gcf,'Position', [5,12,1096,704])
        %     figure(figi), set(gcf,'Position', [358,61,722,532])
        figure(figi), set(gcf,'Position', [97,49,852,630])
        %     figure(figi), suptitle(sprintf('Cly family zonal VMR, %i-%i', yearsin(1), yearsin(end)))
        %     return
        axx = nan(nalt,1);
        axx1 = nan(nalt,1);
        axx2 = nan(nalt,1);
        for i = 1:nalt
            figure(figi), axx(i) = subplot(nalt,1,i); hold(axx(i),'on')
            if n == 1
                    disp('plotting instrument data')
                    errorbar(axx(i),sdates, vmrzon_clo(ipplot(i),:), vmrzon_clo_error(ipplot(i),:),'go-', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_hocl(ipplot(i),:), vmrzon_hocl_error(ipplot(i),:), 'bo-', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_hcl(ipplot(i),:), vmrzon_hcl_error(ipplot(i),:), 'mo-', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_clono2(ipplot(i),:), vmrzon_clono2_error(ipplot(i),:), 'ro-', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_cly(ipplot(i),:), vmrzon_cly_error(ipplot(i),:), 'ko-', 'Linewidth', lw, 'MarkerSize', ms )
            elseif n == 2
                    disp('plotting model data')
                    errorbar(axx(i),sdates, vmrzon_clo(ipplot(i),:), vmrzon_clo_error(ipplot(i),:),'g^--', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_hocl(ipplot(i),:), vmrzon_hocl_error(ipplot(i),:), 'b^--', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_hcl(ipplot(i),:), vmrzon_hcl_error(ipplot(i),:), 'm^--', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_clono2(ipplot(i),:), vmrzon_clono2_error(ipplot(i),:), 'r^--', 'Linewidth', lw, 'MarkerSize', ms )
                    errorbar(axx(i),sdates, vmrzon_cly(ipplot(i),:), vmrzon_cly_error(ipplot(i),:), 'k^--', 'Linewidth', lw, 'MarkerSize', ms )
            end
            
            ylim([0,4.0])
            text(sdates(1), ytext, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ));
            ylabel('pressure [hPa]')
            if i == 1
                title(sprintf('Cly family zonal VMR, %i-%i%c latitude', latmin, latmax, char(176)))
                %             legend('ClO', 'HOCl', 'HCl', 'ClONO2', 'Cly', 'Location', 'eastoutside' )
            end
            if i == 4 && n == 2
%                 n
                dynamicDateTicks([],'x','mm')
%                 dynamicDateTicks
                xlabel('year')
                linkaxes(axx,'xy')
            else
%                 datetick('x','mmyy')
                set(gca,'xticklabel',{[]})
            end
        end
        
        %% contour plot
        ngas = 6;
        axx1 = nan(ngas,1);
        axx2 = nan(ngas,1);
        fs = 12;
        fs_title = 8;
        c = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5 1 2 3 4];
        cutofflow = c(1);
        cgrid = log(c);
%         vmrzon_clo(vmrzon_clo < 0) = 10^-12;
%         vmrzon_hocl(vmrzon_hocl < 0) = 10^-12;
%         vmrzon_hcl(vmrzon_hcl < 0) = 10^-12;
%         vmrzon_clono2(vmrzon_clono2 < 0) = 10^-12;
%         vmrzon_cly(vmrzon_cly < 0) = 10^-12;
%         vmrzon_cly_nohocl(vmrzon_cly_nohocl < 0) = 10^-12;
        
        vmrzon_clo(vmrzon_clo < cutofflow) = cutofflow;
        vmrzon_hocl(vmrzon_hocl < cutofflow) = cutofflow;
        vmrzon_hcl(vmrzon_hcl < cutofflow) = cutofflow;
        vmrzon_clono2(vmrzon_clono2 < cutofflow) = cutofflow;
        vmrzon_cly(vmrzon_cly < cutofflow) = cutofflow;
        vmrzon_cly_edit(vmrzon_cly_edit < cutofflow) = cutofflow;
        vmrzon_clo = log(vmrzon_clo);
        vmrzon_hocl = log(vmrzon_hocl);
        vmrzon_hcl = log(vmrzon_hcl);
        vmrzon_clono2 = log(vmrzon_clono2);
        vmrzon_cly = log(vmrzon_cly);
        vmrzon_cly_edit = log(vmrzon_cly_edit);
        
        ytickspace = [10^-4 10^-2 1 10^2];
%         cgrid = linspace(0,4,17);
% %         cgrid = [linspace(0,1,5),linspace(1.2,4,5)] % the values corresponding to the contour colours
         % Do for the zonal vmr
        %     figi = randi(100);
        figii = figi + 2;
        figure(figii), set(gcf,'Position', [97,49,852,630])

            ylim2 = 300;
            ylim1 = 10^-1;
            if n == 1
                figure(figii),
                fignames = {'ClO\_meas', 'HOCl\_meas', 'HCl\_meas', 'ClONO2\_meas', 'Cly\_meas', 'Cly\_meas (no HOCl, ClONO2 above 40 km)'};
%                                 fignames = {'ClO\_cmam', 'HOCl\_cmam', 'HCl\_cmam', 'ClONO2\_cmam', 'Cly\_cmam', 'Cly\_cmam (no HOCl, ClONO2 above 40 km)'};
                disp('plotting instrument data')
                i = 1; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i), sdates, pace, vmrzon_clo, cgrid); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 2; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i), sdates, pace, vmrzon_hocl, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 3; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i),sdates, pace, vmrzon_hcl, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]);  % colormap(parula(length(cgrid) - 1));
                i = 4; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i),sdates, pace, vmrzon_clono2, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 5; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i),sdates, pace, vmrzon_cly, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 6; axx1(i) = subplot(ngas,1,i); hold(axx1(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx1(i),sdates, pace, vmrzon_cly_edit, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); % colormap(parula(length(cgrid) - 1));
%                 colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','southoutside')
                xlabel('year')
                ylabel('pressure [hPa]')
                ylim([ylim1 ylim2])
                linkaxes(axx1,'xy')
                dynamicDateTicks([],'x','mmyy')
            elseif n == 2
                figure(figii + 1), set(gcf,'Position', [97,49,852,630])
                fignames = {'ClO\_cmam', 'HOCl\_cmam', 'HCl\_cmam', 'ClONO2\_cmam', 'Cly\_cmam', 'Cly\_cmam (no HOCl, ClONO2 above 40 km)'};
                disp('plotting model data')
                i = 1; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i), sdates, pace, vmrzon_clo, cgrid); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 2; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i), sdates, pace, vmrzon_hocl, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 3; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i),sdates, pace, vmrzon_hcl, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 4; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i),sdates, pace, vmrzon_clono2, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 5; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i),sdates, pace, vmrzon_cly, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); set(gca,'xticklabel',{[]}); ylim([ylim1 ylim2]); % colormap(parula(length(cgrid) - 1));
                i = 6; axx2(i) = subplot(ngas,1,i); hold(axx2(i),'on'), box on, title(fignames{i},'FontSize',fs_title)
                contourf(axx2(i),sdates, pace, vmrzon_cly_edit, cgrid ); caxis([cgrid(1),cgrid(end)])
                set(gca, 'Ydir','reverse', 'YScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
                yticks(ytickspace); % colormap(parula(length(cgrid) - 1));
                colorbar('YTick', cgrid, 'Yticklabel', c, 'Location','eastoutside')
                xlabel('year')
                ylabel('pressure [hPa]')
                ylim([ylim1 ylim2])
                linkaxes(axx2,'xy')
                dynamicDateTicks([],'x','mmyy')
            end
            
%             ylim([0,4.0])
%             text(sdates(1), ytext, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ));
            
%             if i == 1
%                 title(sprintf('Cly family zonal VMR, %i-%i%c latitude', latmin, latmax, char(176)))
%                 %             legend('ClO', 'HOCl', 'HCl', 'ClONO2', 'Cly', 'Location', 'eastoutside' )
%             end
%             if i == 6 && n == 1
% %                 dynamicDateTicks([],'x','mmyy')
% %                 datetick
%                 xlabel('year')
%                 ylabel('pressure [hPa]')
%                 linkaxes(axx1,'xy')
%             else
%                 datetick('x','mmyy')
%                 set(gca,'xticklabel',{[]})
%             end
    end
end
%
end

