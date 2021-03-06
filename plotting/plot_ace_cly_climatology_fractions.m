function [ vmrzon_cly_error, vmrzon_hcl_error ] = plot_ace_cly_climatology_fractions( lat_minmax, do_plot )
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define which species you will be plotting
%  datasource = 'instrument';
% datasource = 'model';
datasource = 'both';
file_pre = 'ACEFTS_CLIM_v3_lat_'; % ACEFTS_CLIM_v3_lat_O3_DJF.mat

file_post = '_JJA.mat';
% monthnames = {'DJF', 'MAM', 'JJA', 'SON'};
cgrey = 0.5; % for plotting in grey
newcolours = get(groot,'DefaultAxesColorOrder');
blue = newcolours(1,:);
red = newcolours(7, :);
green = newcolours(5, :);
purple = newcolours(4, :);

switch datasource
    case 'instrument'
%         clo = {'ClOmlspratlatnegfixampmvortex'};
        clo = {'ClOmlsonly'};
        hocl = {'HOClmls_sap'};
        hcl = {'HClv4'};
        clono2 = {'ClONO2v4'};
        %         cly = 'ClOy';
    case 'model'
        clo = {'ClOcmam'};
        hocl = {'HOClcmam'};
        hcl = {'HClcmam'};
        clono2 = {'ClONO2cmam'};
        %         cly = 'Clycmam' % uses more than the 4 gases listed above
    case 'both'
        clo = {'ClOv4withmlssmiles_sap','ClOcmam'};
        hocl = {'HOClmlsfrac10lim4ppb_sap','HOClcmam'};
        hcl = {'HClv4','HClcmam'};
        clono2 = {'ClONO2v4','ClONO2cmam'};
end
figi = randi(100);
%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
clim_dir = 'C:\Users\ryann\ACE\climdata_testing\time_matched_climatology\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\MAESTRO\climdata\';

latmin = lat_minmax(1);
latmax = lat_minmax(2);
if mod(latmin,5) ~= 0 || mod(latmax,5) ~= 0
    error('you must choose latitude limits that are a multiple of 5 degrees, within [-90, 90]')
end
% %         20km, 30km, 40km, 50km
% %         50hPa, 10hPa, 2hPa, .5hPa
% ipplot = [17, 21, 25, 29]; % the indexes of the pressure levels on which to plot
% ipplot = flip(ipplot); % flip around the pressure levels for plotting
% yearsin = years_in;
% lyears = length(yearsin);

% vmrzon_clo = nan(48,36,12*lyears);
% vmrzon_hocl = nan(48,36,12*lyears);
% vmrzon_hcl = nan(48,36,12*lyears);
% vmrzon_clono2 = nan(48,36,12*lyears);

if nargin > 1
    yplot = do_plot;
else
    yplot = 1;
end

%% read in the data
for n = 1:length(clo)
    file_clo = strcat( clim_dir, clo{n},'/', file_pre, clo{n}, file_post);
    file_hocl = strcat( clim_dir, hocl{n},'/', file_pre, hocl{n}, file_post);
    file_hcl = strcat( clim_dir, hcl{n},'/', file_pre, hcl{n}, file_post);
    file_clono2 = strcat( clim_dir, clono2{n},'/', file_pre, clono2{n}, file_post);
    % for total Cly from CMAM
    file_clytot = strcat( clim_dir, 'clytotcmam','/', file_pre, 'clytotcmam', file_post);
    
    if exist(file_clo,'file') ~= 2 || isempty(clo{n})
        fprintf('There is no file for %s%s. Moving on...\n', clo{n}, file_post(1:end-4))
        vmrzon_clo = nan(48,36); % this is hardcoded here, for now.
        vmrzon_clo_error = nan(48,36);
    else
        clim = load(file_clo); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon_clo = clim.vmr_zonal;
%         vmrzon_clo_error = sqrt(clim.vmr_zonal_var);
        vmrzon_clo_error = clim.vmr_zonal_standard_error;
    end
    if exist(file_hocl,'file') ~= 2 || isempty(hocl{n})
        fprintf('There is no file for %s%s. Moving on...\n', hocl{n}, file_post(1:end-4))
        vmrzon_hocl = nan(48,36); % this is hardcoded here, for now.
        vmrzon_hocl_error = nan(48,36);
    else
        clim = load(file_hocl); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon_hocl = clim.vmr_zonal;
%         vmrzon_hocl_error = sqrt(clim.vmr_zonal_var);
        vmrzon_hocl_error = clim.vmr_zonal_standard_error;
    end
    if exist(file_hcl,'file') ~= 2 || isempty(hcl{n})
        fprintf('There is no file for %s%s. Moving on...\n', hcl{n}, file_post(1:end-4))
        vmrzon_hcl = nan(48,36); % this is hardcoded here, for now.
        vmrzon_hcl_error = nan(48,36);
    else
        clim = load(file_hcl); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon_hcl = clim.vmr_zonal;
%         vmrzon_hcl_error = sqrt(clim.vmr_zonal_var);
        vmrzon_hcl_error = clim.vmr_zonal_standard_error;
    end
    if exist(file_clono2,'file') ~= 2 || isempty(clono2{n})
        fprintf('There is no file for %s%s. Moving on...\n', clono2{n}, file_post(1:end-4))
        vmrzon_clono2 = nan(48,36); % this is hardcoded here, for now.
        vmrzon_clono2_error = nan(48,36);
    else
        clim = load(file_clono2); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon_clono2 = clim.vmr_zonal;
%         vmrzon_clono2_error = sqrt(clim.vmr_zonal_var);
                vmrzon_clono2_error = clim.vmr_zonal_standard_error;
    end
    if exist(file_clytot,'file') ~= 2
        fprintf('There is no file for %s%s. Moving on...\n', 'clytotcmam', file_post(1:end-4))
        vmrzon_clytotcmam = nan(48,36); % this is hardcoded here, for now.
        vmrzon_clytotcmam_error = nan(48,36);
    else
        clim = load(file_clytot); clim = clim.climstruct; % the variable is called climstruct in the new data
        vmrzon_clytotcmam = clim.vmr_zonal;
%         vmrzon_clytotcmam_error = sqrt(clim.vmr_zonal_var);
                vmrzon_clytotcmam_error = clim.vmr_zonal_standard_error;
    end
    
%     pace = clim.pressure_hPa;
    zace = clim.altitude_km_mean;
    izace_30up = find(zace>30);
    
    %% replace nans with zeros in the vmr data so that we can ignore points
    % get locations of all-nan profiles
    vmrzon_hocl_nan2zero = vmrzon_hocl;
    vmrzon_clono2_nan2zero = vmrzon_clono2;
    Jnanprofile_hocl = find(nansum(vmrzon_hocl_nan2zero,1) == 0);
    for i = 1:length(Jnanprofile_hocl)
       vmrzon_hocl_nan2zero(:,Jnanprofile_hocl(i)) = 999; % change all-nan profiles to all-999 profiles 
    end
    Jnanprofile_clono2 = find(nansum(vmrzon_clono2_nan2zero,1) == 0);
    for i = 1:length(Jnanprofile_clono2)
       vmrzon_clono2_nan2zero(Jnanprofile_clono2(i)) = 999; % change all-nan profiles to all-999 profiles 
    end
    %hocl: assume missing values values are negligible
    vmrzon_hocl_nan2zero(isnan(vmrzon_hocl_nan2zero)) = 0;
    vmrzon_hocl_error_nan2zero = vmrzon_hocl_error;
    vmrzon_hocl_error_nan2zero(isnan(vmrzon_hocl_error_nan2zero)) = 0;
    %clono2: upper scaled a priori sucks so we cant use it. assume missing
    %values above 30km are negligible.
    dummy = vmrzon_clono2_nan2zero(izace_30up,:);
    dummy(isnan(dummy)) = 0;
    vmrzon_clono2_nan2zero(izace_30up,:) = dummy;
    % restore the 999 values to nans.
    vmrzon_hocl_nan2zero(vmrzon_hocl_nan2zero == 999) = nan;
    vmrzon_clono2_nan2zero(vmrzon_clono2_nan2zero == 999) = nan;
    vmrzon_clono2_error_nan2zero = vmrzon_clono2_error;
    vmrzon_clono2_error_nan2zero(isnan(vmrzon_clono2_error_nan2zero)) = 0;
    
    %% make Cly
    %decide which types of data to use in the sum
    cly_clo = vmrzon_clo; % ignore missing data except for when there is no data in a profile at all
    cly_hocl = vmrzon_hocl_nan2zero; % ignore missing data except for when there is no data in a profile at all
    cly_hcl = vmrzon_hcl; % don't ignore missing data because it is always relevent
    cly_clono2 = vmrzon_clono2_nan2zero; % ignore missing data except for when there is no data in a profile at all
    
    
    cly_clo_error = vmrzon_clo_error;
    cly_hocl_error = vmrzon_hocl_error_nan2zero;
    cly_hcl_error = vmrzon_hcl_error;
    cly_clono2_error = vmrzon_clono2_error_nan2zero;
    
    vmrzon_cly = cly_clo + cly_hocl + cly_hcl + cly_clono2;
    vmrzon_cly_error = sqrt(cly_clo_error.^2 + cly_hocl_error.^2 + cly_hcl_error.^2 + cly_clono2_error.^2);
    
    %%
%     cutofflow = 0.01e-9;
%     vmrzon_hocl(vmrzon_hocl < cutofflow) = cutofflow;
    vmrzon_hocl = vmrzon_hocl_nan2zero;
    %%
    
    
%     vmrzon_clo_nan2zero = vmrzon_clo;
%     vmrzon_clo_nan2zero(isnan(vmrzon_clo_nan2zero)) = 0;
%     vmrzon_hocl_nan2zero = vmrzon_hocl;
%     vmrzon_hocl_nan2zero(isnan(vmrzon_hocl_nan2zero)) = 0;
% %     vmrzon_hcl_nan2zero = vmrzon_hcl;
% %     vmrzon_hcl_nan2zero(isnan(vmrzon_hcl_nan2zero)) = 0;
%     vmrzon_clono2_nan2zero = vmrzon_clono2;
%     vmrzon_clono2_nan2zero(isnan(vmrzon_clono2_nan2zero)) = 0;
%     
%     vmrzon_clo_error_nan2zero = vmrzon_clo_error;
%     vmrzon_clo_error_nan2zero(isnan(vmrzon_clo_error_nan2zero)) = 0;
%     vmrzon_hocl_error_nan2zero = vmrzon_hocl_error;
%     vmrzon_hocl_error_nan2zero(isnan(vmrzon_hocl_error_nan2zero)) = 0;
% %     vmrzon_hcl_error_nan2zero = vmrzon_hcl_error;
% %     vmrzon_hcl_error_nan2zero(isnan(vmrzon_hcl_error_nan2zero)) = 0;
%     vmrzon_clono2_error_nan2zero = vmrzon_clono2_error;
%     vmrzon_clono2_error_nan2zero(isnan(vmrzon_clono2_error_nan2zero)) = 0;
%     
% %     vmrzon_cly = vmrzon_clo + vmrzon_hocl + vmrzon_hcl + vmrzon_clono2;
% %     vmrzon_cly_error = sqrt(vmrzon_clo_error.^2 + vmrzon_hocl_error.^2 + vmrzon_hcl_error.^2 + vmrzon_clono2_error.^2);
%     vmrzon_cly = vmrzon_clo_nan2zero + vmrzon_hocl_nan2zero + vmrzon_hcl + vmrzon_clono2_nan2zero;
%     vmrzon_cly_error = sqrt(vmrzon_clo_error_nan2zero.^2 + vmrzon_hocl_error_nan2zero.^2 + vmrzon_hcl_error.^2 + vmrzon_clono2_error_nan2zero.^2);
% %     vmrzon_cly = vmrzon_clo + vmrzon_hcl + vmrzon_clono2;
%     vmrzon_cly(vmrzon_cly == 0) = nan;
%     vmrzon_cly_error(vmrzon_cly_error == 0) = nan; % need to replace the nans for averaging over latitude later

    %% subset according to the chosen lat limits
    % average over some latitude bins if needed
    lat_bounds = clim.lat_bounds;
    if latmax < latmin
        latmin_old = latmin;
        latmax_old = latmax;
        latmin = latmax_old;
        latmax = latmin_old;
    end
    ilatmin = find(lat_bounds == latmin);
    ilatmax = find(lat_bounds == latmax) - 1;
    
    vmrzon_clo = nanmean(vmrzon_clo(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_hocl = nanmean(vmrzon_hocl(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_hcl = nanmean(vmrzon_hcl(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_clono2 = nanmean(vmrzon_clono2(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_cly = nanmean(vmrzon_cly(:,ilatmin:ilatmax), 2) * 1e9;
    if n == 2
        vmrzon_clytotcmam = nanmean(vmrzon_clytotcmam(:,ilatmin:ilatmax), 2) * 1e9;
    end
    
    vmrzon_clo_error = nanmean(vmrzon_clo_error(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_hocl_error = nanmean(vmrzon_hocl_error(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_hcl_error = nanmean(vmrzon_hcl_error(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_clono2_error = nanmean(vmrzon_clono2_error(:,ilatmin:ilatmax), 2) * 1e9;
    vmrzon_cly_error = nanmean(vmrzon_cly_error(:,ilatmin:ilatmax), 2) * 1e9;
    if n == 2
        vmrzon_clytotcmam_error = nanmean(vmrzon_clytotcmam_error(:,ilatmin:ilatmax), 2) * 1e9;
    end
    
    %% make the fractional data
    clo_frac = vmrzon_clo./vmrzon_cly;
    hocl_frac = vmrzon_hocl./vmrzon_cly;
    hcl_frac = vmrzon_hcl./vmrzon_cly;
    clono2_frac = vmrzon_clono2./vmrzon_cly;
    
    clo_frac_error = vmrzon_clo_error./vmrzon_cly;
    hocl_frac_error = vmrzon_hocl_error./vmrzon_cly;
    hcl_frac_error = vmrzon_hcl_error./vmrzon_cly;
    clono2_frac_error = vmrzon_clono2_error./vmrzon_cly;
    
    
    %% Make the plots if you want
    pace = clim.pressure_hPa;
    zace = clim.altitude_km_mean;
    lw = 2;
    ms = 5;
    fs = 16;
%     ytickspace = [10^-4 10^-2 1 10^2];
%     ylim2 = 300;
%     ylim1 = 2;
    plim2 = pace(8);
    plim1 = pace(31);
    zlim1 = zace(8);
    zlim2 = zace(31);
    % whos
    %% line plot of fractions and of absolute vmrs
    if yplot == 1
%         figpos = [-1227 367 655 544];
%         figpos = [-1227 438 565 473];
        figpos = [3         380        1855         465]; % office monitor
%         figpos = [268   385   565   473];


%         figpos = [97,49,852,630];
        % Do for the zonal vmr
        %     figi = randi(100);
%         figi = 60;
        %     figure(figi), set(gcf,'Position', [5,12,1096,704])
        %     figure(figi), set(gcf,'Position', [358,61,722,532])
        figure(figi), set(gcf,'Position', figpos), box on
% %         figii = figi+1;
% %         figure(figii), set(gcf,'Position', figpos), box on
% %         figiii = figi+2;
% %         figure(figiii), set(gcf,'Position', figpos), box on
        %     figure(figi), suptitle(sprintf('Cly family zonal VMR, %i-%i', yearsin(1), yearsin(end)))
        %     return
            figure(figi), subplot(1, 3, 3), box on, hold on;
            if n == 1
                disp('plotting instrument data')
                yyaxis left
                errorbar(clo_frac, pace, clo_frac_error, 'horizontal', '^-', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(hocl_frac, pace, hocl_frac_error, 'horizontal', 'x-', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(hcl_frac, pace, hcl_frac_error, 'horizontal', 'd-', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(clono2_frac, pace, clono2_frac_error, 'horizontal', 's-', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
            elseif n == 2
                disp('plotting model data')
                yyaxis left
%                 errorbar(clo_frac, pace, clo_frac_error, 'horizontal', 'g^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(hocl_frac, pace, hocl_frac_error, 'horizontal', 'bx--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(hcl_frac, pace, hcl_frac_error, 'horizontal', 'md--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(clono2_frac, pace, clono2_frac_error, 'horizontal', 'rs--', 'Linewidth', lw, 'MarkerSize', ms )
                
                plot(clo_frac, pace, 'g^--', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                plot(hocl_frac, pace, 'x--', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                plot(hcl_frac, pace, 'd--', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                plot(clono2_frac, pace, 's--', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
                
                yyaxis right
                dummy = errorbar(clono2_frac, zace, clono2_frac_error, 'horizontal', 'r--', 'Linewidth', lw, 'MarkerSize', ms );
                delete(dummy)
            end
%             title(sprintf('Cly family zonal VMR fraction, %i-%i%c latitude', latmin, latmax, char(176)))
            yyaxis left
            yticklabels([])
            set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
            ylim([plim1 plim2]);
            yyaxis right
            ylim([zlim1 zlim2]);
            ylabel('altitude [km]')
            xlabel('VMR fraction of Cly [ppbv/ppbv]')
            xlim([-0.1, 1.1])
            
            %absolute vmrs
            %positive error bar
            vmrzon_clo_error_pos = vmrzon_clo_error;
            vmrzon_hocl_error_pos = vmrzon_hocl_error;
            vmrzon_hcl_error_pos = vmrzon_hcl_error;
            vmrzon_clono2_error_pos = vmrzon_clono2_error;
            vmrzon_cly_error_pos = vmrzon_cly_error;
            if n == 2
                vmrzon_clytotcmam_error_pos = vmrzon_clytotcmam_error;
            end
            % negative error bar
            vmrzon_clo_error_neg = vmrzon_clo_error;
            vmrzon_hocl_error_neg = vmrzon_hocl_error;
            vmrzon_hcl_error_neg = vmrzon_hcl_error;
            vmrzon_clono2_error_neg = vmrzon_clono2_error;
            vmrzon_cly_error_neg = vmrzon_cly_error;
            if n == 2
                vmrzon_clytotcmam_error_neg = vmrzon_clytotcmam_error;
            end
            
% %             % with a logarithmic scale
% %             c = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 4];
% %             cgrid = log(c);
% %             vmrzon_clo(vmrzon_clo < 0) = 10^-12;
% %             vmrzon_hocl(vmrzon_hocl < 0) = 10^-12;
% %             vmrzon_hcl(vmrzon_hcl < 0) = 10^-12;
% %             vmrzon_clono2(vmrzon_clono2 < 0) = 10^-12;
% %             vmrzon_cly(vmrzon_cly < 0) = 10^-12;
% %             vmrzon_clo = log(vmrzon_clo);
% %             vmrzon_hocl = log(vmrzon_hocl);
% %             vmrzon_hcl = log(vmrzon_hcl);
% %             vmrzon_clono2 = log(vmrzon_clono2);
% %             vmrzon_cly = log(vmrzon_cly);
% %             
% %             %positive error bar
% %             vmrzon_clo_error_pos = log(1 + vmrzon_clo_error./vmrzon_clo); % errors are transformed this way because they must be added to the vmr data in log space
% %             vmrzon_hocl_error_pos = log(1 + vmrzon_hocl_error./vmrzon_hocl);
% %             vmrzon_hcl_error_pos = log(1 + vmrzon_hcl_error./vmrzon_hcl);
% %             vmrzon_clono2_error_pos = log(1 + vmrzon_clono2_error)./vmrzon_clono2;
% %             vmrzon_cly_error_pos = log(1 + vmrzon_cly_error./vmrzon_cly);
% %             % negative error bar
% %             vmrzon_clo_error_neg = log(vmrzon_clo./(vmrzon_clo - vmrzon_clo_error));
% %             vmrzon_hocl_error_neg = log(vmrzon_hocl./(vmrzon_hocl - vmrzon_hocl_error));
% %             vmrzon_hcl_error_neg = log(vmrzon_hcl./(vmrzon_hcl - vmrzon_hcl_error));
% %             vmrzon_clono2_error_neg = log(vmrzon_clono2./(vmrzon_clono2 - vmrzon_clono2_error));
% %             vmrzon_cly_error_neg = log(vmrzon_cly./(vmrzon_cly - vmrzon_cly_error));
            
%             vmrzon_clo_error = log(vmrzon_clo_error);
%             vmrzon_hocl_error = log(vmrzon_hocl_error);
%             vmrzon_hcl_error = log(vmrzon_hcl_error);
%             vmrzon_clono2_error = log(vmrzon_clono2_error);
%             vmrzon_cly_error = log(vmrzon_cly_error);
            
%             figure(figii), hold on;
            figure(figi), subplot(1, 3, 1), box on, hold on;
            if n == 1
                disp('plotting instrument data')
                yyaxis left
                errorbar(vmrzon_clo, pace, vmrzon_clo_error_neg, vmrzon_clo_error_pos, 'horizontal', '^-', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_hocl, pace, vmrzon_hocl_error_neg, vmrzon_hocl_error_pos, 'horizontal', 'x-', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_hcl, pace, vmrzon_hcl_error_neg, vmrzon_hcl_error_pos, 'horizontal', 'd-', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_clono2, pace, vmrzon_clono2_error_neg, vmrzon_clono2_error_pos, 'horizontal', 's-', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar( vmrzon_cly, pace, vmrzon_cly_error_pos,vmrzon_cly_error_neg, 'horizontal', 'ko-', 'Linewidth', lw, 'MarkerSize', ms )
            elseif n == 2
                disp('plotting model data')
                yyaxis left
%                 errorbar(vmrzon_clo, pace, vmrzon_clo_error_neg, vmrzon_clo_error_pos, 'horizontal', 'g^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_hocl, pace, vmrzon_hocl_error_neg, vmrzon_hocl_error_pos, 'horizontal', 'b^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_hcl, pace, vmrzon_hcl_error_neg, vmrzon_hcl_error_pos, 'horizontal', 'm^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_clono2, pace, vmrzon_clono2_error_neg, vmrzon_clono2_error_pos, 'horizontal', 'r^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_cly, pace, vmrzon_cly_error_neg, vmrzon_cly_error_pos, 'horizontal', 'k^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_clytotcmam, pace, vmrzon_clytotcmam_error_neg, vmrzon_clytotcmam_error_pos, 'horizontal',  'color', [0 0 0]+cgrey , 'Linewidth', lw, 'MarkerSize', ms )
                
                
                plot(vmrzon_clo, pace, '^--', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_hocl, pace, 'x--', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_hcl, pace, 'd--', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_clono2, pace, 's--', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_cly, pace, 'ko--', 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_clytotcmam, pace, '--', 'color', [0 0 0]+cgrey , 'Linewidth', lw, 'MarkerSize', ms )
                
%                 yyaxis right
%                 dummy = errorbar(vmrzon_cly, zace, vmrzon_cly_error, 'horizontal', 'k^--', 'Linewidth', lw, 'MarkerSize', ms );
%                 delete(dummy)
            end
%             title(sprintf('Cly family zonal VMR, %i-%i%c latitude', latmin, latmax, char(176)))
            yyaxis left
% %             set(gca, 'Ydir','reverse', 'YScale', 'log', 'Xtick', cgrid, 'XTicklabel', c, 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
            set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);

            ylabel('pressure [hPa]')
            ylim([plim1 plim2]);
            yyaxis right
            yticklabels([])
%             ylim([zlim1 zlim2]);
%             ylabel('altitude [km]')
            xlabel('VMR [ppbv]')
            xlim([0.001 4])
            
            %% the errors have to be fixed for the logarithmic plot
            %logarithmic x scale
            %positive error bar
            vmrzon_clo_error_pos = vmrzon_clo_error;
            vmrzon_hocl_error_pos = vmrzon_hocl_error;
            vmrzon_hcl_error_pos = vmrzon_hcl_error;
            vmrzon_clono2_error_pos = vmrzon_clono2_error;
            vmrzon_cly_error_pos = vmrzon_cly_error;
            if n == 2
                vmrzon_clytotcmam_error_pos = vmrzon_clytotcmam_error;
            end
            % negative error bar
            addon = -1e-5;
            if vmrzon_clo - vmrzon_clo_error < 0
                vmrzon_clo_error_neg = vmrzon_clo + addon;
            else
                vmrzon_clo_error_neg = vmrzon_clo_error;
            end
            if vmrzon_hocl - vmrzon_hocl_error < 0
                vmrzon_hocl_error_neg = vmrzon_hocl + addon;
            else
                vmrzon_hocl_error_neg = vmrzon_hocl_error;
            end
            if vmrzon_hcl - vmrzon_hcl_error < 0
                vmrzon_hcl_error_neg = vmrzon_hcl + addon;
            else
                vmrzon_hcl_error_neg = vmrzon_hcl_error;
            end
            if vmrzon_clono2 - vmrzon_clono2_error < 0
                vmrzon_clono2_error_neg = vmrzon_clono2 + addon;
            else
                vmrzon_clono2_error_neg = vmrzon_clono2_error;
            end
            if vmrzon_cly - vmrzon_cly_error < 0
                vmrzon_cly_error_neg = vmrzon_cly + addon;
            else
                vmrzon_cly_error_neg = vmrzon_cly_error;
            end
            if n == 2
                if vmrzon_clytotcmam - vmrzon_clytotcmam_error < 0
                    vmrzon_clytotcmam_error_neg = vmrzon_clytotcmam + addon;
                else
                    vmrzon_clytotcmam_error_neg = vmrzon_clytotcmam_error;
                end
            end
%             figure(figiii), hold on;
            figure(figi), subplot(1, 3, 2), box on, hold on;
            if n == 1
                disp('plotting instrument data')
                yyaxis left
                errorbar(vmrzon_clo, pace, vmrzon_clo_error_neg, vmrzon_clo_error_pos, 'horizontal', '^-', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_hocl, pace, vmrzon_hocl_error_neg, vmrzon_hocl_error_pos, 'horizontal', 'x-', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_hcl, pace, vmrzon_hcl_error_neg, vmrzon_hcl_error_pos, 'horizontal', 'd-', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar(vmrzon_clono2, pace, vmrzon_clono2_error_neg, vmrzon_clono2_error_pos, 'horizontal', 's-', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
                errorbar( vmrzon_cly, pace, vmrzon_cly_error_pos,vmrzon_cly_error_neg, 'horizontal', 'ko-', 'Linewidth', lw, 'MarkerSize', ms )
            elseif n == 2
                disp('plotting model data')
                yyaxis left
%                 errorbar(vmrzon_clo, pace, vmrzon_clo_error_neg, vmrzon_clo_error_pos, 'horizontal', 'g^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_hocl, pace, vmrzon_hocl_error_neg, vmrzon_hocl_error_pos, 'horizontal', 'b^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_hcl, pace, vmrzon_hcl_error_neg, vmrzon_hcl_error_pos, 'horizontal', 'm^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_clono2, pace, vmrzon_clono2_error_neg, vmrzon_clono2_error_pos, 'horizontal', 'r^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_cly, pace, vmrzon_cly_error_neg, vmrzon_cly_error_pos, 'horizontal', 'k^--', 'Linewidth', lw, 'MarkerSize', ms )
%                 errorbar(vmrzon_clytotcmam, pace, vmrzon_clytotcmam_error_neg, vmrzon_clytotcmam_error_pos, 'horizontal', 'color', [0 0 0]+cgrey , 'Linewidth', lw, 'MarkerSize', ms )
                
                plot(vmrzon_clo, pace, '^--', 'color', green, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_hocl, pace, 'x--', 'color', blue, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_hcl, pace, 'd--', 'color', purple, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_clono2, pace, 's--', 'color', red, 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_cly, pace, 'ko--', 'Linewidth', lw, 'MarkerSize', ms )
                plot(vmrzon_clytotcmam, pace, '--', 'color', [0 0 0]+cgrey , 'Linewidth', lw, 'MarkerSize', ms )
                
                
%                 yyaxis right
%                 dummy = errorbar(vmrzon_cly, zace, vmrzon_cly_error, 'horizontal', 'k^--', 'Linewidth', lw, 'MarkerSize', ms );
%                 delete(dummy)
            end
            title(sprintf('Cly family zonal VMR, %i-%i%c latitude', latmin, latmax, char(176)))
            yyaxis left
            yticklabels([])
% %             set(gca, 'Ydir','reverse', 'YScale', 'log', 'Xtick', cgrid, 'XTicklabel', c, 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
            set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'log', 'YMinorTick','on','YMinorGrid','OFF', 'XGrid', 'ON', 'FontSize',fs);
            
%             ylabel('pressure [hPa]')
            ylim([plim1 plim2]);
            yyaxis right
            yticklabels([])
%             ylim([zlim1 zlim2]);
%             ylabel('altitude [km]')
            xlabel('VMR [ppbv]')
            xlim([0.001 4])
            set(gca,'Xtick',[1e-3,1e-2,1e-1,1e0],'XTickLabel',[0.001, 0.01, 0.1, 1])
            
            
            
            
% %             xlim([cgrid(1) cgrid(end)])

%             set(ax1,'units','normalized','position',[0.1 0.1 0.4 0.8]);
%             set(ax2,'units','normalized','position',[0.5 0.1 0.4 0.8]);
%             set(ax1,'xscale','log','xlim',[0 0.5]);
%             set(ax2,'xscale','linear','xlim',[0.5 4]);
%             set([ax1 ax2],'ylim',[plim1 plim2],'box','off','Ydir','reverse', 'YScale', 'log');
%             set([ax1 ax2],'box','off');
            
    end
    
end
%
end

