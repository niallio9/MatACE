function [ trend, test2, fit_stats ] = plot_ace_cly_trend( years_in, lat_minmax, do_plot )
%A funcion to compare the climatology made by Jaho, with the current
%version of the climatology. The assumption is that the two versions are
%made on the same latitude and altitude grid, and ar for the same gas.

% *INPUT*
%           gasname_in: STRING - the name of the gas for which you want to
%           compare the climatologies.
%
%           years_in: STRING - the years over which you want to plot.
%
%           lat_minmax: STRING - the minimum and maximum of the latitude
%           range that you want to plot: [min, max].
%
% *OUTPUT*
%           fit_stats: STRUCTURE - the properties of the fit to the data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define which species you will be plotting
%  datasource = 'instrument';
% datasource = 'model';
datasource = 'both';
monthnames = {'DJF', 'MAM', 'JJA', 'SON'};
cgrey = 0.5; % for plotting in grey
newcolours = get(groot,'DefaultAxesColorOrder');
blue = newcolours(1,:);
red = newcolours(7, :);
green = newcolours(5, :);
purple = newcolours(4, :);

yearsin = years_in;
lyears = length(yearsin);
months = 1:12;
sdates = nan(1,12*lyears);
latmin = lat_minmax(1);
latmax = lat_minmax(2);
if mod(latmin,5) ~= 0 || mod(latmax,5) ~= 0
    error('you must choose latitude limits that are a multiple of 5 degrees, within [-90, 90]')
end
%         20km, 30km, 40km, 50km
%         50hPa, 10hPa, 2hPa, .5hPa
ipplot = [17, 21, 25, 29]; % the indexes of the pressure levels on which to plot
ipplot = [21]; % the indexes of the pressure levels on which to plot
% ipplot = 17:29;
ipplot = flip(ipplot); % flip around the pressure levels for plotting
nalt = length(ipplot);
if nargin > 2
    yplot = do_plot;
else
    yplot = 0;
end
file_pre = 'ACEFTS_CLIM_v3_lat_'; % ACEFTS_CLIM_v3_O3_2008_12.mat
% newfile_pre = strcat('ACEMAESTRO_CLIM_v1_lat_',gasname,'_'); % ACEFTS_CLIM_v3_O3_2008_12.mat
file_post = '.mat';

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
        %         clo = {'ClOmlsfrac10lim4ppb_sap','ClOcmam'};
        hocl = {'HOClmlsfrac10lim4ppb_sap','HOClcmam'};
        hcl = {'HClv4','HClcmam'};
        clono2 = {'ClONO2v4','ClONO2cmam'};
        n2o = {'N2O', 'N2Ocmam'};
end
figi = randi(100);

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
clim_dir = 'C:\Users\ryann\ACE\climdata_testing\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\MAESTRO\climdata\';

%% read in the data
for n = 1:length(clo)
    vmrzon_clo = nan(48,36,12*lyears);
    vmrzon_hocl = nan(48,36,12*lyears);
    vmrzon_hcl = nan(48,36,12*lyears);
    vmrzon_clono2 = nan(48,36,12*lyears);
    vmrzon_n2o = nan(48,36,12*lyears);
    vmrzon_clytot = nan(48,36,12*lyears);
    %
    vmrzon_clo_error = nan(48,36,12*lyears);
    vmrzon_hocl_error = nan(48,36,12*lyears);
    vmrzon_hcl_error = nan(48,36,12*lyears);
    vmrzon_clono2_error = nan(48,36,12*lyears);
    vmrzon_n2o_error = nan(48,36,12*lyears);
    vmrzon_clytot_error = nan(48,36,12*lyears);
    %
    vmrzon_clo_obscount = nan(48,36,12*lyears);
    vmrzon_hocl_obscount = nan(48,36,12*lyears);
    vmrzon_hcl_obscount = nan(48,36,12*lyears);
    vmrzon_clono2_obscount = nan(48,36,12*lyears);
    vmrzon_n2o_obscount = nan(48,36,12*lyears);
    vmrzon_clytot_obscount = nan(48,36,12*lyears);
    ij = 0;
    n
    for j = 1:lyears
        
        for i = 1:12 % do all months
            ij = ij+1;
            %             ij
            sdates(ij) = datenum(yearsin(j),i,15); % the 15th of each month for each year;
            
            file_cloi = strcat( clim_dir, clo{n},'/','serial_month','/', file_pre, clo{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_hocli = strcat( clim_dir, hocl{n},'/','serial_month','/', file_pre, hocl{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_hcli = strcat( clim_dir, hcl{n},'/','serial_month','/', file_pre, hcl{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_clono2i = strcat( clim_dir, clono2{n},'/','serial_month','/', file_pre, clono2{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_n2oi = strcat( clim_dir, n2o{n},'/','serial_month','/', file_pre, n2o{n}, sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            file_clytoti = strcat( clim_dir, 'clytotcmam','/','serial_month','/', file_pre, 'clytotcmam', sprintf('_%i_%0.2i',yearsin(j), months(i)), file_post);
            
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
            if exist(file_n2oi,'file') ~= 2 || isempty(n2o{n})
                fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', n2o{n}, yearsin(j), months(i))
                vmrzon_n2o(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                vmrzon_n2o_error(:,:,ij) = nan(48,36);
                vmrzon_n2o_obscount(:,:,ij) = nan(48,36);
            else
                clim = load(file_n2oi); clim = clim.climstruct; % the variable is called climstruct in the new data
                vmrzon_n2o(:,:,ij) = clim.vmr_zonal;
                vmrzon_n2o_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                vmrzon_n2o_obscount(:,:,ij) = clim.obs_count;
            end
            if n == 2
                if exist(file_clytoti,'file') ~= 2
                    fprintf('There is no file for %s_%i_%0.2i. Moving on...\n', 'clytotcmam', yearsin(j), months(i))
                    vmrzon_clytot(:,:,ij) = nan(48,36); % this is hardcoded here, for now.
                    vmrzon_clytot_error(:,:,ij) = nan(48,36);
                    vmrzon_clytot_obscount(:,:,ij) = nan(48,36);
                else
                    clim = load(file_clytoti); clim = clim.climstruct; % the variable is called climstruct in the new data
                    vmrzon_clytot(:,:,ij) = clim.vmr_zonal;
                    vmrzon_clytot_error(:,:,ij) = sqrt(clim.vmr_zonal_var);
                    vmrzon_clytot_obscount(:,:,ij) = clim.obs_count;
                end
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
    vmrzon_hocl_error_nan2zero = vmrzon_hocl_error;
    vmrzon_clono2_error_nan2zero = vmrzon_clono2_error;
    [Jnanprofile_hocl, Knanprofile_hocl] = find(squeeze(nansum(vmrzon_hocl_nan2zero,1) == 0));
    for i = 1:length(Jnanprofile_hocl)
        vmrzon_hocl_nan2zero(:,Jnanprofile_hocl(i),Knanprofile_hocl(i)) = 999; % change all-nan profiles to all-999 profiles
        vmrzon_hocl_error_nan2zero(:,Jnanprofile_hocl(i),Knanprofile_hocl(i)) = 999; % change all-nan profiles to all-999 profiles
    end
    [Jnanprofile_clono2, Knanprofile_clono2] = find(squeeze(nansum(vmrzon_clono2_nan2zero,1) == 0));
    for i = 1:length(Jnanprofile_clono2)
        vmrzon_clono2_nan2zero(:,Jnanprofile_clono2(i),Knanprofile_clono2(i)) = 999; % change all-nan profiles to all-999 profiles
        vmrzon_clono2_error_nan2zero(:,Jnanprofile_clono2(i),Knanprofile_clono2(i)) = 999; % change all-nan profiles to all-999 profiles
    end
    %hocl: assume missing values values are negligible
    vmrzon_hocl_nan2zero(isnan(vmrzon_hocl_nan2zero)) = 0;
    vmrzon_hocl_error_nan2zero(isnan(vmrzon_hocl_error_nan2zero)) = 0;
    %clono2: upper scaled a priori sucks so we cant use it. assume missing
    %values above 30km are negligible.
    dummy = vmrzon_clono2_nan2zero(izace_30up,:,:);
    dummy(isnan(dummy)) = 0;
    vmrzon_clono2_nan2zero(izace_30up,:,:) = dummy;
    dummy = vmrzon_clono2_error_nan2zero(izace_30up,:,:);
    dummy(isnan(dummy)) = 0;
    vmrzon_clono2_error_nan2zero(izace_30up,:,:) = dummy;
    % restore the 999 values to nans.
    vmrzon_hocl_nan2zero(vmrzon_hocl_nan2zero == 999) = nan;
    vmrzon_clono2_nan2zero(vmrzon_clono2_nan2zero == 999) = nan;
    vmrzon_hocl_error_nan2zero(vmrzon_hocl_error_nan2zero == 999) = nan;
    vmrzon_clono2_error_nan2zero(vmrzon_clono2_error_nan2zero == 999) = nan;
    
    %% make Cly
    %decide which types of data to use in the sum
    cly_clo = vmrzon_clo;
    cly_hocl = vmrzon_hocl;
    cly_hcl = vmrzon_hcl;
    cly_clono2 = vmrzon_clono2;
    
    cly_clo_edit = vmrzon_clo;
    cly_hocl_edit = vmrzon_hocl_nan2zero; % ignore missing data except for when there is no data in a profile at all
    cly_hcl_edit = vmrzon_hcl; % don't ignore missing data because it is always relevent
    cly_clono2_edit = vmrzon_clono2_nan2zero; % ignore missing data except for when there is no data in a profile at all
    
    cly_clo_error_edit = vmrzon_clo_error;
    cly_hocl_error_edit = vmrzon_hocl_error_nan2zero; % ignore missing data except for when there is no data in a profile at all
    cly_hcl_error_edit = vmrzon_hcl_error; % don't ignore missing data because it is always relevent
    cly_clono2_error_edit = vmrzon_clono2_error_nan2zero; % ignore missing data except for when there is no data in a profile at all
    
    cly_clo_error = vmrzon_clo_error;
    cly_hocl_error = vmrzon_hocl;
    cly_hcl_error = vmrzon_hcl_error;
    cly_clono2_error = vmrzon_clono2_error;
    
    vmrzon_cly = cly_clo + cly_hocl + cly_hcl + cly_clono2;
    vmrzon_cly_error = sqrt(cly_clo_error.^2 + cly_hocl_error.^2 + cly_hcl_error.^2 + cly_clono2_error.^2);
    
    vmrzon_cly_edit = cly_clo_edit + cly_hcl_edit + cly_clono2_edit + cly_hocl_edit;
    vmrzon_cly_error_edit = sqrt(cly_clo_error_edit.^2 + cly_hcl_error_edit.^2 + cly_clono2_error_edit.^2 + cly_hocl_error_edit.^2);
% % %     vmrzon_cly_edit = vmrzon_hcl ./ (cly_clo_edit + cly_hcl_edit + cly_clono2_edit + cly_hocl_edit);
% % %     vmrzon_cly_error_edit = vmrzon_hcl_error ./ sqrt(cly_clo_error_edit.^2 + cly_hcl_error_edit.^2 + cly_clono2_error_edit.^2 + cly_hocl_error_edit.^2);
    %%
    %     vmrzon_hocl = vmrzon_hocl_nan2zero;
    %%
    
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
    
    vmrzon_clo = squeeze(nanmean(vmrzon_clo(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hocl = squeeze(nanmean(vmrzon_hocl(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hcl = squeeze(nanmean(vmrzon_hcl(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clono2 = squeeze(nanmean(vmrzon_clono2(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly = squeeze(nanmean(vmrzon_cly(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_edit = squeeze(nanmean(vmrzon_cly_edit(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clytot = squeeze(nanmean(vmrzon_clytot(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_n2o = squeeze(nanmean(vmrzon_n2o(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    
    vmrzon_clo_error = squeeze(nanmean(vmrzon_clo_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hocl_error = squeeze(nanmean(vmrzon_hocl_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_hcl_error = squeeze(nanmean(vmrzon_hcl_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clono2_error = squeeze(nanmean(vmrzon_clono2_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_error = squeeze(nanmean(vmrzon_cly_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_cly_edit_error = squeeze(nanmean(vmrzon_cly_error_edit(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_clytot_error = squeeze(nanmean(vmrzon_clytot_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    vmrzon_n2o_error = squeeze(nanmean(vmrzon_n2o_error(:,ilatmin:ilatmax,:), 2 )) * 1e9;
    
    %% get the seasonal data
    [ vmrzon_cly_edit_seasonal_average, vmrzon_cly_edit_seasonal_timeseries, vmrzon_cly_edit_seasonal_timeseries_detrended, mid_dates_mjd  ] = get_seasonality_2d( vmrzon_cly_edit, datenum2mjd(sdates));
    [ vmrzon_clo_seasonal_average, vmrzon_clo_seasonal_timeseries, vmrzon_clo_seasonal_timeseries_detrended] = get_seasonality_2d( vmrzon_clo, datenum2mjd(sdates));
    [ vmrzon_hocl_seasonal_average, vmrzon_hocl_seasonal_timeseries, vmrzon_hocl_seasonal_timeseries_detrended] = get_seasonality_2d( vmrzon_hocl, datenum2mjd(sdates));
    [ vmrzon_hcl_seasonal_average, vmrzon_hcl_seasonal_timeseries, vmrzon_hcl_seasonal_timeseries_detrended] = get_seasonality_2d( vmrzon_hcl, datenum2mjd(sdates));
    [ vmrzon_clono2_seasonal_average, vmrzon_clono2_seasonal_timeseries, vmrzon_clono2_seasonal_timeseries_detrended ] = get_seasonality_2d( vmrzon_clono2, datenum2mjd(sdates));
    [ vmrzon_cly_seasonal_average, vmrzon_cly_seasonal_timeseries, vmrzon_cly_seasonal_timeseries_detrended ] = get_seasonality_2d( vmrzon_cly, datenum2mjd(sdates));
    [ vmrzon_clytot_seasonal_average, vmrzon_clytot_seasonal_timeseries, vmrzon_clytot_seasonal_timeseries_detrended ] = get_seasonality_2d( vmrzon_clytot, datenum2mjd(sdates));
    [ vmrzon_n2o_seasonal_average, vmrzon_n2o_seasonal_timeseries, vmrzon_n2o_seasonal_timeseries_detrended ] = get_seasonality_2d( vmrzon_n2o, datenum2mjd(sdates));
    
    [ ~, vmrzon_cly_edit_error_seasonal_timeseries  ] = get_seasonality_2d( vmrzon_cly_edit_error, datenum2mjd(sdates));
    [ ~, vmrzon_clo_error_seasonal_timeseries] = get_seasonality_2d( vmrzon_clo_error, datenum2mjd(sdates));
    [ ~, vmrzon_hocl_error_seasonal_timeseries] = get_seasonality_2d( vmrzon_hocl_error, datenum2mjd(sdates));
    [ ~, vmrzon_hcl_error_seasonal_timeseries] = get_seasonality_2d( vmrzon_hcl_error, datenum2mjd(sdates));
    [ ~, vmrzon_clono2_error_seasonal_timeseries] = get_seasonality_2d( vmrzon_clono2_error, datenum2mjd(sdates));
    [ ~, vmrzon_cly_error_seasonal_timeseries] = get_seasonality_2d( vmrzon_cly_error, datenum2mjd(sdates));
    [ ~, vmrzon_clytot_error_seasonal_timeseries ] = get_seasonality_2d( vmrzon_clytot_error, datenum2mjd(sdates));
    [ ~, vmrzon_n2o_error_seasonal_timeseries ] = get_seasonality_2d( vmrzon_n2o_error, datenum2mjd(sdates));
    vmrzon_clo_error_seasonal_timeseries(vmrzon_clo_error_seasonal_timeseries == 0) = nan;
    vmrzon_hcl_error_seasonal_timeseries(vmrzon_hcl_error_seasonal_timeseries == 0) = nan;
    vmrzon_hocl_error_seasonal_timeseries(vmrzon_hocl_error_seasonal_timeseries == 0) = nan;
    vmrzon_clono2_error_seasonal_timeseries(vmrzon_clono2_error_seasonal_timeseries == 0) = nan;
    vmrzon_cly_error_seasonal_timeseries(vmrzon_cly_error_seasonal_timeseries == 0) = nan;
    vmrzon_cly_edit_error_seasonal_timeseries(vmrzon_cly_edit_error_seasonal_timeseries == 0) = nan;
    vmrzon_clytot_error_seasonal_timeseries(vmrzon_clytot_error_seasonal_timeseries == 0) = nan;
    
    %% remove the +2.8% per decade linear trend in N2O
    % find the quarterly increase of i * 0.07/100 * [n2o]_1 and subtract
    i_n2o = 0 : length(vmrzon_n2o_seasonal_timeseries_detrended(1,:)) - 1; % a column vector
    i_n2o = repmat(i_n2o, [length(vmrzon_n2o_seasonal_timeseries_detrended(:,1)), 1]); % repeat the vector to match the size of the n2o array
    trend07 = (0.07/100) * i_n2o .* vmrzon_n2o_seasonal_timeseries_detrended(:,1); % mulitply by the values in the first column of the timeseries
    vmrzon_n2o_seasonal_timeseries_detrended_07 = vmrzon_n2o_seasonal_timeseries_detrended - trend07;
    
    %% get the fits to the data using [Y] = a + b*[X] + c*[n2o]
    warning off
    fit_clo = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_clo_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_clo_error_seasonal_timeseries );
    fit_hocl = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_hocl_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_hocl_error_seasonal_timeseries );
    fit_hcl = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_hcl_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_hcl_error_seasonal_timeseries );
    fit_clono2 = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_clono2_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_clono2_error_seasonal_timeseries );
    fit_cly_edit = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_cly_edit_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_cly_edit_error_seasonal_timeseries );
%     fit_cly = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_cly_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_cly_error_seasonal_timeseries );
    fit_clytot = get_trend_2d_poly11( mid_dates_mjd, vmrzon_n2o_seasonal_timeseries_detrended_07, vmrzon_clytot_seasonal_timeseries_detrended, 1, 1 ./ vmrzon_clytot_error_seasonal_timeseries );
    warning on
    %%
    
    %% make the line of the data using the a and b coefficients
    sdates = mjd2datenum(mid_dates_mjd);
    lw = 1;
    ms = 4;
    fs = 16;
    for i = 1 : nalt
        figure(figi),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = blue;
        else
            ci = red;
        end
        X0X1 = fit_clo.p00(ipplot(i), 1) + fit_clo.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_clo.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_clo_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_clo_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'color', ci), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci)
        plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci)
        
        approx_trend = [(fit_clo.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_clo.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('ClO trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')
    %
    figp = figi + 1;
    for i = 1 : nalt
        figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = blue;
        else
            ci = red;
        end
        X0X1 = fit_hocl.p00(ipplot(i), 1) + fit_hocl.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_hocl.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_hocl_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_hocl_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'color', ci), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci)
        plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci)
        
        approx_trend = [(fit_hocl.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_hocl.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('HOCl trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')
    %
    figp = figp + 1;
    for i = 1 : nalt
        figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = blue;
        else
            ci = red;
        end
        X0X1 = fit_clono2.p00(ipplot(i), 1) + fit_clono2.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_clono2.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_clono2_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_clono2_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'color', ci), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci)
        plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci)
        
        approx_trend = [(fit_clono2.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_clono2.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('ClONO2 trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')
    %
    figp = figp + 1;
    for i = 1 : nalt
        figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = 'k';
        else
            ci = red;
        end
        X0X1 = fit_hcl.p00(ipplot(i), 1) + fit_hcl.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_hcl.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_hcl_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_hcl_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'o-', 'color', ci, 'Linewidth', lw, 'MarkerSize', ms), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci, 'Linewidth', lw)
%         plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci, 'Linewidth', lw)
        
        approx_trend = [(fit_hcl.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_hcl.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('HCl trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')
    set(gca, 'FontSize',fs);
% % %     %
% % %     figp = figp + 1;
% % %     for i = 1 : nalt
% % %         figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
% % %         if n == 1
% % %             ci = blue;
% % %         else
% % %             ci = red;
% % %         end
% % %         X0X1 = fit_cly.p00(ipplot(i), 1) + fit_cly.p10(ipplot(i), 1) .* mid_dates_mjd;
% % %         X2 = fit_cly.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
% % %         X0X1X2 = X0X1 + X2;
% % %         total_mean = nanmean(vmrzon_cly_seasonal_average(ipplot(i),:), 2);
% % %         plot(sdates, vmrzon_cly_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'color', ci), hold on
% % %         plot(sdates,  X0X1 + total_mean, 'color', ci)
% % %         plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
% % %         plot(sdates, X2 + total_mean, '-.', 'color', ci)
% % %         
% % %         approx_trend = [(fit_cly.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_cly.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
% % %         text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
% % %         text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
% % %         title('Cly trend')
% % %         ylabel('VMR [ppb]')
% % %         xlabel('date')
% % %         
% % %     end
% % %     dynamicDateTicks([],'x','mm')
    %
    figp = figp + 1;
    for i = 1 : nalt
        figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = 'k';
        else
            ci = red;
        end
        X0X1 = fit_cly_edit.p00(ipplot(i), 1) + fit_cly_edit.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_cly_edit.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_cly_edit_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_cly_edit_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'o-', 'color', ci, 'Linewidth', lw, 'MarkerSize', ms), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci, 'Linewidth', lw)
%         plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci, 'Linewidth', lw)
        
        approx_trend = [(fit_cly_edit.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_cly_edit.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('Cly_e_d_i_t trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')
    set(gca, 'FontSize',fs);
    %
    figp = figp + 1;
    for i = 1 : nalt
        figure(figp),% colorOrder = get(gca, 'ColorOrder'); ci = colorOrder(i, :);
        if n == 1
            ci = blue;
        else
            ci = red;
        end
        X0X1 = fit_clytot.p00(ipplot(i), 1) + fit_clytot.p10(ipplot(i), 1) .* mid_dates_mjd;
        X2 = fit_clytot.p01(ipplot(i), 1) .* vmrzon_n2o_seasonal_timeseries_detrended_07(ipplot(i), :);
        X0X1X2 = X0X1 + X2;
        total_mean = nanmean(vmrzon_clytot_seasonal_average(ipplot(i),:), 2);
        plot(sdates, vmrzon_clytot_seasonal_timeseries_detrended(ipplot(i), :) + total_mean, 'color', ci), hold on
        plot(sdates,  X0X1 + total_mean, 'color', ci)
        plot(sdates, X0X1X2 + total_mean, '--', 'color', ci)
        plot(sdates, X2 + total_mean, '-.', 'color', ci)
        
        approx_trend = [(fit_clytot.p10(ipplot(i), 1) * (365.25*10)), 100 * (fit_clytot.p10(ipplot(i), 1) * (365.25*10)) ./ total_mean]
        text(sdates(end), X0X1(end) + total_mean, sprintf('%0.1f hPa (~%0.1f km)', pace(ipplot(i)), zace(ipplot(i)) ), 'Color', ci);
        text(sdates(1), X0X1(1) + 0.9 * total_mean, sprintf('%0.2fppb, %0.2f%%', approx_trend(1), approx_trend(2) ), 'Color', ci);
        title('Cly_t_o_t_a_l trend')
        ylabel('VMR [ppb]')
        xlabel('date')
        
    end
    dynamicDateTicks([],'x','mm')

% % 
% %     % test = fit_hcl;
% %     cly_edit_trend = 100 * (fit_cly_edit.p10 * (365.25*10)) ./ nanmean(vmrzon_cly_edit_seasonal_average, 2);
% % %     cly_trend = (fit_cly.p10 * (365.25*10)) ./ nanmean(vmrzon_cly_seasonal_average, 2);
% %     clytot_trend = 100 * (fit_clytot.p10 * (365.25*10)) ./ nanmean(vmrzon_clytot_seasonal_average, 2);
% %     hcl_trend = 100 * (fit_hcl.p10 * (365.25*10)) ./ nanmean(vmrzon_hcl_seasonal_average, 2);
%     test.x = cly_trend;
    if n == 1
        trend.ace_clo = fit_clo.p10 * (365.25*10);
        trend.ace_clo_p = 100 * trend.ace_clo ./ nanmean(vmrzon_clo_seasonal_average, 2);
        trend.ace_hocl = fit_hocl.p10 * (365.25*10);
        trend.ace_hocl_p = 100 * trend.ace_hocl ./ nanmean(vmrzon_hocl_seasonal_average, 2);
        trend.ace_clono2 = fit_clono2.p10 * (365.25*10);
        trend.ace_clono2_p = 100 * trend.ace_clono2 ./ nanmean(vmrzon_clono2_seasonal_average, 2);
        trend.ace_hcl = fit_hcl.p10 * (365.25*10);
        trend.ace_hcl_p = 100 * trend.ace_hcl ./ nanmean(vmrzon_hcl_seasonal_average, 2);
        trend.ace_cly = fit_cly_edit.p10 * (365.25*10);
        trend.ace_cly_p = 100 * trend.ace_cly ./ nanmean(vmrzon_cly_edit_seasonal_average, 2);
    else
        trend.cmam_clo = fit_clo.p10 * (365.25*10);
        trend.cmam_clo_p = 100 * trend.cmam_clo ./ nanmean(vmrzon_clo_seasonal_average, 2);
        trend.cmam_hocl = fit_hocl.p10 * (365.25*10);
        trend.cmam_hocl_p = 100 * trend.cmam_hocl ./ nanmean(vmrzon_hocl_seasonal_average, 2);
        trend.cmam_clono2 = fit_clono2.p10 * (365.25*10);
        trend.cmam_clono2_p = 100 * trend.cmam_clono2 ./ nanmean(vmrzon_clono2_seasonal_average, 2);
        trend.cmam_hcl = fit_hcl.p10 * (365.25*10);
        trend.cmam_hcl_p = 100 * trend.cmam_hcl ./ nanmean(vmrzon_hcl_seasonal_average, 2);
        trend.cmam_cly = fit_cly_edit.p10 * (365.25*10);
        trend.cmam_cly_p = 100 * trend.cmam_cly ./ nanmean(vmrzon_cly_edit_seasonal_average, 2);
    end
    trend.pgrid = pace;
    trend.zgrid = zace;
    %     test.y = cly_edit_trend;
    %     test.z = clytot_trend;
    
end
%% trend plots
lw = 2;
ms = 5;
fs = 16;
ylimlow = 17;
ylimhigh = 31;
% absolute trend
confi_ace = abs(trend.ace_clo(:,1) - trend.ace_clo(:,2));
confi_cmam = abs(trend.cmam_clo(:,1) - trend.cmam_clo(:,2));
figure, set(gcf, 'Position', [26         500        1217         474])
subplot(1, 2, 1), [hl,hp] = boundedline(trend.ace_clo(:,1), pace, confi_ace,'-ko', ...
                      trend.cmam_clo(:,1), pace, confi_cmam, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('ClO trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [ppb / decade]')
% percentage trend
confi_ace_p = abs(trend.ace_clo_p(:,1) - trend.ace_clo_p(:,2));
confi_cmam_p = abs(trend.cmam_clo_p(:,1) - trend.cmam_clo_p(:,2));
% figure
subplot(1, 2, 2), [hl,hp] = boundedline(trend.ace_clo_p(:,1), pace, confi_ace_p,'-ko', ...
                      trend.cmam_clo_p(:,1),pace, confi_cmam_p, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('ClO trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [% / decade]')

% absolute trend
confi_ace = abs(trend.ace_hocl(:,1) - trend.ace_hocl(:,2));
confi_cmam = abs(trend.cmam_hocl(:,1) - trend.cmam_hocl(:,2));
figure, set(gcf, 'Position', [26         500        1217         474])
subplot(1, 2, 1), [hl,hp] = boundedline(trend.ace_hocl(:,1), pace, confi_ace,'-ko', ...
                      trend.cmam_hocl(:,1), pace, confi_cmam, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('HOCl trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [ppb / decade]')
% percentage trend
confi_ace_p = abs(trend.ace_hocl_p(:,1) - trend.ace_hocl_p(:,2));
confi_cmam_p = abs(trend.cmam_hocl_p(:,1) - trend.cmam_hocl_p(:,2));
subplot(1, 2, 2), [hl,hp] = boundedline(trend.ace_hocl_p(:,1), pace, confi_ace_p,'-ko', ...
                      trend.cmam_hocl_p(:,1),pace, confi_cmam_p, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('HOCl trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [% / decade]')

% absolute trend
confi_ace = abs(trend.ace_clono2(:,1) - trend.ace_clono2(:,2));
confi_cmam = abs(trend.cmam_clono2(:,1) - trend.cmam_clono2(:,2));
figure, set(gcf, 'Position', [26         500        1217         474])
subplot(1, 2, 1), [hl,hp] = boundedline(trend.ace_clono2(:,1), pace, confi_ace,'-ko', ...
                      trend.cmam_clono2(:,1), pace, confi_cmam, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('ClONO_2 trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [ppb / decade]')
% percentage trend
confi_ace_p = abs(trend.ace_clono2_p(:,1) - trend.ace_clono2_p(:,2));
confi_cmam_p = abs(trend.cmam_clono2_p(:,1) - trend.cmam_clono2_p(:,2));
subplot(1, 2, 2), [hl,hp] = boundedline(trend.ace_clono2_p(:,1), pace, confi_ace_p,'-ko', ...
                      trend.cmam_clono2_p(:,1),pace, confi_cmam_p, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('ClONO_2 trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [% / decade]')

% absolute trend
confi_ace = abs(trend.ace_hcl(:,1) - trend.ace_hcl(:,2));
confi_cmam = abs(trend.cmam_hcl(:,1) - trend.cmam_hcl(:,2));
figure, set(gcf, 'Position', [26         500        1217         474])
subplot(1, 2, 1), [hl,hp] = boundedline(trend.ace_hcl(:,1), pace, confi_ace,'-ko', ...
                      trend.cmam_hcl(:,1), pace, confi_cmam, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('HCl trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [ppb / decade]')
% percentage trend
confi_ace_p = abs(trend.ace_hcl_p(:,1) - trend.ace_hcl_p(:,2));
confi_cmam_p = abs(trend.cmam_hcl_p(:,1) - trend.cmam_hcl_p(:,2));
subplot(1, 2, 2), [hl,hp] = boundedline(trend.ace_hcl_p(:,1), pace, confi_ace_p,'-ko', ...
                      trend.cmam_hcl_p(:,1),pace, confi_cmam_p, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('HCl trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(ylimhigh) pace(ylimlow)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(ylimlow) zace(ylimhigh)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [% / decade]')

% absolute trend
confi_ace = abs(trend.ace_cly(:,1) - trend.ace_cly(:,2));
confi_cmam = abs(trend.cmam_cly(:,1) - trend.cmam_cly(:,2));
figure, set(gcf, 'Position', [26         500        1217         474])
subplot(1, 2, 1), [hl,hp] = boundedline(trend.ace_cly(:,1), pace, confi_ace,'-ko', ...
                      trend.cmam_cly(:,1), pace, confi_cmam, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('Cl_y trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(31) pace(17)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(17) zace(31)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [ppb / decade]')



% percentage trend
confi_ace_p = abs(trend.ace_cly_p(:,1) - trend.ace_cly_p(:,2));
confi_cmam_p = abs(trend.cmam_cly_p(:,1) - trend.cmam_cly_p(:,2));
subplot(1, 2, 2), [hl,hp] = boundedline(trend.ace_cly_p(:,1), pace, confi_ace_p,'-ko', ...
                      trend.cmam_cly_p(:,1),pace, confi_cmam_p, '-ro', ...
                      'nan', 'gap', 'orientation', 'horiz', 'alpha');
title('Cl_y trend 60S - 60N')
yyaxis left
set(gca, 'Ydir','reverse', 'YScale', 'log', 'XScale', 'linear', 'YMinorTick','on','YMinorGrid','OFF', 'XMinorGrid', 'ON', 'FontSize',fs);
%             ylabel('pressure [hPa]')
ylim([pace(31) pace(17)]);
ylabel('pressure [hPa]')
yyaxis right
ylim([zace(17) zace(31)]);
ylabel('approximate altitude [km]')
xlabel('\Delta VMR [% / decade]')


%
end

