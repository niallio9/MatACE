function [ vmrzon_cly_error, vmrzon_hcl_error ] = make_ace_cly_climatology( appendix )
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
file_post = strcat('_', appendix, '.mat');

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
        clo = {'ClOmlsfrac10lim4ppb_sap'};
        hocl = {'HOClmlsfrac10lim4ppb_sap'};
        hcl = {'HCl'};
        clono2 = {'ClONO2'};
        %         cly = 'ClOy';
    case 'model'
        clo = {'ClOcmam'};
        hocl = {'HOClcmam'};
        hcl = {'HClcmam'};
        clono2 = {'ClONO2cmam'};
        %         cly = 'Clycmam' % uses more than the 4 gases listed above
    case 'both'
        clo = {'ClOmlsfrac10lim4ppb_sap','ClOcmam'};
        hocl = {'HOClmlsfrac10lim4ppb_sap','HOClcmam'};
        hcl = {'HCl','HClcmam'};
        clono2 = {'ClONO2','ClONO2cmam'};
end

%% define some things
home_linux = '/home/niall/Dropbox/climatology/'; %#ok<NASGU>
home_mac = '/Users/niall/Dropbox/climatology/'; %#ok<NASGU>
home_windows = 'C:\Users\ryann\jaho\'; %#ok<NASGU>
% newclim_dir = strcat(home_windows,'climdata_v3p5_nr\');
clim_dir = 'C:\Users\ryann\ACE\climdata_testing\time_matched_climatology\';
% newclim_dir = 'C:\Users\ryann\MLS\climdata\';
% newclim_dir = 'C:\Users\ryann\ACE\MAESTRO\climdata\';


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
    
    %%  make and save the climstruct
    clim = load(file_clo); clim = clim.climstruct; % the variable is called climstruct in the new data
    climstruct_out = clim;
    climstruct_out.gas(1:3) = 'Cly';
    climstruct_out.vmr = vmrzon_cly;
    climstruct_out.vmr_error = vmrzon_cly_error;
    
    
    
end
%
end

