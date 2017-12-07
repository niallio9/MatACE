function [ ] = make_ace_climatology( tanstruct, out_directory)
%A function to create zonally averaged climatologies of ACE measurements.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           climstruct: STRUCTURE - with fields corresponding to
%           zonally averaged ace data. The fields of the structure are
%           listed in (A), below.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/17

tic
%% Things that may be changed often
if nargin < 2
    home_linux = '/home/niall/Dropbox/climatology/nryan/'; %#ok<NASGU>
    home_mac = '/Users/niall/Dropbox/climatology/nryan/';
    climdirectory = strcat(home_mac,'matclim/');
else
    climdirectory = out_directory;
end
% 48 vertical (adding mesosphere). These grid levels were taken from Jaho's
% code
pgrid = [1000 850 700 500 400 300 250 200 170 150 130 115 100 90 80 70 50 30 20 ...
    15 10 7 5 3 2 1.5 1 0.7 0.5 0.3 0.2 0.15 0.1 0.08 0.05 0.03 0.02 0.01 0.007 0.004 ...
    0.003 0.002 0.001 0.0008 0.0005 0.0003 0.0002 0.0001]';
% Standard 5 degree latitude bins for for ACE zonal climatology
lat_bounds = -90:5:90;
%the name of the output files
savename_pre = 'ACEFTS_CLIM_v3_';
% cells with the names of the month
monthnames = {'January', 'February', 'March', 'April', 'May', 'June', ...
    'July', 'August', 'September', 'October', 'November', 'December'};

%% Define some things
gas = tanstruct;
climstruct = []; %#ok<NASGU>

%% Only continue of the output directory exists
if isdir(climdirectory)
    
    %% Apply the flags to the ACE data
    fprintf('\nApplying the data flags...\n')
    [ gas ] = apply_ace_flags(gas);
    
    %% Interpolate the data and dmps to the model pressure grid (defined above)
    % fprintf('\nInterpolating the gas data to the pre-defined pressure grid...')
    [ gas ] = interpolate_ace_to_pgrid( gas, pgrid );
    % fprintf('Done\n')
    
    %% Loop through the months and create climatologies for each.
    %Within the loop, separate out the contributions to each latitude band
    %defined above
    fprintf('\nPreparing climatology by month...')
    for i = 1:12
        %subset the ace data and DMPs by month
        warning off % supress warnings about reducing the data to zero here. There is output below if this is the case
        gas_monthi = subset_ace_by_month(gas,i);
        warning on
        %bin the monthly data by latitude. this creates a structure with new fields that are explained below in (A).
        fprintf('\nBinning the %s data by latitude\n', monthnames{i});
        gas_monthi_latbin  = bin_ace_by_lat( gas_monthi, lat_bounds );
        % add
        
        % save the file as a matlab structure for now. In the chosen directory
        if  nansum(gas_monthi_latbin.date_mjd_mean) ~= 0
            climstruct = gas_monthi_latbin; %#ok<NASGU>
            savename_post = sprintf('_%02d',i);
            savedest = strcat(climdirectory,savename_pre,gas.gas,savename_post);
            fprintf('Saving %s climatology to %s\n', monthnames{i}, savedest);
            save(savedest,'climstruct');
        else
            fprintf('There are no climatology data for %s\n', monthnames{i});
        end
    end
%     testout = climstruct;
    fprintf('\nDone :)\n')
    %
    toc
else
    fprintf('\nThe directory, %s, does not exist. Please create it first\n', climdirectory)
end
end

%% (A)
% fieldnames of the structure produced by 'bin_ace_by_lat'
% 
%     'source_file'
%     'beta_angle_mean'
%     'beta_angle_var'
%     'date_mjd_mean'
%     'gas'
%     'lon_tangent_mean'
%     'lat_tangent_mean'
%     'pressure_hPa'
%     'altitude_km_mean'
%     'vmr_zonal'
%     'vmr_zonal_var'
%     'vmr_zonal_error'
%     'vmr_zonal_standard_error'
%     'lat'
%     'lat_bounds'
%     'obs_count'
%     'obs_location'
%     'lst_median'
%     'lst_mean'
%     'lst_max'
%     'lst_min'
%     'lst_var'
