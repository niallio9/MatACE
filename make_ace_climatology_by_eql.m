function [ climstruct ] = make_ace_climatology_by_eql( tanstruct )
%A function to create zonally averaged climatologies of ACE measurements
%using equivalent latitude. 
%The function assumes that the ACE Derived Meteorological Product(DMP)data
%has been added to the gas data. This can be done using'merge_ace_dmp.m'. 

% *INPUT*
%           tanstruct: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. The DMP data must also be added to
%           the tanstruct so that it has the equivalent latitude
%           information.
%
% *OUTPUT*
%           climstruct: STRUCTURE - with fields corresponding to
%           zonally averaged ace data. The fields of the structure are
%           listed in (A), below.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 02/18

%tic
%% Things that may be changed often
% 48 vertical (adding mesosphere). These grid levels were taken from Jaho's
% code
pgrid = [1000 850 700 500 400 300 250 200 170 150 130 115 100 90 80 70 50 30 20 ...
    15 10 7 5 3 2 1.5 1 0.7 0.5 0.3 0.2 0.15 0.1 0.08 0.05 0.03 0.02 0.01 0.007 0.004 ...
    0.003 0.002 0.001 0.0008 0.0005 0.0003 0.0002 0.0001]';
% Standard 5 degree latitude bins for for ACE zonal climatology
eql_bounds = -90:5:90;

%% Define some things
gas = tanstruct;
climstruct = []; %#ok<NASGU>

%% check if the glc information is included in the input structure
if ~isfield(gas,'eql')
    error('There is no DMP equivalent latitude information present in the input ACE data');
end

%% Apply the flags to the ACE data
fprintf('\nApplying the data flags...')
[ gas ] = apply_ace_flags(gas);
fprintf('done')

% only perform the rest of the filtering if there is data left after
% applying the flags
if isempty(gas.occultation) == 0
%     %% Replace erroneaous GLC latitude data with the tangent latitude
%     [ gas ] = filter_ace_bad_lat(gas, 10); % using a limit of 10 for the difference in GLC and tangent latitudes. This is a bit ad hoc.
%     [ gas ] = filter_ace_bad_lon(gas, 40); % using a limit of 40 for the difference in GLC and tangent longitudes. This is a bit ad hoc.
    %% Interpolate the data to the model pressure grid (defined above)
    % fprintf('\nInterpolating the gas data to the pre-defined pressure grid...')
    [ gas ] = interpolate_ace_to_pgrid( gas, pgrid ); % this also filters out bad pressure measurements (ones with values of zero).
    % fprintf('Done\n')
end

%% create climatology
%Within this, separate out the contributions to each equivalent latitude band
%defined above
%Bin the data by equivalent latitude. this creates a structure with new fields that are explained below in (A).
fprintf('\nBinning the data by equivalent latitude...');
gas_eqlbin  = bin_ace_by_eql( gas, eql_bounds );
fprintf('done\n')
climstruct_out = gas_eqlbin;
if  nansum(gas_eqlbin.date_mjd_mean) == 0
    fprintf('There are no climatology data for the input. Too few observations');
end
climstruct = climstruct_out;
%
%toc
end

%% (A)
% fieldnames of the structure produced by 'bin_ace_by_eql'
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
%     'eql'
%     'eql_bounds'
%     'obs_count'
%     'obs_location'
%     'lst_median'
%     'lst_mean'
%     'lst_max'
%     'lst_min'
%     'lst_var'
