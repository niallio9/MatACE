function [ tanstruct_out ] = merge_ace_clim_data( climstruct1_in, climstruct2_in )
%A function to merge two ace climatology files so that the data is a
%combination of the data in both files.

% *INPUT*
%           climstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with
%           'make_ace_climatology_serialmonth.m' or
%           'make_ace_climatology_multiple.m'
%
%           climstruct2_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with
%           'make_ace_climatology_serialmonth.m' or
%           'make_ace_climatology_multiple.m'
%
% *OUTPUT*
%           climstruct_out: STRUCTURE - output has the same fields as the
%           input, but with a combination of the non-intersecting data in
%           the two input files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 07/2018
%% Define some things
gas1 = climstruct1_in;
gas2 = climstruct2_in;
gasout = [];
% time_dim2 = 2;
time_dim3 = 3; % the dimension in which the lat x alt arrays will be concatenated. It will correspond to time
%% check if one of the inputs is empty
if isempty(gas1) && isempty(gas2) % if the variables are empty
    fprintf('both of the inputs are empty. stoppping...\n')
    tanstruct_out = [];
    return
elseif isempty(gas2)
    disp('the second input is empty, returning the first')
    tanstruct_out = climstruct1_in;
    return
elseif isempty(gas1)
    disp('the first input is empty, returning the second')
    tanstruct_out = climstruct2_in;
    return
elseif isempty(gas1.start_date) && isempty(gas2.time) % if the structures have no occultations
    fprintf('both of the input structures are empty. stoppping...\n')
    tanstruct_out = [];
    return
elseif isempty(gas2.start_date)
    disp('the second input structure is empty, returning the first')
    tanstruct_out = climstruct1_in;
    return
elseif isempty(gas1.start_date)
    disp('the first input structure is empty, returning the second')
    tanstruct_out = climstruct2_in;
    return
end
% so both structures should contain data from here
lalt1 = length(gas1.pressure_hPa);
lalt2 = length(gas2.pressure_hPa);
llat1 = length(gas1.lat);
llat2 = length(gas2.lat);
%% check if the gasnames of the two file match and the latitudes and pressures
if ~strcmp(gas1.gas, gas2.gas)
    error('The gas names for each input file don''t match')
end
if isfield(gas1, 'climatology_type') && isfield(gas2, 'climatology_type')
    if ~strcmp(gas1.climatology_type, gas2.climatology_type)
        error('The climatology_type for each input file doesn''t match')
    end
end
if ~isequal(gas1.pressure_hPa, gas2.pressure_hPa)
    error('The pressures for each input file don''t match')
end
if ~isequal(gas1.lat, gas2.lat)
    error('The latitudes for each input file don''t match')
end

%% merge the data
if isfield(gas1,'source_file') && isfield(gas2,'source_file') 
    gasout.source_file = strcat(gas1.source_file, ', ', gas2.source_file);
else
    fprintf('the ''source_file'' fields were not merged because it is not present in both inputs\n')
end
gasout.beta_angle_mean = [gas1.beta_angle_mean; gas2.beta_angle_mean];
gasout.beta_angle_var = [gas1.beta_angle_var; gas2.beta_angle_var];
gasout.doy_mean = [gas1.doy_mean; gas2.doy_mean];
gasout.start_date = [gas1.start_date; gas2.start_date];
gasout.end_date = [gas1.end_date; gas2.end_date];
gasout.gas = gas1.gas;
gasout.lat_tangent_mean = [gas1.lat_tangent_mean; gas2.lat_tangent_mean];
gasout.lon_tangent_mean = [gas1.lon_tangent_mean; gas2.lon_tangent_mean];
gasout.altitude_km_mean = gas1.altitude_km_mean;
gasout.pressure_hPa = gas1.pressure_hPa;
gasout.vmr_zonal = cat(time_dim3, gas1.vmr_zonal, gas2.vmr_zonal);
gasout.vmr_zonal_var = cat(time_dim3, gas1.vmr_zonal_var, gas2.vmr_zonal_var);
gasout.vmr_zonal_error = cat(time_dim3, gas1.vmr_zonal_error, gas2.vmr_zonal_error);
gasout.vmr_zonal_standard_error = cat(time_dim3, gas1.vmr_zonal_standard_error, gas2.vmr_zonal_standard_error);
gasout.lat = gas1.lat;
gasout.lat_bounds = gas1.lat_bounds;
gasout.lon_mean = cat(time_dim3, gas1.lon_mean, gas2.lon_mean);
gasout.obs_count = cat(time_dim3, gas1.obs_count, gas2.obs_count);
gasout.obs_location = cat(time_dim3, gas1.obs_location, gas2.obs_location);
gasout.lst_median = [gas1.lst_median; gas2.lst_median];
gasout.lst_mean = [gas1.lst_mean; gas2.lst_mean];
gasout.lst_max = [gas1.lst_max; gas2.lst_max];
gasout.lst_min = [gas1.lst_min; gas2.lst_min];
gasout.lst_var = [gas1.lst_var; gas2.lst_var];
if isfield(gas1, 'climatology_type') && isfield(gas2, 'climatology_type')
    gasout.climatology_type = gas1.climatology_type;
else
    fprintf('the ''climatology_type'' fields were not merged because it is not present in both inputs\n')
end
if isfield(gas1, 'time') && isfield(gas2, 'time')
    gasout.time = [gas1.time; gas2.time];
else
    fprintf('the ''time'' fields were not merged because it is not present in both inputs\n')
end


% %% remove duplicate occultations
% disp('removing duplicate occultations...')
% % orbitnames = get_ace_occultation_names(gasout);
% %the next line also automatically sorts the data
% [~,igood] = unique(gasout.date_mjd(1,:)); % included an indexing here because the sampled MLS data has an 'alt x occ' array of times 
% % size(igood)
% gasout = reduce_tanstruct_by_rowindex(gasout, igood);
% disp(gasout)

%%
tanstruct_out = gasout;
%
end

