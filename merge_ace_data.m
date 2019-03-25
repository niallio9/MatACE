function [ tanstruct_out ] = merge_ace_data( tanstruct1_in, tanstruct2_in )
%A function to merge two ace files so that the data is a combination of the
%data in both files. Duplicate occultations are removed

% *INPUT*
%           tanstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct2_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with a combination of the non-intersecting data in
%           the two input files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 07/2018
%% Define some things
gas1 = tanstruct1_in;
gas2 = tanstruct2_in;
gasout = [];
%% check if one of the inputs is empty
if isempty(gas1) && isempty(gas2) % if the variables are empty
    fprintf('both of the inputs are empty. stoppping...\n')
    tanstruct_out = [];
    return
elseif isempty(gas2)
    disp('the second input is empty, returning the first')
    tanstruct_out = tanstruct1_in;
    return
elseif isempty(gas1)
    disp('the first input is empty, returning the second')
    tanstruct_out = tanstruct2_in;
    return
elseif isempty(gas1.occultation) && isempty(gas2.occultation) % if the structures have no occultations
    fprintf('both of the input structures are empty. stoppping...\n')
    tanstruct_out = [];
    return
elseif isempty(gas2.occultation)
    disp('the second input structure is empty, returning the first')
    tanstruct_out = tanstruct1_in;
    return
elseif isempty(gas1.occultation)
    disp('the first input structure is empty, returning the second')
    tanstruct_out = tanstruct2_in;
    return
end
% so both structures should contain data from here
lgas1 = length(gas1.occultation);
lgas2 = length(gas2.occultation);
%% check if the gasnames of the two file match
if ~strcmp(gas1.gas, gas2.gas)
    error('The gas names for each input file don''t match')
end

%% merge the data
if isfield(gas1,'source_file') && isfield(gas2,'source_file') % when there is DMP data included in the tanstruct
    gasout.source_file = sstrcat(gas1.source_file, ', ', gas2.source_file);
else
    fprintf('the ''source_file'' fields were not merged because it is not present in both inputs\n')
end
gasout.occultation = [gas1.occultation, gas2.occultation];
gasout.sr1ss0 = [gas1.sr1ss0, gas2.sr1ss0];
if isfield(gas1,'beta_angle') && isfield(gas2,'beta_angle')  % when there is glc data included
    gasout.beta_angle = [gas1.beta_angle, gas2.beta_angle];
else
    fprintf('the ''beta_angle'' fields were not merged because it is not present in both inputs\n')
end
gasout.date_mjd = [gas1.date_mjd, gas2.date_mjd];
gasout.gas = gas1.gas;
if isequal(gas1.altitude_km, gas2.altitude_km) % if the altitude arrays are equal
    if size(gas1.altitude_km, 2) == 1 && size(gas1.altitude_km, 2) == 1 % if they are both vectors
        gasout.altitude_km = gas1.altitude_km;
    else % if they are arrays
        gasout.altitude_km = [gas1.altitude_km, gas2.altitude_km];
    end
else % if the arrays are not equal
    if size(gas1.altitude_km, 2) == 1 && size(gas1.altitude_km, 2) == 1 % if they are both vectors
        gasout.altitude_km = [repmat(gas1.altitude_km,[1,lgas1]), repmat(gas2.altitude_km,[1,lgas2])];
    else % if they are arrays
        gasout.altitude_km = [gas1.altitude_km, gas2.altitude_km];
    end
end
%% the vmr field can have more than 2 dimensions in special cases.
% size(gas1.vmr)
% size(gas2.vmr)
gasout.vmr = cat(length(size(gas1.vmr)), gas1.vmr, gas2.vmr); % concatenate the vmr array ilong the last dimension, which should be time
% gasout.vmr = [gas1.vmr, gas2.vmr];
%%
if isfield(gas1,'vmr_error') && isfield(gas2,'vmr_error')  % when there is glc data included
    gasout.vmr_error = [gas1.vmr_error, gas2.vmr_error];
else
    fprintf('the ''vmr_error'' fields were not merged because it is not present in both inputs\n')
end
gasout.lat_tangent = [gas1.lat_tangent, gas2.lat_tangent];
gasout.lon_tangent = [gas1.lon_tangent, gas2.lon_tangent];
if isfield(gas1,'quality_flags') && isfield(gas2,'quality_flags')
    gasout.quality_flags = [gas1.quality_flags, gas2.quality_flags];
else
    fprintf('the ''quality_flags'' field was not merged because it is not present in both inputs\n')
end
if isfield(gas1,'pressure_hPa') && isfield(gas2,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure
    if isequal(gas1.pressure_hPa, gas2.pressure_hPa) % if the altitude arrays are equal
        if size(gas1.pressure_hPa, 2) == 1 && size(gas1.pressure_hPa, 2) == 1 % if they are both vectors
            gasout.pressure_hPa = gas1.pressure_hPa;
        else % if they are arrays
            gasout.pressure_hPa = [gas1.pressure_hPa, gas2.pressure_hPa];
        end
    else % if the arrays are not equal
        if size(gas1.pressure_hPa, 2) == 1 && size(gas1.pressure_hPa, 2) == 1 % if they are both vectors
            gasout.pressure_hPa = [repmat(gas1.pressure_hPa,[1,lgas1]), repmat(gas2.pressure_hPa,[1,lgas2])];
        else % if they are arrays
            gasout.pressure_hPa = [gas1.pressure_hPa, gas2.pressure_hPa];
        end
    end
else
    warning('one or both of the inputs don''t contain the ''pressure_hPa'' field')
end
if isfield(gas1,'lon') && isfield(gas2,'lon')  % when there is glc data included
    gasout.lon = [gas1.lon, gas2.lon];
    gasout.lat = [gas1.lat, gas2.lat];
else
    fprintf('the ''lat'' and ''lon'' fields were not respectively merged because they are not present in both inputs\n')
end
if isfield(gas1,'eql') && isfield(gas2,'eql') % when there is DMP data included in the tanstruct
    gasout.eql = [gas1.eql, gas2.eql];
else
    fprintf('the ''eql'' fields were not merged because it is not present in both inputs\n')
end
if isfield(gas1,'spv') && isfield(gas2,'spv') % when there is DMP data included in the tanstruct
    gasout.spv = [gas1.spv, gas2.spv];
else
    fprintf('the ''spv'' fields were not merged because it is not present in both inputs\n')
end
if isfield(gas1,'distance') && isfield(gas2,'distance')
    gasout.distance = [gas1.distance, gas2.distance];
else
    fprintf('the ''distance'' fields were not merged because it is not present in both inputs\n')
end
if isfield(gas1,'time_diff') && isfield(gas2,'time_diff')
    gasout.time_diff = [gas1.time_diff, gas2.time_diff];
else
    fprintf('the ''time_diff'' fields were not merged because it is not present in both inputs\n')
end
if isfield(gas1,'lst_ratio') && isfield(gas2,'lst_ratio')
    gasout.lst_ratio = [gas1.lst_ratio, gas2.lst_ratio];
else
    fprintf('the ''lst_ratio'' fields were not merged because it is not present in both inputs\n');
end
if isfield(gas1,'lst') && isfield(gas2,'lst')
    gasout.lst = [gas1.lst, gas2.lst];
else
    fprintf('the ''lst'' fields were not merged because it is not present in both inputs\n');
end

%% remove duplicate occultations
disp('removing duplicate occultations...')
% orbitnames = get_ace_occultation_names(gasout);
%the next line also automatically sorts the data
[~,igood] = unique(gasout.date_mjd(1,:)); % included an indexing here because the sampled MLS data has an 'alt x occ' array of times 
% size(igood)
gasout = reduce_tanstruct_by_rowindex(gasout, igood);
disp(gasout)

%%
tanstruct_out = gasout;
%
end

