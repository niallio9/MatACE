function [tanstruct_inside, tanstruct_outside, tanstruct_edge, dmpstruct_inside, dmpstruct_outside, dmpstruct_edge ] = ...
        subset_ace_by_vortex( tanstruct_in, dmpstruct_in )
%A function to subset ace data corresponding to whether the location of a
%measurement point with respect to the polar vortex: inside, outside, or
%within the edge of the vortex

% *INPUT*
%           tanstruct_in: STRUCTURE - a .MAT structure containing ACE
%           data and metadata. It is usually created using
%           'read_ace_ncdata', 'read_ace_ncdata_for_mat'.
%
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           dmp info.
%
%           dmpstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the data that coincides with the
%           gas info.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Define some things
gas = tanstruct_in;
dmp = dmpstruct_in;
dmp = remove999_ace_dmp(dmp);
[gas, dmp] = match_ace_data_dmp(gas, dmp);
spv_in = 1.6e-4;
spv_out = 1.2e-4;


%% subset the ace data to match the altitude range of the dmp data
disp('subsetting the ACE data to the altitude range of the DMP data...')
gas = subset_ace_by_alt(gas, dmp.altitude_km(1), dmp.altitude_km(end));
izlow = find(gas.altitude_km == dmp.altitude_km(1)); % where the gas data altitude matches the min alt' of the dmp 
izhigh = find(gas.altitude_km == dmp.altitude_km(end)); % where the gas data altitude matches the max alt' of the dmp
if isempty(izlow) || isempty(izhigh)
    error('the altitude fields of the measurement and DMP data don''t match. Stopping...')
end
% make sure that the altitudes match 
if ~isequal(gas.altitude_km(izlow:izhigh,:), dmp.altitude_km)
    error('the altitude fields of the subsetted measurement and DMP data don''t match. Stopping...')
end
disp('done')

% fill the spv data into an array that matches the size of the gas data.
% Will use this for indexing later.
[gas, dmp] = match_ace_data_dmp(gas, dmp);
spv = dmp.spv;
spv_acegrid = nan(size(gas.vmr)); % will want to use this for creating an index later
% spvsign = sign(spv);
% spv = abs(spv);
spv_acegrid(izlow:izhigh,:) = spv;
% spv_acegridsign = sign(spv_acegrid);


%% pick out the data that corresponds to the vortex relative location
% reduce and subset the DMP data
ispv_in = find(abs(spv) > spv_in); % get the indices of the data that lie in the vortex
ispv_out = find(abs(spv) < spv_out); % get the indices of the latitudes that lie within the chosen range
ispv_edge = find(abs(spv) >= spv_out & abs(spv) <= spv_in); % get the indices of the latitudes that lie within the chosen range
% subset the dmp data
dmpout_in = reduce_dmpstruct_data_by_index(dmp,ispv_in);
% dmpout_in.spv = spvsign.*dmpout_in.spv;
dmpout_out = reduce_dmpstruct_data_by_index(dmp,ispv_out);
dmpout_edge = reduce_dmpstruct_data_by_index(dmp,ispv_edge);

% same for the gas data
ispv_in = find(abs(spv_acegrid) > spv_in); % get the indices of the data that lie in the vortex
ispv_out = find(abs(spv_acegrid) < spv_out); % get the indices of the latitudes that lie within the chosen range
ispv_edge = find(abs(spv_acegrid) >= spv_out & abs(spv_acegrid) <= spv_in); % get the indices of the latitudes that lie within the chosen range
%Subset the data
gasout_in = reduce_tanstruct_data_by_index(gas,ispv_in);
gasout_out = reduce_tanstruct_data_by_index(gas,ispv_out);
gasout_edge = reduce_tanstruct_data_by_index(gas,ispv_edge);

%% match the data again so that the rows are the same size
% some of the rows that didn't have any data left may have been deleted
disp('subsetting data by vortex-relative position...')
[~,gasout_in, dmpout_in] = evalc('match_ace_data_dmp(gasout_in, dmpout_in);');
[~,gasout_out, dmpout_out] = evalc('match_ace_data_dmp(gasout_out, dmpout_out);');
[~,gasout_edge, dmpout_edge] = evalc('match_ace_data_dmp(gasout_edge, dmpout_edge);');
% ass the spv field to the measurement data
gasout_in.spv = nan(size(gasout_in.vmr));
gasout_in.spv(izlow:izhigh,:) = dmpout_in.spv;
gasout_out.spv = nan(size(gasout_out.vmr));
gasout_out.spv(izlow:izhigh,:) = dmpout_out.spv;
gasout_edge.spv = nan(size(gasout_edge.vmr));
gasout_edge.spv(izlow:izhigh,:) = dmpout_edge.spv;
disp('Done :)')
%
tanstruct_inside = gasout_in;
tanstruct_outside = gasout_out;
tanstruct_edge = gasout_edge;

dmpstruct_inside = dmpout_in;
dmpstruct_outside = dmpout_out;
dmpstruct_edge = dmpout_edge;
%
end

