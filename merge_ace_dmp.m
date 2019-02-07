function [ tanstruct_out ] = merge_ace_dmp( tanstruct_in, dmpstruct_in)
%A function to read ACE DMP data, according to the occultations of
%the input, and add the latitude and longitutde information to the input
%structure. The input ACE data should use the original ACE 1km-grid. The
%tangent latitude and longitude are used when there is no information in
%the DMP file.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           dmpstruct_in: Structure - contains the DMP data for the
%           ace measurements.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with ACE DMP information added: sp, eql, and theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 12/2018
tic
%% Only do if the lat and lon info isnt already there
gas = tanstruct_in;
dmp = dmpstruct_in;
dmp = remove999_ace_dmp(dmp);

tanstruct_out = gas;
tanstruct_out.spv = nan(size(gas.vmr));
tanstruct_out.eql = nan(size(gas.vmr));
tanstruct_out.theta = nan(size(gas.vmr));

izlow = find(gas.altitude_km == dmp.altitude_km(1)); % where the gas data altitude matches the min alt' of the dmp 
izhigh = find(gas.altitude_km == dmp.altitude_km(end)); % where the gas data altitude matches the max alt' of the dmp
if isempty(izlow) || isempty(izhigh)
    error('the altitude fields of the measurement and DMP data don''t match. Stopping...')
end
% make sure that the altitudes match 
if ~isequal(gas.altitude_km(izlow:izhigh,:), dmp.altitude_km)
    error('the altitude fields of the subsetted measurement and DMP data don''t match. Stopping...')
end

[ ~, ~, ygas, ydmp ] = match_ace_data_dmp( gas, dmp );

if isempty(ygas)
    disp('there are no matching data found')
    tanstruct_out = [];
    return
end

tanstruct_out.spv(izlow:izhigh, ygas) = dmp.spv(:,ydmp);
tanstruct_out.eql(izlow:izhigh, ygas) = dmp.eql(:,ydmp);
tanstruct_out.theta(izlow:izhigh, ygas) = dmp.theta(:,ydmp);

%
toc
end

