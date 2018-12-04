function [ tanstruct_maestro_out ] = merge_maestro_ace_pressure( tanstruct_maestro, tanstruct_ace)
%A function to read ACE GLC data files, according to the occultations of
%the input, and add the latitude and longitutde information to the input
%structure. The input ACE data should use the original ACE 1km-grid. The
%tangent latitude and longitude are used when there is no information in
%the GLC file.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           glcstruct_in: STRUCTURE - contains the location data for the
%           maestro measurements. This is also refered to a a SunriseTable
%           or SunsetTable or SunriseSunsetTable or sr1ss0_table
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with MAESTRO GLC information added: latitude and
%           longitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 05/2018
tic
%% define some things
interptype = 'pchip';
maestro = tanstruct_maestro;
ace = tanstruct_ace;

[maestro, ace] = match_ace_data(maestro, ace);
% maestro
logpace = log(ace.pressure_hPa);
% interpolate the ace pressures to the maestro altitudes

for i = 1:length(maestro.occultation)
    logpmaestroi = interp1(ace.altitude_km(:,1), logpace(:,i), maestro.altitude_km(:,i), interptype, nan);
    maestro.pressure_hPa(:,i) = exp(logpmaestroi);
    % maestro.altitude_km = maestro.altitude_km(:,1);
end

%% deal with NaN pressures created from the interpolation
[badI] = find(isnan(maestro.pressure_hPa));
maestro.vmr(badI) = nan;
maestro.vmr_error(badI) = nan;
if length(maestro.altitude_km(1,:)) > 1
    maestro.altitude_km(badI) = nan;
end
if isfield(maestro,'lon')
   maestro.lon(badI) = nan;
   maestro.lat(badI) = nan;
end

%% remove any columns that are only nans now
goodcol = find(nansum(maestro.vmr) ~= 0);
maestro = reduce_tanstruct_by_rowindex(maestro, goodcol);
goodcol = find(nansum(maestro.vmr_error) ~= 0); % same for the errors columns that are all nans
maestro = reduce_tanstruct_by_rowindex(maestro, goodcol);

%% filter out any data for which the vmr error is larger than 50% of the vmr
if ~strcmp(maestro.gas,'H2O') % don't do this or the H2O data because the error is massive for that product
    vmrdif = (abs(maestro.vmr) - abs(maestro.vmr_error))./abs(maestro.vmr);
    goodcol = find(vmrdif > 0.5 & maestro.vmr_error > 0 & maestro.vmr < 20e-6);
    % whos
    maestro = reduce_tanstruct_data_by_index(maestro, goodcol);
end

%% for the version 1.2 data, SPARC book says that the altitude range is 5-60km.
badrow = find(maestro.altitude_km(:,1) > 60);
maestro.vmr(badrow,:) = nan;
maestro.vmr_error(badrow,:) = nan;
if length(maestro.altitude_km(1,:)) > 1
    maestro.altitude_km(badrow,:) = nan;
end
if isfield(maestro,'lon')
   maestro.lon(badrow,:) = nan;
   maestro.lat(badrow,:) = nan;
end


tanstruct_maestro_out = maestro;
%
toc
end

