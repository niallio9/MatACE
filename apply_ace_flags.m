function [ tanstruct_out_flagged ] = apply_ace_flags(tanstruct_in)
%A function to apply the flags that are supplied for the ACE data. The
%profiles that are flagged as 'not to be used' are removed. Data points
%that are not be used are or replaced with NaNs. 
% The ACE flags are descibed in Sheese et al., 2015: Detecting physically
% unrealistic outliers in ACE-FTS atmospheric measurements.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out_flagged: STRUCTURE - output has the same
%           fields as the input, but with the flagged data removed or
%           edited.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017

%% Remove flagged values
datain = tanstruct_in;
if ~isfield(datain,'quality_flags')
    error('There is no quality flag information here. Has this data been interoplated already?')    
end
out = tanstruct_in; %this structure will be edited later on

% find where 4 <= quality flag <= 7. The whole profile is to be removed in this case.
[~, badj] = find(datain.quality_flags >= 4 & datain.quality_flags <= 7); % i -> altitudes. j -> occultations
badj = unique(badj); % remove multiple listings of the same column
out.occultation(badj) = []; % remove the values corresponding to flagged occultation(s)
out.sr1ss0(badj) = [];
out.beta_angle(badj) = [];
out.date_mjd(:,badj) = [];
out.lat_tangent(badj) = [];
out.lon_tangent(badj) = [];
out.vmr(:,badj) = [];
out.vmr_error(:,badj) = [];
out.quality_flags(:,badj) = [];
if length(out.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    out.altitude_km(:,badj) = [];
    fprintf('\nWarning: The input structure appears to have already been interpolated. Applying the flags here is a bit iffy...\n')
end
if isfield(out,'pressure_hPa') % for the data that has been interpolated to a pressure grid
   if length(out.pressure_hPa(1,:)) > 1
      out.pressure_hPa(:,badj) = [];
   end
end
if isfield(out,'lon')
   out.lon(:,badj) = [];
   out.lat(:,badj) = [];
end
if isfield(out,'eql')
   out.eql(:,badj) = []; 
end
if isfield(out,'spv')
   out.spv(:,badj) = []; 
end
if isfield(out,'date_mjd_mls')
   out.date_mjd_mls(:,badj) = []; 
end

%% Change any other values with flags > 1 to NaNs
[badI] = find(out.quality_flags > 1);
out.vmr(badI) = nan;
out.vmr_error(badI) = nan;
if isfield(out,'lon')
    out.lon(badI) = nan;
    out.lat(badI) = nan;
end
if isfield(out,'eql')
    out.eql(badI) = nan;
end
if isfield(out,'spv')
    out.spv(badI) = nan;
end
% remove any columns that are only nans now
goodcol = find(nansum(out.vmr) ~= 0);
out = reduce_tanstruct_by_rowindex(out, goodcol);

%%
tanstruct_out_flagged = out;

end

