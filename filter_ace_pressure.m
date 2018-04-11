function [ data_structure_flagged ] = filter_ace_pressure(tanstruct_in)
%A function to find data points in an ACE pressure measurement that are
%zeros, and remove those measurements.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           data_structure_flagged: STRUCTURE - output has the same
%           fields as the input, but with occultations removed that have
%           zero pressure values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 10/2017
% NJR - 04/2018 - use 'reduce-tanstruct_by_rowindex.m'

%% Remove measurements with zero values
datain = tanstruct_in;
[~, goodj] = find(sum(datain.pressure_hPa == 0) == 0); % i -> altitudes. j -> occultations
goodj = unique(goodj); % remove multiple listings of the same column
out = reduce_tanstruct_by_rowindex(datain, goodj);

% datain = tanstruct_in;
% out = tanstruct_in; %this structure will be edited later on
% 
% [~, badj] = find(datain.pressure_hPa == 0); % i -> altitudes. j -> occultations
% badj = unique(badj); % remove multiple listings of the same column
% out.occultation(badj) = []; % remove the values corresponding to flagged occultation(s)
% out.sr1ss0(badj) = [];
% out.beta_angle(badj) = [];
% out.date_mjd(badj) = [];
% out.lat_tangent(badj) = [];
% out.lon_tangent(badj) = [];
% out.vmr(:,badj) = [];
% out.vmr_error(:,badj) = [];
% if isfield(out,'quality_flags')
%     out.quality_flags(:,badj) = [];
% end
% if length(out.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
%     out.altitude_km(:,badj) = [];
%     fprintf('\nWarning: The input structure appears to have already been interpolated. Applying the filter here is a bit iffy...\n')
% end
% if isfield(out,'pressure_hPa') % for the data that has been interpolated to a pressure grid
%    if length(out.pressure_hPa(1,:)) > 1
%       out.pressure_hPa(:,badj) = [];
%    end
% end
% if isfield(out,'lon')
%    out.lon(:,badj) = [];
%    out.lat(:,badj) = [];
% end
% if isfield(gas,'eql')
%    out.eql(:,badj) = []; 
% end

%%
data_structure_flagged = out;
%
end

