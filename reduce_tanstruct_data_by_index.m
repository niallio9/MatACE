function [ tanstruct_out ] = reduce_tanstruct_data_by_index( tanstruct_in, indices, delete_rows )
%A function to reduce the ace data according to the provided
%indicies. Data that do not correspond to the indices are changed to NaNs.
%If a whole row (occultation) of data is converted to NaNs, then that 
%occultation is removed from the structure.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'. 
%
%           indices: VECTOR - the indices of the data that you would like
%           to keep. Linear indexing is used (see matlab doc). 
%
% *OUTPUT*
%           tanstruct_out: MATLAB STRUCTURE with the same fields as the
%           input but with the data reduced according to the input indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/2017
%   NJR - 02/2018 include 'eql' field
%   NJR - 05/2018 can use for tanstructs that have more than one row of
%   date or lat_tangent / lon_tangent info
%   NJR - 07/2018 can use for structures with 'distance', 'time_diff',
%   and/or 'lst_ratio' fields.

gasin = tanstruct_in;
gasout = gasin;
ygas = indices;
if isvector(ygas)
    datasize = size(gasout.vmr); 
    gasout.vmr = nan(datasize);
    gasout.vmr(ygas) = gasin.vmr(ygas);
    gasout.vmr_error = nan(datasize);
    gasout.vmr_error(ygas) = gasin.vmr_error(ygas);
    if isequal(size(gasin.date_mjd), datasize) % for a case where the tanstruct has been modified to have a different info for each data point
        gasout.date_mjd = nan(datasize);
        gasout.date_mjd(ygas) = gasin.date_mjd(ygas);
    end
    if isequal(size(gasin.lat_tangent), datasize) % for a case where the tanstruct has been modified to have a different info for each data point
        gasout.lat_tangent = nan(datasize);
        gasout.lat_tangent(ygas) = gasin.lat_tangent(ygas);
        gasout.lon_tangent = nan(datasize);
        gasout.lon_tangent(ygas) = gasin.lon_tangent(ygas);
    end
    if isfield(gasout,'quality_flags')
        gasout.quality_flags = nan(datasize);
        gasout.quality_flags(ygas) = gasin.quality_flags(ygas);
%     else
%         fprintf('\njust to let you know, there''s no quality flags being indexed in here\n')
    end
    if length(gasout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
        gasout.altitude_km = nan(datasize);
        gasout.altitude_km(ygas) = gasin.altitude_km(ygas);
    end
    if isfield(gasout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
        if length(gasout.pressure_hPa(1,:)) > 1
            gasout.pressure_hPa = nan(datasize);
            gasout.pressure_hPa(ygas) = gasin.pressure_hPa(ygas);
        end
    end
    if isfield(gasout,'lon') % when there is glc data included in tanstruct
        gasout.lon = nan(datasize);
        gasout.lat = nan(datasize);
        gasout.lon(ygas) = gasin.lon(ygas);
        gasout.lat(ygas) = gasin.lat(ygas);
    end
    if isfield(gasout,'eql') % when there is DMP data included in the tanstruct
        gasout.eql = nan(datasize);
        gasout.eql(ygas) = gasin.eql(ygas);
    end
    if isfield(gasout,'spv') % when there is DMP data included in the tanstruct
        gasout.spv = nan(datasize);
        gasout.spv(ygas) = gasin.spv(ygas);
    end
    if isfield(gasout,'theta') % when there is DMP data included in the tanstruct
        gasout.theta = nan(datasize);
        gasout.theta(ygas) = gasin.theta(ygas);
    end
    if isfield(gasout,'distance') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.distance = nan(datasize);
        gasout.distance(ygas) = gasin.distance(ygas);
    end
    if isfield(gasout,'time_diff') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.time_diff = nan(datasize);
        gasout.time_diff(ygas) = gasin.time_diff(ygas);
    end
    if isfield(gasout,'lst_ratio') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.lst_ratio = nan(datasize);
        gasout.lst_ratio(ygas) = gasin.lst_ratio(ygas);
    end
    if isfield(gasout,'lst') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.lst = nan(datasize);
        gasout.lst(ygas) = gasin.lst(ygas);
    end
    
    if nargin < 3 || delete_rows == 1
        % now remove any colums that are all NaNs
        goodcol = find(nansum(gasout.vmr) ~= 0); % find the columns that are all nans
        gasout = reduce_tanstruct_by_rowindex(gasout,goodcol);
    end
    
    tanstruct_out = gasout;
else
    error('The input indices should be in a vector. We are using linear indexing here')
end

end

