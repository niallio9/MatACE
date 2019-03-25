function [ tanstruct_out ] = remove_ace_data_by_index( tanstruct_in, indices, delete_rows )
%A function to reduce the ace data according to the provided
%indices. Data that do not correspond to the indices are changed to NaNs.
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
%   NJR - 02/2019

gasin = tanstruct_in;
gasout = gasin;
to_remove = indices;
if isvector(to_remove)
    datasize = size(gasout.vmr); 
    gasout.vmr(to_remove) = nan;
    gasout.vmr_error(to_remove) = nan;
    if isequal(size(gasin.date_mjd), datasize) % for a case where the tanstruct has been modified to have a different info for each data point
        gasout.date_mjd(to_remove) = nan;
    end
    if isequal(size(gasin.lat_tangent), datasize) % for a case where the tanstruct has been modified to have a different info for each data point
        gasout.lat_tangent(to_remove) = nan;
        gasout.lon_tangent(to_remove) = nan;
    end
    if isfield(gasout,'quality_flags')
        gasout.quality_flags(to_remove) = nan;
%     else
%         fprintf('\njust to let you know, there''s no quality flags being indexed in here\n')
    end
    if length(gasout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
        gasout.altitude_km(to_remove) = nan;
    end
    if isfield(gasout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
        if length(gasout.pressure_hPa(1,:)) > 1
            gasout.pressure_hPa(to_remove) = nan;
        end
    end
    if isfield(gasout,'lon') % when there is glc data included in tanstruct
        gasout.lon(to_remove) = nan;
        gasout.lat(to_remove) = nan;
    end
    if isfield(gasout,'eql') % when there is DMP data included in the tanstruct
        gasout.eql(to_remove) = nan;
    end
    if isfield(gasout,'spv') % when there is DMP data included in the tanstruct
        gasout.spv(to_remove) = nan;
    end
    if isfield(gasout,'theta') % when there is DMP data included in the tanstruct
        gasout.theta(to_remove) = nan;
    end
    if isfield(gasout,'distance') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.distance(to_remove) = nan;
    end
    if isfield(gasout,'time_diff') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.time_diff(to_remove) = nan;
    end
    if isfield(gasout,'lst_ratio') % for structures that have been created from other datasets using coincidence criteria, etc.
        gasout.lst_ratio(to_remove) = nan;
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

