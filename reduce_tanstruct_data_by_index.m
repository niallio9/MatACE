function [ tanstruct_out ] = reduce_tanstruct_data_by_index( tanstruct_in, indices )
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
%   NJR - 02/2018

gasin = tanstruct_in;
gasout = gasin;
ygas = indices;
if isvector(ygas)
    datasize = size(gasout.vmr); 
    gasout.vmr = nan(datasize);
    gasout.vmr(ygas) = gasin.vmr(ygas);
    gasout.vmr_error = nan(datasize);
    gasout.vmr_error(ygas) = gasin.vmr_error(ygas);
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
    % now remove any colums that are all NaNs
    goodcol = find(nansum(gasout.vmr) ~= 0); % find the columns that are all nans
    gasout = reduce_tanstruct_by_rowindex(gasout,goodcol);
    
    tanstruct_out = gasout;
else
    error('The input indices should be in a vector. We are using linear indexing here')
end

end

