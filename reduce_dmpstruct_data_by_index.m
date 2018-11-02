function [ dmpstruct_out ] = reduce_dmpstruct_data_by_index( dmpstruct_in, indices, delete_rows )
%A function to reduce the ace dmp data according to the provided
%indices. Data that do not correspond to the indices are changed to NaNs.
%If a whole row (occultation) of data is converted to NaNs, then that 
%occultation is removed from the structure.

% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'. 
%
%           indices: VECTOR - the indices of the data that you would
%           like to keep. Linear indexing is used (see matlab doc).
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE - output has almost the same fields as
%           the input, but with the data reduced according to the input
%           indices.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/17
%   NJR - 10/18 - edited for new DMP format

dmpin = dmpstruct_in;
dmpout = dmpin;
ydmp = indices;
if isvector(ydmp)
    datasize = size(dmpout.T); 
    dmpout.T = nan(datasize);
    dmpout.T(ydmp) = dmpin.T(ydmp);
    dmpout.lon = nan(datasize);
    dmpout.lon(ydmp) = dmpin.lon(ydmp);
    dmpout.lat = nan(datasize);
    dmpout.lat(ydmp) = dmpin.lat(ydmp);
    dmpout.Theta = nan(datasize);
    dmpout.Theta(ydmp) = dmpin.Theta(ydmp);
    dmpout.spv = nan(datasize);
    dmpout.spv(ydmp) = dmpin.spv(ydmp);
    dmpout.eql = nan(datasize);
    dmpout.eql(ydmp) = dmpin.eql(ydmp);
    if length(dmpout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
        dmpout.altitude_km = nan(datasize);
        dmpout.altitude_km(ydmp) = dmpin.altitude_km(ydmp);
    end
    if length(dmpout.pressure_hPa(1,:)) > 1
        dmpout.pressure_hPa = nan(datasize);
        dmpout.pressure_hPa(ydmp) = dmpin.pressure_hPa(ydmp);
    end
    
    if nargin < 3 || delete_rows == 1
        % now remove any colums that are all NaNs
        goodcol = find(nansum(dmpout.lon) ~= 0); % find the columns that are all nans
        dmpout = reduce_dmpstruct_by_rowindex(dmpout,goodcol);
    end
    
    dmpstruct_out = dmpout;
    
else
    error('The input indices should be in a vector. We are using linear indexing here')
end

end

