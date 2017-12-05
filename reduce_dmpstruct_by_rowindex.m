function [ dmpstruct_out ] = reduce_dmpstruct_by_rowindex( dmpstruct_in, rowindices )
%A function to reduce the ace dmp data according to the
%provided indicies.

% *INPUT*
%           dmpstruct_in: STRUCTURE - containins the ACE dmp data. It is
%           usually created using 'read_ace_dmp' or 'read_ace_dmp_for_mat'.
%
%           rowindices: VECTOR - the indices of the occultations that
%           you would like to keep.
%
% *OUTPUT*
%           dmpstruct_out: STRUCTURE with the same fields as the
%           input but with a size reduced according to the input indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 10/17

dmpout = dmpstruct_in;
ydmp = rowindices;

dmpout.source_file = dmpout.source_file(ydmp);
dmpout.occultation = dmpout.occultation(ydmp);
dmpout.sr1ss0 = dmpout.sr1ss0(ydmp);
dmpout.date_mjd = dmpout.date_mjd(ydmp);
dmpout.lat_tangent = dmpout.lat_tangent(ydmp);
dmpout.lon_tangent = dmpout.lon_tangent(ydmp);
dmpout.version = dmpout.version(ydmp);
dmpout.T = dmpout.T(:,ydmp);
if length(dmpout.pressure_hPa(1,:)) > 1
    dmpout.pressure_hPa = dmpout.pressure_hPa(:,ydmp);
end
dmpout.lon = dmpout.lon(:,ydmp);
dmpout.lat = dmpout.lat(:,ydmp);
dmpout.Theta = dmpout.Theta(:,ydmp);
dmpout.spv = dmpout.spv(:,ydmp);
dmpout.eql = dmpout.eql(:,ydmp);
if length(dmpout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    dmpout.altitude_km = dmpout.altitude_km(:,ydmp);
end

dmpstruct_out = dmpout;

if isempty(ydmp)
   warning('The ace dmp structure has been reduced to zero entries')
end
%
end

