function [ tanstruct_out ] = reduce_tanstruct_by_rowindex( tanstruct_in, rowindices )
%A function to reduce the ace data according to the provided
%indicies. The indices indicate the occultations that you would like to
%keep.

% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           rowindices: VECTOR - the indices of the occultations that
%           you would like to keep.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same fields as the
%           input, but with a size reduced according to the input indices.
%           Only the occultations that correspond to the row indices
%           remain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NJR - 11/2017
%   NJR - 02/2018 - include 'eql' field

gasout = tanstruct_in;
ygas = rowindices;
gasout.occultation = gasout.occultation(ygas);
gasout.sr1ss0 = gasout.sr1ss0(ygas);
gasout.beta_angle = gasout.beta_angle(ygas);
gasout.date_mjd = gasout.date_mjd(ygas);
gasout.vmr = gasout.vmr(:,ygas);
gasout.vmr_error = gasout.vmr_error(:,ygas);
gasout.lat_tangent = gasout.lat_tangent(ygas);
gasout.lon_tangent = gasout.lon_tangent(ygas);
if isfield(gasout,'quality_flags')
    gasout.quality_flags = gasout.quality_flags(:,ygas);
    %     else
    %         fprintf('\njust to let you know, there''s no quality flags in here\n')
end
if length(gasout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    gasout.altitude_km = gasout.altitude_km(:,ygas);
end
if isfield(gasout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
    if length(gasout.pressure_hPa(1,:)) > 1
        gasout.pressure_hPa = gasout.pressure_hPa(:,ygas);
    end
end
if isfield(gasout,'lon') % when there is glc data included
    gasout.lon = gasout.lon(:,ygas);
    gasout.lat = gasout.lat(:,ygas);
end
if isfield(gasout,'eql') % when there is DMP data included in the tanstruct
    gasout.eql = gasout.eql(:,ygas);
end

tanstruct_out = gasout;

if isempty(ygas)
   warning('The ace structure has been reduced to zero entries')
end
%
end

