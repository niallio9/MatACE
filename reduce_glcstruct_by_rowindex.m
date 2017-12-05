function [ tanstruct_out ] = reduce_glcstruct_by_rowindex( glcstruct_in, rowindices )
%A function to reduce the ace glc data according to the provided
%indicies. The indices indicate the data that you would like to keep.

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
%   NJR - 10/17

glcout = glcstruct_in;
yglc = rowindices;

glcout.occultation = glcout.occultation(yglc);
glcout.sr1ss0 = glcout.sr1ss0(yglc);
glcout.beta_angle = glcout.beta_angle(yglc);
glcout.date_mjd = glcout.date_mjd(yglc);
glcout.lat_tangent = glcout.lat_tangent(yglc);
glcout.lon_tangent = glcout.lon_tangent(yglc);
glcout.lon = glcout.vmr(:,yglc);
glcout.lat = glcout.vmr_error(:,yglc);
if length(glcout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    glcout.altitude_km = glcout.altitude_km(:,yglc);
end
glcout.temperature = glcout.temperature(yglc);
if isfield(glcout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
    if length(glcout.pressure_hPa(1,:)) > 1
        glcout.pressure_hPa = glcout.pressure_hPa(:,yglc);
    end
end
glcout.lon = glcout.lon(:,yglc);
glcout.lat = glcout.lat(:,yglc);
glcout.sun_heading = glcout.sun_heading(:,yglc);

tanstruct_out = glcout;

if isempty(yglc)
   warning('The ace structure has been reduced to zero entries')
end
%
end

