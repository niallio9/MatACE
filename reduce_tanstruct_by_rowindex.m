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
%   NJR - 07/2018 can use for structures with 'distance', 'time_diff',
%   and/or 'lst_ratio' fields.
%   NJR - 11/2018 can use with source_file and version: included in MAESTRO
%   structures, and others ***** THIS DOESNT REALLY WORK. USE ONE SOURCE FILE...************FIX************************* 

gasout = tanstruct_in;
ygas = rowindices;
gasout.occultation = gasout.occultation(:,ygas);
gasout.sr1ss0 = gasout.sr1ss0(:,ygas);
if isfield (gasout,'beta_angle')
    gasout.beta_angle = gasout.beta_angle(:,ygas);
end
gasout.date_mjd = gasout.date_mjd(:,ygas);
%% the vmr field can have more than 2 dimensions in special cases.
if length(size(gasout.vmr)) == 2
    gasout.vmr = gasout.vmr(:,ygas);
elseif length(size(gasout.vmr)) == 3
    gasout.vmr = gasout.vmr(:,:,ygas);
else
    error('The scripts are set up to allow vmr arrays of 2 or 3 dimensions. There are %i dimensions here', length(size(gasout.vmr)))
end
%%
if isfield (gasout,'vmr_error')
    gasout.vmr_error = gasout.vmr_error(:,ygas);
end
gasout.lat_tangent = gasout.lat_tangent(:,ygas);
gasout.lon_tangent = gasout.lon_tangent(:,ygas);
% if isfield (gasout,'source_file')
%         gasout.source_file = gasout.source_file(:,ygas);
% end
if isfield (gasout,'version')
    gasout.version = gasout.version(:,ygas);
end
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
if isfield(gasout,'spv') % when there is DMP data included in the tanstruct
    gasout.spv = gasout.spv(:,ygas);
end
if isfield(gasout,'theta') % when there is DMP data included in the tanstruct
    gasout.theta = gasout.theta(:,ygas);
end
if isfield(gasout,'distance')
    gasout.distance = gasout.distance(:,ygas);
end
if isfield(gasout,'time_diff')
    gasout.time_diff = gasout.time_diff(:,ygas);
end
if isfield(gasout,'lst_ratio')
    gasout.lst_ratio = gasout.lst_ratio(:,ygas);
end
if isfield(gasout,'lst')
    gasout.lst = gasout.lst(:,ygas);
end

tanstruct_out = gasout;

if isempty(ygas)
   warning('The ace structure has been reduced to zero entries')
end
%
end

