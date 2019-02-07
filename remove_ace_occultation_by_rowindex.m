function [ tanstruct_out ] = remove_ace_occultation_by_rowindex( tanstruct_in, rowindices )
%A function to reduce the ace data according to the provided
%indicies. The indices indicate the occultations that you would like to
%delete.

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
%   NJR - 02/2019 

gasout = tanstruct_in;
to_remove = rowindices;
gasout.occultation(:, to_remove) = [];
gasout.sr1ss0(:, to_remove) = [];
gasout.beta_angle(:, to_remove) = [];
gasout.date_mjd(:, to_remove) = [];
gasout.vmr(:, to_remove) = [];
gasout.vmr_error(:, to_remove) = [];
gasout.lat_tangent(:, to_remove) = [];
gasout.lon_tangent(:, to_remove) = [];
% if isfield (gasout,'source_file')
%         gasout.source_file = gasout.source_file(:,ygas);
% end
if isfield (gasout,'version')
    gasout.version(:, to_remove) = [];
end
if isfield(gasout,'quality_flags')
    gasout.quality_flags(:, to_remove) = [];
    %     else
    %         fprintf('\njust to let you know, there''s no quality flags in here\n')
end
if length(gasout.altitude_km(1,:)) > 1 % for the data that has been interpolated to a pressure grid
    gasout.altitude_km(:, to_remove) = [];
end
if isfield(gasout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
    if length(gasout.pressure_hPa(1,:)) > 1
        gasout.pressure_hPa(:, to_remove) = [];
    end
end
if isfield(gasout,'lon') % when there is glc data included
    gasout.lon(:, to_remove) = [];
    gasout.lat(:, to_remove) = [];
end
if isfield(gasout,'eql') % when there is DMP data included in the tanstruct
    gasout.eql(:, to_remove) = [];
end
if isfield(gasout,'spv') % when there is DMP data included in the tanstruct
    gasout.spv(:, to_remove) = [];
end
if isfield(gasout,'theta') % when there is DMP data included in the tanstruct
    gasout.theta(:, to_remove) = [];
end
if isfield(gasout,'distance')
    gasout.distance(:, to_remove) = [];
end
if isfield(gasout,'time_diff')
    gasout.time_diff(:, to_remove) = [];
end
if isfield(gasout,'lst_ratio')
    gasout.lst_ratio(:, to_remove) = [];
end

tanstruct_out = gasout;

if isempty(gasout.occultation)
   warning('The ace structure has been reduced to zero entries')
end
%
end

