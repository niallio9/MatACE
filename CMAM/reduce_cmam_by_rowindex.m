function [ cmam_out ] = reduce_cmam_by_rowindex( cmam_in, rowindices )
%A function to reduce the ace data according to the provided
%indicies. The indices indicate the occultations that you would like to
%keep.

% *INPUT*
%           cmam_in: STRUCTURE - contains the gas specific CMAM data.
%           This structure can be created with 'read_cmam_ncdata.m'. The
%           structure created by 'sample_cmam_for_ace.m' can also be used
%           as input here.
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
%   NJR - 11/17

gasout = cmam_in;
ygas = rowindices;
if isfield(gasout,'occultation') %this should be present for the cmam data that is sampled to ACE locations
    gasout.occultation = gasout.occultation(ygas);
    gasout.sr1ss0 = gasout.sr1ss0(ygas);
end
gasout.date_mjd = gasout.date_mjd(ygas);
gasout.vmr = gasout.vmr(:,ygas);
if isfield(gasout,'pressure_hPa') % for the data that has been interpolated to a pressure grid, or the new version that includes pressure included
    if length(gasout.pressure_hPa(1,:)) > 1
        gasout.pressure_hPa = gasout.pressure_hPa(:,ygas);
    end
end
if isfield(gasout,'lon') % this should always be true for the inputs
    gasout.lon = gasout.lon(:,ygas);
    gasout.lat = gasout.lat(:,ygas);
end

cmam_out = gasout;

if isempty(ygas)
   warning('The cmam structure has been reduced to zero entries')
end
%
end

