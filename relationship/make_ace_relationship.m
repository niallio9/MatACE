function [ relstruct_out ] = make_ace_relationship( tanstruct_x_in, tanstruct_y_in )
%A function to calculate the relationship between two ACE gas measurements.
%The correlation and linear relationship of coincident data points is
%calculated.

% *INPUT*
%           tanstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
%           tanstruct1_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has the same
%           fields as the input, but with the coincident VMR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 09/2018

%% Define some things
gasx = apply_ace_flags(tanstruct_x_in);
gasy = apply_ace_flags(tanstruct_y_in);
min_alt = 0;
max_alt = 35;

if ~isequal(gasx.altitude_km, gasy.altitude_km)
    error('the altitude grids of each gas are not equal. Stopping...')
end

%% subset the data to a certain altitude
gasx = subset_ace_by_alt(gasx, min_alt, max_alt);
gasy = subset_ace_by_alt(gasy, min_alt, max_alt);

%% match the data
[gasx, gasy] = match_ace_data(gasx,gasy);
norbit = length(gasx.occultation);
nalt = length(gasx.altitude_km(:,1));

%% reshape data for correlation/slope
column_size = [norbit.*nalt, 1];
vmrcol1 = reshape(gasx.vmr, column_size);
vmrerrorcol1 = reshape(gasx.vmr_error, column_size);
vmrcol2 = reshape(gasy.vmr, column_size);
vmrerrorcol2 = reshape(gasy.vmr_error, column_size);
% if length(gas1.altitude_km(1,:)) == 1
%     altcol = repmat(gas1.altitude_km,[norbit,1]);
% else
%     altcol = reshape(gas1.altitude_km, column_size);
% end

%% get the correlation and the slope of the gas vmrs
Inonan = find(~isnan(vmrcol1) & ~isnan(vmrcol2));
vmrcol1 = vmrcol1(Inonan);
vmrcol2 = vmrcol2(Inonan);
vmrerrorcol1 = vmrerrorcol1(Inonan);
vmrerrorcol2 = vmrerrorcol2(Inonan);
% altcol = altcol(Inonan);
% correlation matrix
% size(vmrcol1)
% size(vmrcol2)
[R, P] = corrcoef(vmrcol1, vmrcol2);
% slope and intercept
slope = nan(1,2);
intercept = nan(1,2);
[intercept(1), slope(1), intercept(2), slope(2)] = york_fit(vmrcol1',vmrcol2',vmrerrorcol1',vmrerrorcol2');
% slope
% intercept

out.source_file = strcat(gasx.source_file,', ', gasy.source_file);
out.occultation = gasx.occultation;
out.sr1ss0 = gasx.sr1ss0;
out.beta_angle = gasx.beta_angle;
out.date_mjd = gasx.date_mjd;
out.gas_x = gasx.gas;
out.gas_y = gasy.gas;
out.altitude_km = gasx.altitude_km;
out.vmr_x = gasx.vmr;
out.vmr_y = gasy.vmr;
out.vmr_error_x = gasx.vmr_error;
out.vmr_error_y = gasy.vmr_error;
out.slope = slope(1);
out.slope_error = slope(2);
out.intercept = intercept(1);
out.intercept_error = intercept(2);
if length(R(:,1)) > 1
    out.correlation = R(1,2);
    out.p_value = P(1,2);
else
    out.correlation = nan;
    out.p_value = nan;
end
out.lat_tangent = gasx.lat_tangent;
out.lon_tangent = gasx.lon_tangent;
out.quality_flags_x = gasx.quality_flags;
out.quality_flags_y = gasy.quality_flags;
out.pressure_hPa = gasx.pressure_hPa;
if isfield(gasx,'lon')
    out.lon = gasx.lon;
    out.lat = gasx.lat;
elseif isfield(gasy,'lon')
    out.lon = gasy.lon;
    out.lat = gasy.lat;
end
if isfield(gasx,'eql')
    out.eql = gasx.eql;
elseif isfield(gasy,'eql')
    out.eql = gasy.eql;
end
if isfield(gasx,'spv')
    out.spv = gasx.spv;
elseif isfield(gasy,'spv')
    out.spv = gasy.spv;
end


relstruct_out = out;
%
end

