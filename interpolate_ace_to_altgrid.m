function [ tanstruct_out ] = interpolate_ace_to_altgrid( tanstruct_in, altitude_grid )
%A function to interpolate the variables of ACE an data strcture to a
%user-defined altitude grid. It is recommended to apply the quality flags
%before using this function, with 'aplly_ace_flags.m'. The interpolation
%will be off if the -999 values are kept and the flags should be applied
%anyway. The quality flags cannot be interpolated here.
%
% *INPUT*
%           tanstruct_in: STRUCTURE - contains the gas specific ACE data.
%           This structure can be created with 'read_ace_ncdata.m' or with
%           'read_ace_ncdata_for_mat.m'.
%            
%           altitude_grid: ARRAY - a vector or matrix containing the
%           altitude grid onto which you want to interpolate the ace
%           variables. This altitude should be in units of km.
%
% *OUTPUT*
%           tanstruct_out: STRUCTURE - output has almost the same
%           fields as the input, but with the data interpolated onto the
%           input altitude grid. The new altitude grid is included. The
%           quality flags are removed as they cannot be interpolated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NJR - 11/2017
% NJR - 02/2018

%% define some things
interptype = 'pchip'; %the spline interpolation method is giving bad results near the boundaries
req_points = 4; % the number of points required for an interpolation
gas = tanstruct_in;
zgrid = altitude_grid;
if isvector(zgrid)
   if ~iscolumn(zgrid)
      error('The interpolation grid should be an array or a column vector') 
   end
end
lorbit = length(gas.occultation);% this is the same as the length of p
zace = gas.altitude_km; % this is the same altitude grid for all ace measurements
lgrid = length(zgrid(:,1));
colgrid = length(zgrid(1,:));
if colgrid > 1 && colgrid ~= lorbit
    error('the number of interpolation grids does not match the number of occultations')
end

% the fields of the new data file that are the same
gasout.source_file = gas.source_file;
gasout.occultation = gas.occultation;
gasout.sr1ss0 = gas.sr1ss0;
gasout.beta_angle = gas.beta_angle;
gasout.date_mjd = gas.date_mjd;
gasout.gas = gas.gas;
% the fields of the new tanstruct file that will be interpolated
gasout.vmr = nan(lgrid,lorbit);
% size(gasout.vmr)
gasout.vmr_error = nan(lgrid,lorbit);
%these two also stay the same
gasout.lon_tangent = gas.lon_tangent;
gasout.lat_tangent = gas.lat_tangent;
%if there is a pressure field
if isfield(gas,'pressure_hPa')
    gasout.pressure_hPa = nan(lgrid,lorbit);
    if length(gas.pressure_hPa(1,:)) == 1 || length(gas.altitude_km(1,:)) > 1
        warning('It seems like this structure has already been interpolated');
    end
end
if isfield(gas,'eql')
   gasout.eql = nan(lgrid,lorbit); 
end
gasout.altitude_km = zgrid;

%% Interpolate the fields of the ace structure
fprintf('\ninterpolating the ACE data\n')
zgridi = zgrid(:,1); % set this initially outside the loop so it doesn't have to reset each time if zgrid is a vector
zacei = zace(:,1);
if isfield(gas,'pressure_hPa')
  pacei = gas.pressure_hPa(:,1); 
  logpacei = log(pacei);
end
for i = 1:lorbit
    %i
    if length(zgrid(1,:)) > 1 % If there are multiple grids to which to interpolate, i.e. zgrid is an array
       zgridi = zgrid(:,i); % this is only changed if zgrid is not a vector
    end
    if length(zace(1,:)) > 1 % If the structure already has multiple altitude grids for some reason i.e. zace is an array
       zacei = zace(:,i); % this is only changed if zace is not a vector
    end
    lgood = ~isnan(gas.vmr(:,i)); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
    if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
        vmri = gas.vmr(lgood,i); % reduce the size of the vector
        vmr_errori = gas.vmr_error(lgood,i);
        zacei_s = zacei(lgood); % reduce the size as well
        %interpolate the fields in log-pressure space
        gasout.vmr(:,i) = interp1(zacei_s,vmri,zgridi,interptype,nan);
        gasout.vmr_error(:,i) = interp1(zacei_s,vmr_errori,zgridi,interptype,nan);
    end
    %interpolating the pressure grid might get a bit messy...but it seems
    %ok for now
    if isfield(gas,'pressure_hPa')
        if length(gas.pressure_hPa(1,:)) > 1 % this should really be an array if you are interpolating in altitude space. You shouldn't be interpolating already interpolated data
            logpacei = log(gas.pressure_hPa(:,i)); % this is only changed if it is an array
        end
        lgood = ~isnan(pacei); % the indices of the P field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
            logpacei = logpacei(lgood); % reduce the size of the vector
            zacei_s = zacei(lgood); % reduce the size as well
            logpaceouti = interp1(zacei_s,logpacei,zgridi,interptype,nan);
            gasout.pressure_hPa(:,i) = exp(logpaceouti);
        end
    end
    if isfield(gas,'lon')
        lgood = ~isnan(gas.lon(:,i)); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
            loni = gas.lon(lgood,i); % reduce the size of the vector
            zacei_s = zacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.lon(:,i) = interp1(zacei_s,loni,zgridi,interptype,nan);
        end
        lgood = ~isnan(gas.lat(:,i)); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
            lati = gas.lat(lgood,i); % reduce the size of the vector
            zacei_s = zacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.lat(:,i) = interp1(zacei_s,lati,zgridi,interptype,nan);
        end
    end
    if isfield(gas,'eql')
        lgood = ~isnan(gas.eql(:,i)); % the indices of the vmr field that do not contain nans. nans will be placed where there is no data using 'apply_ace_flags.m'
        if sum(lgood) >= req_points % need at least 4 points to perform the interpolation
            eqli = gas.eql(lgood,i); % reduce the size of the vector
            zacei_s = zacei(lgood); % reduce the size as well
            %interpolate the fields in log-pressure space
            gasout.eql(:,i) = interp1(zacei_s,eqli,zgridi,interptype,nan);
        end
    end
end
fprintf('Done\n')
tanstruct_out = gasout;
end

